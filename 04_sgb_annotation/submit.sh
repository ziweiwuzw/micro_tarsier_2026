#!/usr/bin/env bash
set -Eeuo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
export ANNOT_SCRIPT_DIR="$SCRIPT_DIR"
command -v sbatch >/dev/null || die 'sbatch unavailable; submit on Spartan'
for job_script in "$SCRIPT_DIR"/*.slurm; do
    grep -Fq 'source "${ANNOT_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}/common.sh"' "$job_script" ||
        die "Old Slurm script still present: $job_script; upload the updated .slurm files before submitting"
done
if [[ "${1:-}" == fegenie ]]; then
    require_file "$ANN_OUT/input_manifest.tsv"
    require_dir "$ANN_OUT/01_orf"
    for file in install_fegenie.slurm test_fegenie.slurm run_fegenie_batch.slurm summary_job.slurm; do
        require_file "$SCRIPT_DIR/$file"
    done
    mag_count=$(awk 'END {print NR-1}' "$ANN_OUT/input_manifest.tsv")
    (( mag_count >= 2 )) || die 'FeGenie pilot requires at least two MAGs'
    while IFS=$'\t' read -r mag sgb rep genome genome_sha256; do
        [[ "$mag" == mag_id ]] && continue
        require_file "$ANN_OUT/01_orf/$mag/.done"
        require_file "$ANN_OUT/01_orf/$mag/$mag.fa"
    done < "$ANN_OUT/input_manifest.tsv"
    stride="${FEGENIE_BATCH_STRIDE:-500}"
    [[ "$stride" =~ ^[1-9][0-9]*$ ]] || die 'FEGENIE_BATCH_STRIDE must be positive'
    array_count=$(( mag_count < stride ? mag_count : stride ))
    mkdir -p "$ANN_OUT/logs" "$ANN_OUT/diagnostics"
    submit() {
        local response
        response=$(sbatch --parsable --output="$ANN_OUT/logs/%x-%A_%a.out" "$@")
        printf '%s\n' "${response%%;*}"
    }
    if [[ -s "$ANN_ROOT/software/fegenie_install/.done" ]]; then
        install=REUSE
        test=$(submit "$SCRIPT_DIR/test_fegenie.slurm")
    else
        install=$(submit "$SCRIPT_DIR/install_fegenie.slurm")
        test=$(submit --dependency="afterok:$install" "$SCRIPT_DIR/test_fegenie.slurm")
    fi
    batch=$(submit --dependency="afterok:$test" --array="1-${array_count}" \
        --export="ALL,FEGENIE_BATCH_STRIDE=$stride" "$SCRIPT_DIR/run_fegenie_batch.slurm")
    tables=$(submit --dependency="afterok:$batch" --job-name=fegenie_tables --partition=sapphire \
        --cpus-per-task=1 --mem=12G --time=12:00:00 "$SCRIPT_DIR/summary_job.slurm" fegenie_export)
    printf 'stage\tjob_id\ninstall\t%s\npilot\t%s\nbatch\t%s\ntables\t%s\n' \
        "$install" "$test" "$batch" "$tables" | tee "$ANN_OUT/diagnostics/fegenie_jobs.tsv"
    exit 0
fi
if [[ -z "${SGB_LIST:-}" && -s "$ANN_OUT/selected_sgbs.tsv" ]]; then
    SGB_LIST="$ANN_OUT/selected_sgbs.tsv"
fi
if [[ -n "${SGB_LIST:-}" ]]; then
    python3 "$SCRIPT_DIR/annotation.py" --pipeline-root "$PIPELINE_ROOT" --out "$ANN_OUT" prepare --sgb-list "$SGB_LIST"
else
    python3 "$SCRIPT_DIR/annotation.py" --pipeline-root "$PIPELINE_ROOT" --out "$ANN_OUT" prepare
fi
n=$(awk 'END {print NR-1}' "$ANN_OUT/input_manifest.tsv")
(( n > 0 )) || die 'Empty MAG manifest'
log_dir="$ANN_OUT/logs"
mkdir -p "$log_dir"
submit() {
    local response
    response=$(sbatch --parsable --export=ALL --output="$log_dir/%x-%j.out" "$@") ||
        die 'sbatch failed; no further annotation jobs were submitted'
    printf '%s\n' "${response%%;*}"
}
submit_array_chunks() {
    local total="$1" concurrency="$2" script="$3" offset=0 size end job
    shift 3
    local extra=("$@")
    local chunk_limit="${SLURM_ARRAY_CHUNK_SIZE:-500}"
    [[ "$total" =~ ^[1-9][0-9]*$ ]] || die "Invalid array total: $total"
    [[ "$concurrency" =~ ^[1-9][0-9]*$ ]] || die "Invalid array concurrency: $concurrency"
    [[ "$chunk_limit" =~ ^[1-9][0-9]*$ ]] || die "SLURM_ARRAY_CHUNK_SIZE must be positive"
    local jobs=()
    local previous=""
    while (( offset < total )); do
        size=$(( total - offset < chunk_limit ? total - offset : chunk_limit ))
        end="$size"
        local chunk_args=()
        if [[ -n "$previous" ]]; then
            chunk_args=("--dependency=afterok:$previous")
        elif ((${#extra[@]})); then
            chunk_args=("${extra[@]}")
        fi
        if ((${#chunk_args[@]})); then
            job=$(sbatch --parsable --export="ALL,ARRAY_OFFSET=$offset" \
                --output="$log_dir/%x-%j.out" --array="1-${end}%${concurrency}" \
                "${chunk_args[@]}" "$script") ||
                die "sbatch failed at array offset $offset; no later chunks were submitted"
        else
            job=$(sbatch --parsable --export="ALL,ARRAY_OFFSET=$offset" \
                --output="$log_dir/%x-%j.out" --array="1-${end}%${concurrency}" \
                "$script") ||
                die "sbatch failed at array offset $offset; no later chunks were submitted"
        fi
        job="${job%%;*}"
        [[ "$job" =~ ^[1-9][0-9]*$ ]] || die "Unexpected sbatch job ID: $job"
        jobs+=("$job")
        previous="$job"
        offset=$((offset + size))
    done
    local joined
    joined="$(IFS=:; printf '%s' "${jobs[*]}")"
    printf '%s\n' "$joined"
}
case "${1:-}" in
    genes)
        submit_array_chunks "$n" "${MAX_CONCURRENT:-10}" "$SCRIPT_DIR/run_prodigal.slurm" ;;
    dram2)
        submit "$SCRIPT_DIR/run_dram2.slurm" ;;
    dbcan)
        submit_array_chunks "$n" "${MAX_CONCURRENT:-10}" "$SCRIPT_DIR/run_dbcan.slurm" ;;
    homology-index)
        submit "$SCRIPT_DIR/build_homology_db.slurm" ;;
    homology)
        submit_array_chunks "$n" "${MAX_HOMOLOGY_CONCURRENT:-3}" "$SCRIPT_DIR/screen_homology.slurm" ;;
    tables)
        require_file "$ANN_OUT/03_dram2_all_members/.done"
        while IFS=$'\t' read -r mag sgb rep genome genome_sha256; do
            [[ "$mag" == mag_id ]] && continue
            require_file "$ANN_OUT/04_dbcan/$mag/.done"
            require_file "$ANN_OUT/06_homology/$mag/.done"
        done < "$ANN_OUT/input_manifest.tsv"
        submit "$SCRIPT_DIR/export_trait_tables.slurm" ;;
    summarize)
        submit --job-name=sgb_traits --cpus-per-task=1 --mem=16G --time=08:00:00 "$SCRIPT_DIR/summary_job.slurm" summarize ;;
    integrate)
        submit --job-name=sgb_integrate --cpus-per-task=1 --mem=32G --time=08:00:00 "$SCRIPT_DIR/summary_job.slurm" integrate ;;
    modules|all)
        genes_job=$(submit_array_chunks "$n" "${MAX_CONCURRENT:-10}" "$SCRIPT_DIR/run_prodigal.slurm")
        dram_job=$(submit --dependency="afterok:$genes_job" "$SCRIPT_DIR/run_dram2.slurm")
        dbcan_job=$(submit_array_chunks "$n" "${MAX_CONCURRENT:-10}" "$SCRIPT_DIR/run_dbcan.slurm" --dependency="afterok:$genes_job")
        ref_job=$(submit "$SCRIPT_DIR/build_homology_db.slurm")
        homology_job=$(submit_array_chunks "$n" "${MAX_HOMOLOGY_CONCURRENT:-3}" "$SCRIPT_DIR/screen_homology.slurm" --dependency="afterok:$genes_job:$ref_job")
        printf 'genes=%s DRAM2=%s dbCAN=%s reference=%s homology=%s\n' \
            "$genes_job" "$dram_job" "$dbcan_job" "$ref_job" "$homology_job"
        if [[ "$1" == all ]]; then
            tables_job=$(submit --dependency="afterok:$dram_job:$dbcan_job:$homology_job" \
                "$SCRIPT_DIR/export_trait_tables.slurm")
            printf 'source_tables=%s\n' "$tables_job"
        fi ;;
    *) die 'Usage: bash submit.sh genes|dram2|dbcan|homology-index|homology|tables|summarize|integrate|modules|all|fegenie' ;;
esac
