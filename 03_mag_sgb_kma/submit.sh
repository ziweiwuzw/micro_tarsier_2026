#!/usr/bin/env bash
# Submit a new full run with Pan_troglodytes MetaWRAP bins.
# The historical SemiBin2-based KMA/SGB output is never modified.
set -Eeuo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PIPELINE_SCRIPT_DIR="$HERE"

export PAN_SOURCE=metawrap
export SOURCE_ROOT="${SOURCE_ROOT:-/path/to/scratch/08_mag_pipeline/all_species_3binners_v2}"
export PIPELINE_ROOT="${PIPELINE_ROOT:-/path/to/scratch/08_mag_pipeline/all_species_3binners_v2_pan_metawrap_kma}"
export LOG_ROOT="${LOG_ROOT:-/path/to/project/08_vamb_semibin2/logs/mag_all_species_3binners_v2_pan_metawrap_kma}"
export ANN_ROOT="${ANN_ROOT:-/path/to/scratch/annotation_stack}"
export ANN_OUT="${ANN_OUT:-$ANN_ROOT/results/sgb_annotation_all_species_3binners_v2_pan_metawrap_kma}"
export ANNOT_SCRIPT_DIR="${ANNOT_SCRIPT_DIR:-$HERE/../04_sgb_annotation}"

# Derive every run-specific path from the new root, even if the submitter has
# variables from an earlier pipeline in the shell environment.
export READS_MANIFEST="$PIPELINE_ROOT/reads.tsv"
export CANDIDATE_DIR="$PIPELINE_ROOT/candidates/genomes"
export CANDIDATE_MANIFEST="$PIPELINE_ROOT/candidates/candidate_manifest.tsv"
export CHECKM2_OUT="$PIPELINE_ROOT/checkm2"
export HQ_DIR="$PIPELINE_ROOT/hq/genomes"
export HQ_MANIFEST="$PIPELINE_ROOT/hq/hq_manifest.tsv"
export DREP_OUT="$PIPELINE_ROOT/drep"
export SGB_ROOT="$PIPELINE_ROOT/sgb"
export SGB_REP_DIR="$SGB_ROOT/representatives"
export REFERENCE_ROOT="$PIPELINE_ROOT/sgb_reference"
export KMA_ROOT="$PIPELINE_ROOT/kma"
export GTDBTK_MASH_DB="$PIPELINE_ROOT/gtdbtk_cache/gtdb_r220.msh"

# shellcheck disable=SC1091
source "$HERE/common.sh"

mode="${1:-check}"
[[ "$mode" == check || "$mode" == submit ]] || die 'Usage: bash submit.sh check|submit'
[[ "$PIPELINE_ROOT" == /* && "$ANN_OUT" == /* ]] || die 'Output roots must be absolute paths'
[[ "$PIPELINE_ROOT" == *pan_metawrap_kma ]] || die "Expected a dedicated pan_metawrap_kma output root: $PIPELINE_ROOT"
[[ "$ANN_OUT" == *pan_metawrap_kma ]] || die "Expected a dedicated pan_metawrap_kma annotation root: $ANN_OUT"
[[ ! -e "$PIPELINE_ROOT" && ! -L "$PIPELINE_ROOT" ]] || die "New pipeline output already exists: $PIPELINE_ROOT"
[[ ! -e "$ANN_OUT" && ! -L "$ANN_OUT" ]] || die "New annotation output already exists: $ANN_OUT"
require_file "$PILOT_GROUPS_TSV"
require_file "$INPUT_READS_MANIFEST"
require_file "$HERE/collect_candidates.py"
require_file "$HERE/stage_job.slurm"
require_file "$HERE/kma_job.slurm"
require_file "$HERE/build_sgb_source_tables.py"
require_file "$ANNOT_SCRIPT_DIR/launch_annotation.slurm"
require_file "$ANNOT_SCRIPT_DIR/submit.sh"
grep -q 'PAN_SOURCE' "$HERE/collect_candidates.py" ||
    die 'Collector has not been updated to accept PAN_SOURCE=metawrap'
for job_script in "$ANNOT_SCRIPT_DIR"/*.slurm; do
    grep -Fq 'source "${ANNOT_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}/common.sh"' "$job_script" ||
        die "Outdated annotation Slurm script: $job_script"
done
[[ -x "$KMA_ENV/bin/kma" && -x "$KMA_ENV/bin/Rscript" ]] ||
    die "KMA/MicrobeCensus environment incomplete: $KMA_ENV"
reads_count="$(awk 'END {print NR-1}' "$INPUT_READS_MANIFEST")"
(( reads_count == 547 )) || die "Expected 547 reads manifest rows, found $reads_count"

groups=()
while IFS= read -r group; do
    groups+=("$group")
done < <(find "$SOURCE_ROOT/groups" -mindepth 1 -maxdepth 1 -type d -print | sort)
(( ${#groups[@]} == 32 )) || die "Expected 32 groups, found ${#groups[@]}"
pan_count=0
for group in "${groups[@]}"; do
    species="${group##*/}"
    refinement="$group/refinement"
    marker="$refinement/.metawrap_refinement.done"
    bins="$refinement/metawrap_50_10_bins"
    [[ -s "$marker" && -d "$bins" ]] || die "Incomplete MetaWRAP refinement: $species"
    bad_link="$(find -L "$bins" -maxdepth 1 -type l -name '*.fa' -print -quit)"
    [[ -z "$bad_link" ]] || die "Broken final-bin symlink: $bad_link"
    empty_bin="$(find -L "$bins" -maxdepth 1 -type f -name '*.fa' -size 0 -print -quit)"
    [[ -z "$empty_bin" ]] || die "Empty final-bin FASTA: $empty_bin"
    count="$(find -L "$bins" -maxdepth 1 -type f -name '*.fa' -size +0c | wc -l | tr -d '[:space:]')"
    (( count > 0 )) || die "No nonempty MetaWRAP FASTA in $bins"
    if [[ "$species" == Pan_troglodytes ]]; then pan_count="$count"; fi
done
(( pan_count > 0 )) || die 'Pan_troglodytes group is missing'
printf '[READY] groups=32 Pan_MetaWRAP_bins=%s\nNEW_PIPELINE_ROOT=%s\nNEW_ANN_OUT=%s\n' \
    "$pan_count" "$PIPELINE_ROOT" "$ANN_OUT"
[[ "$mode" == submit ]] || exit 0

command -v sbatch >/dev/null || die 'sbatch unavailable; submit on Spartan'
max="${MAX_CONCURRENT:-8}"
[[ "$max" =~ ^[1-9][0-9]*$ ]] || die 'MAX_CONCURRENT must be positive'
mkdir -p "$LOG_ROOT"
receipt="$LOG_ROOT/submission_$(date +%Y%m%d_%H%M%S).tsv"
printf 'stage\tjob_id\n' > "$receipt"
submit() {
    local stage="$1" response id
    shift
    response="$(sbatch --parsable --export=ALL --output="$LOG_ROOT/%x_%j.out" \
        --error="$LOG_ROOT/%x_%j.err" "$@")"
    id="${response%%;*}"
    [[ "$id" =~ ^[1-9][0-9]*$ ]] || die "Unexpected sbatch response for $stage: $response"
    printf '%s\t%s\n' "$stage" "$id" >> "$receipt"
    printf '%s\n' "$id"
}
submit_array() {
    local stage="$1" response id
    shift
    response="$(sbatch --parsable --export=ALL --output="$LOG_ROOT/%x_%A_%a.out" \
        --error="$LOG_ROOT/%x_%A_%a.err" "$@")"
    id="${response%%;*}"
    [[ "$id" =~ ^[1-9][0-9]*$ ]] || die "Unexpected sbatch response for $stage: $response"
    printf '%s\t%s\n' "$stage" "$id" >> "$receipt"
    printf '%s\n' "$id"
}

collect="$(submit collect --job-name=mag_collect --cpus-per-task=1 --mem=8G --time=14-00:00:00 "$HERE/stage_job.slurm" collect)"
checkm2="$(submit checkm2 --dependency="afterok:$collect" "$HERE/run_checkm2.slurm")"
filter="$(submit filter --dependency="afterok:$checkm2" --job-name=mag_filter --cpus-per-task=1 --mem=8G --time=14-00:00:00 "$HERE/stage_job.slurm" filter)"
drep="$(submit drep --dependency="afterok:$filter" "$HERE/run_drep.slurm")"
sgb="$(submit sgb --dependency="afterok:$drep" --job-name=mag_sgb_ids --cpus-per-task=1 --mem=8G --time=14-00:00:00 "$HERE/stage_job.slurm" sgb)"
source_metadata="$(submit source_metadata --dependency="afterok:$sgb" --job-name=sgb_sources --cpus-per-task=1 --mem=8G --time=1-00:00:00 "$HERE/stage_job.slurm" source_metadata)"
taxonomy="$(submit taxonomy --dependency="afterok:$sgb" "$HERE/run_gtdbtk_taxonomy.slurm")"
tree="$(submit tree --dependency="afterok:$taxonomy" "$HERE/infer_gtdbtk_tree.slurm" "$PIPELINE_ROOT")"
reference="$(submit reference --dependency="afterok:$sgb" --job-name=sgb_reference --cpus-per-task=1 --mem=32G --time=1-00:00:00 "$HERE/stage_job.slurm" reference)"
index="$(submit kma_index --dependency="afterok:$reference" --job-name=kma_index --cpus-per-task=8 --mem=256G --time=2-00:00:00 "$HERE/kma_job.slurm" index)"
samples="$(submit_array kma_samples --dependency="afterok:$index" --job-name=kma_sample --cpus-per-task=16 --mem=64G --time=2-00:00:00 --array="1-547%$max" "$HERE/kma_job.slurm" sample)"
merge="$(submit kma_merge --dependency="afterok:$samples" --job-name=kma_merge --cpus-per-task=1 --mem=16G --time=1-00:00:00 "$HERE/kma_job.slurm" merge)"
annotation="$(submit functional_annotation_launcher --dependency="afterok:$sgb" "$ANNOT_SCRIPT_DIR/launch_annotation.slurm")"

cat "$receipt"
printf 'Receipt: %s\n' "$receipt"
