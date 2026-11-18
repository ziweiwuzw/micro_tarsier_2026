#!/usr/bin/env bash
SCRIPT_VERSION="2026-07-23-auto-workroot-v2"
#SBATCH --job-name=megahit_species
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=7-00:00:00
#SBATCH --output=megahit_logs/%x_%A_%a.out
#SBATCH --error=megahit_logs/%x_%A_%a.err

# ==============================================================================
# MEGAHIT per-sample assembly with per-species Slurm array tasks
# ==============================================================================
# Each Slurm array task handles one species. Samples within that species are
# assembled one-by-one with MEGAHIT.
#
# Workflow for each species:
#   1. Copy one sample's paired nonhost FASTQ files to scratch.
#   2. Assemble that sample independently with MEGAHIT in scratch.
#   3. Validate final.contigs.fa and remove copied FASTQ/intermediate contigs.
#   4. Repeat until every sample from the species has a valid assembly.
#   5. Delete that species' *_1/2.nonhost.fastq.gz files from 02_host_removed_reads.
#   6. Transfer each sample assembly to the project directory using flat names:
#        <sample>_final.contigs.fa
#        <sample>_megahit.log
#        <sample>_options.json
#        <sample>_megahit.done
#   7. Validate all final outputs and only then remove the species scratch folder.
#
# Automatic project detection:
#   - WORK_ROOT defaults to the physical current directory: pwd -P
#   - RUN_NAME is generated from the current path, for example:
#       /path/to/metadata
#       -> the project-specific MEGAHIT output directory
#   - Different project directories therefore use different scratch subfolders.
#
# Important:
#   - Run manifest/submit from the directory containing 02_host_removed_reads.
#   - DELETE_NONHOST_READS=1 deletes only the paired nonhost FASTQ files.
#     host_mapping.flagstat.txt and host_removed.done are retained.
#   - If reads have already been deleted and transfer is interrupted, do not
#     delete scratch or rebuild manifests. Run this script with resume.
#   - Do not submit two copies of this pipeline concurrently from the same
#     WORK_ROOT/RUN_NAME.
# ============================================================================== 

set -Eeuo pipefail
shopt -s nullglob

SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"

log() {
    printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2
}

die() {
    log "ERROR: $*"
    exit 1
}

need() {
    command -v "$1" >/dev/null 2>&1 || die "Command not found: $1"
}

safe_name() {
    printf '%s' "$1" | sed 's#[^A-Za-z0-9._-]#_#g'
}

# ----------------------- Automatic project identification --------------------
# WORK_ROOT can still be explicitly overridden, but otherwise always uses the
# directory from which the user launches manifest/submit/resume.
WORK_ROOT="${WORK_ROOT:-$(pwd -P)}"
WORK_ROOT="$(readlink -f "$WORK_ROOT")"
WORK_ROOT="${WORK_ROOT%/}"

PROJECT_DIR_NAME="$(safe_name "$(basename "$WORK_ROOT")")"
PROJECT_ID="$({
    printf '%s\n' "$WORK_ROOT" | awk -F '/' '
        {
            for (i = 1; i <= NF; i++) {
                if ($i ~ /^punim[0-9]+$/) {
                    id = $i
                }
            }
        }
        END {print id}
    '
} || true)"

if [[ -z "$PROJECT_ID" ]]; then
    PROJECT_ID="$(safe_name "$(basename "$(dirname "$WORK_ROOT")")")"
else
    PROJECT_ID="$(safe_name "$PROJECT_ID")"
fi

AUTO_RUN_NAME="${PROJECT_ID}_${PROJECT_DIR_NAME}_megahit"
RUN_NAME="$(safe_name "${RUN_NAME:-$AUTO_RUN_NAME}")"

# ------------------------------- Main paths -----------------------------------
READS_ROOT="${READS_ROOT:-$WORK_ROOT/02_host_removed_reads/reads}"
OUT_ROOT="${OUT_ROOT:-$WORK_ROOT/03_megahit_sample_assemblies}"
ASSEMBLY_ROOT="${ASSEMBLY_ROOT:-$OUT_ROOT/assemblies}"
REPORT_DIR="${REPORT_DIR:-$OUT_ROOT/reports}"
LOG_DIR="${LOG_DIR:-$WORK_ROOT/megahit_logs}"

SCRATCH_ROOT="${SCRATCH_ROOT:-/path/to/scratch/megahit_sample_assembly_work}"
RUN_SCRATCH="${RUN_SCRATCH:-$SCRATCH_ROOT/$RUN_NAME}"
RUN_IDENTITY_FILE="$RUN_SCRATCH/.megahit_run_identity.tsv"

# ------------------------------ Slurm resources -------------------------------
SBATCH_ACCOUNT="${SBATCH_ACCOUNT:?Set SBATCH_ACCOUNT for your Slurm account}"
SBATCH_PARTITION="${SBATCH_PARTITION:-sapphire}"
MEGAHIT_CPUS="${MEGAHIT_CPUS:-32}"
MEGAHIT_MEM="${MEGAHIT_MEM:-128G}"
MEGAHIT_TIME="${MEGAHIT_TIME:-7-00:00:00}"
MAX_SPECIES_CONCURRENT="${MAX_SPECIES_CONCURRENT:-2}"

# ------------------------------ MEGAHIT options -------------------------------
MEGAHIT_THREADS="${MEGAHIT_THREADS:-$MEGAHIT_CPUS}"
MIN_CONTIG_LEN="${MIN_CONTIG_LEN:-1000}"
MEGAHIT_EXTRA_ARGS="${MEGAHIT_EXTRA_ARGS:-}"

# ------------------------------- Safety options -------------------------------
DELETE_NONHOST_READS="${DELETE_NONHOST_READS:-1}"
VERIFY_INPUT_GZIP="${VERIFY_INPUT_GZIP:-1}"
OVERWRITE="${OVERWRITE:-0}"
PRUNE_INTERMEDIATE="${PRUNE_INTERMEDIATE:-1}"

# --------------------------------- Reports ------------------------------------
SAMPLE_MANIFEST="$REPORT_DIR/megahit_sample_manifest.tsv"
SPECIES_MANIFEST="$REPORT_DIR/megahit_species_manifest.tsv"
PAIRING_ISSUES="$REPORT_DIR/megahit_pairing_issues.tsv"
MANIFEST_SUMMARY="$REPORT_DIR/megahit_manifest_summary.txt"
SAMPLE_STATUS_DIR="$REPORT_DIR/sample_status"
SPECIES_STATUS_DIR="$REPORT_DIR/species_status"
SAMPLE_SUMMARY="$REPORT_DIR/megahit_sample_summary.tsv"
SPECIES_SUMMARY="$REPORT_DIR/megahit_species_summary.tsv"
RUN_SUMMARY="$REPORT_DIR/megahit_run_summary.txt"

mkdir -p "$REPORT_DIR" "$SAMPLE_STATUS_DIR" "$SPECIES_STATUS_DIR"

init_module() {
    if ! command -v module >/dev/null 2>&1; then
        [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
        [[ -f /usr/share/Modules/init/bash ]] && source /usr/share/Modules/init/bash
    fi
    command -v module >/dev/null 2>&1 || die "Environment Modules command is unavailable"
}

load_megahit() {
    init_module
    module load MEGAHIT
    need megahit
}

ensure_run_identity() {
    local existing_root existing_run tmp

    mkdir -p "$RUN_SCRATCH"

    if [[ -s "$RUN_IDENTITY_FILE" ]]; then
        existing_root="$(awk -F '\t' '$1 == "WORK_ROOT" {print $2; exit}' "$RUN_IDENTITY_FILE")"
        existing_run="$(awk -F '\t' '$1 == "RUN_NAME" {print $2; exit}' "$RUN_IDENTITY_FILE")"

        [[ "$existing_root" == "$WORK_ROOT" ]] || die \
            "Scratch collision: $RUN_SCRATCH belongs to WORK_ROOT=$existing_root, not $WORK_ROOT"
        [[ "$existing_run" == "$RUN_NAME" ]] || die \
            "Scratch collision: stored RUN_NAME=$existing_run, current RUN_NAME=$RUN_NAME"
        return 0
    fi

    tmp="${RUN_IDENTITY_FILE}.tmp.$$"
    {
        printf 'WORK_ROOT\t%s\n' "$WORK_ROOT"
        printf 'RUN_NAME\t%s\n' "$RUN_NAME"
        printf 'SCRIPT_PATH\t%s\n' "$SCRIPT_PATH"
        printf 'CREATED\t%s\n' "$(date -Is)"
    } > "$tmp"
    mv "$tmp" "$RUN_IDENTITY_FILE"
}

manifest_row() {
    local file="$1"
    local task_id="$2"
    awk -F '\t' -v row="$((task_id + 1))" 'NR == row {print; exit}' "$file"
}

contig_file_valid() {
    local file="$1"
    [[ -s "$file" ]] || return 1
    grep -q '^>' "$file" || return 1
}

scratch_assembly_valid() {
    local dir="$1"
    contig_file_valid "$dir/final.contigs.fa"
}

write_scratch_done() {
    local dir="$1"
    local done="$dir/megahit.done"

    scratch_assembly_valid "$dir" || return 1
    if [[ ! -s "$done" ]]; then
        {
            printf 'completed=%s\n' "$(date -Is)"
            printf 'final_contigs=%s\n' "$dir/final.contigs.fa"
            printf 'contigs=%s\n' "$(grep -c '^>' "$dir/final.contigs.fa")"
            printf 'bytes=%s\n' "$(stat -c '%s' "$dir/final.contigs.fa")"
        } > "$done"
    fi
}

write_final_done() {
    local final_contigs="$1"
    local final_done="$2"
    local source_note="${3:-unknown}"
    local tmp

    contig_file_valid "$final_contigs" || return 1
    mkdir -p "$(dirname "$final_done")"
    tmp="${final_done}.tmp.${SLURM_JOB_ID:-$$}"
    {
        printf 'completed=%s\n' "$(date -Is)"
        printf 'final_contigs=%s\n' "$final_contigs"
        printf 'contigs=%s\n' "$(grep -c '^>' "$final_contigs")"
        printf 'bytes=%s\n' "$(stat -c '%s' "$final_contigs")"
        printf 'source=%s\n' "$source_note"
    } > "$tmp"
    mv -f "$tmp" "$final_done"
}

final_output_valid() {
    local final_contigs="$1"
    local final_done="$2"

    contig_file_valid "$final_contigs" || return 1

    # A contig file is transferred with an atomic rename. If the marker was lost
    # after that rename, safely reconstruct it instead of rebuilding the sample.
    if [[ ! -s "$final_done" ]]; then
        write_final_done "$final_contigs" "$final_done" "recovered_existing_final_contigs"
    fi

    [[ -s "$final_done" ]]
}

prune_scratch_assembly() {
    local dir="$1"
    [[ "$PRUNE_INTERMEDIATE" == 1 ]] || return 0
    scratch_assembly_valid "$dir" || return 1

    rm -rf "$dir/intermediate_contigs"
    rm -f "$dir"/*.tmp "$dir"/*.tmp.* 2>/dev/null || true
}

copy_input_atomic() {
    local src="$1"
    local dst="$2"
    local tmp

    [[ -s "$src" ]] || die "Missing input: $src"
    mkdir -p "$(dirname "$dst")"

    tmp="${dst}.copying.${SLURM_JOB_ID:-$$}"
    rm -f "$tmp"
    rsync -a --partial "$src" "$tmp"

    [[ -s "$tmp" ]] || die "Copy failed: $src"
    [[ "$(stat -c '%s' "$src")" -eq "$(stat -c '%s' "$tmp")" ]] || \
        die "Input size mismatch: $src"

    if [[ "$VERIFY_INPUT_GZIP" == 1 ]]; then
        gzip -t "$tmp"
    fi

    mv -f "$tmp" "$dst"
}

copy_scratch_to_flat_output() {
    local scratch_dir="$1"
    local sample="$2"
    local final_contigs="$3"
    local final_log="$4"
    local final_options="$5"
    local final_done="$6"

    local final_dir stage_dir stage_contigs stage_log stage_options stage_done

    scratch_assembly_valid "$scratch_dir" || die "Invalid scratch assembly: $scratch_dir"

    final_dir="$(dirname "$final_contigs")"
    stage_dir="$final_dir/.${sample}.transfer.${SLURM_JOB_ID:-$$}"
    stage_contigs="$stage_dir/${sample}_final.contigs.fa"
    stage_log="$stage_dir/${sample}_megahit.log"
    stage_options="$stage_dir/${sample}_options.json"
    stage_done="$stage_dir/${sample}_megahit.done"

    mkdir -p "$final_dir"
    rm -rf "$stage_dir"
    mkdir -p "$stage_dir"

    rsync -a "$scratch_dir/final.contigs.fa" "$stage_contigs"
    contig_file_valid "$stage_contigs" || die "Transferred contig validation failed: $stage_contigs"

    if [[ -f "$scratch_dir/log" ]]; then
        rsync -a "$scratch_dir/log" "$stage_log"
    fi
    if [[ -f "$scratch_dir/options.json" ]]; then
        rsync -a "$scratch_dir/options.json" "$stage_options"
    fi

    {
        printf 'completed=%s\n' "$(date -Is)"
        printf 'final_contigs=%s\n' "$final_contigs"
        printf 'contigs=%s\n' "$(grep -c '^>' "$stage_contigs")"
        printf 'bytes=%s\n' "$(stat -c '%s' "$stage_contigs")"
        printf 'source=%s\n' "$scratch_dir"
    } > "$stage_done"

    if [[ "$OVERWRITE" == 1 ]]; then
        rm -f "$final_contigs" "$final_log" "$final_options" "$final_done"
    fi

    # Files are renamed on the same final filesystem. The done marker is moved
    # last, so it records completion of the flat-output transfer.
    mv -f "$stage_contigs" "$final_contigs"
    [[ -f "$stage_log" ]] && mv -f "$stage_log" "$final_log"
    [[ -f "$stage_options" ]] && mv -f "$stage_options" "$final_options"
    mv -f "$stage_done" "$final_done"
    rmdir "$stage_dir" 2>/dev/null || rm -rf "$stage_dir"

    final_output_valid "$final_contigs" "$final_done" || \
        die "Final flat output validation failed: $final_contigs"
}

build_manifests() {
    need find
    need awk
    need sort
    need dirname
    need basename
    need readlink

    [[ -d "$READS_ROOT" ]] || {
        printf 'Current WORK_ROOT: %s\n' "$WORK_ROOT" >&2
        printf 'Expected input:    %s\n' "$READS_ROOT" >&2
        die "Input directory not found. Run the script from the project directory containing 02_host_removed_reads"
    }

    mkdir -p "$REPORT_DIR" "$SAMPLE_STATUS_DIR" "$SPECIES_STATUS_DIR"

    printf 'species\tproject\tsample_id\tr1\texpected_r2\tissue\n' > "$PAIRING_ISSUES"
    printf 'sample_id\tspecies\tproject\tr1\tr2\tscratch_assembly_dir\tfinal_contigs\tfinal_log\tfinal_options\tfinal_done\n' > "$SAMPLE_MANIFEST"

    local species_dir species r1 r2 rel project sample
    local scratch_dir final_prefix final_contigs final_log final_options final_done

    while IFS= read -r -d '' species_dir; do
        species="$(basename "$species_dir")"

        while IFS= read -r -d '' r1; do
            r2="${r1%_1.nonhost.fastq.gz}_2.nonhost.fastq.gz"
            rel="${r1#${species_dir}/}"
            project="$(dirname "$rel")"
            [[ "$project" != '.' ]] || project="no_project"
            sample="$(basename "${r1%_1.nonhost.fastq.gz}")"

            if [[ ! -s "$r2" ]]; then
                printf '%s\t%s\t%s\t%s\t%s\tmissing_R2\n' \
                    "$species" "$project" "$sample" "$r1" "$r2" >> "$PAIRING_ISSUES"
                continue
            fi

            scratch_dir="$RUN_SCRATCH/$species/results/$project/$sample"
            final_prefix="$ASSEMBLY_ROOT/$species/$project/$sample"
            final_contigs="${final_prefix}_final.contigs.fa"
            final_log="${final_prefix}_megahit.log"
            final_options="${final_prefix}_options.json"
            final_done="${final_prefix}_megahit.done"

            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "$sample" "$species" "$project" "$r1" "$r2" \
                "$scratch_dir" "$final_contigs" "$final_log" \
                "$final_options" "$final_done" >> "$SAMPLE_MANIFEST"
        done < <(
            find "$species_dir" -type f -name '*_1.nonhost.fastq.gz' -print0 | LC_ALL=C sort -z
        )

        while IFS= read -r -d '' r2; do
            r1="${r2%_2.nonhost.fastq.gz}_1.nonhost.fastq.gz"
            [[ -s "$r1" ]] && continue

            rel="${r2#${species_dir}/}"
            project="$(dirname "$rel")"
            [[ "$project" != '.' ]] || project="no_project"
            sample="$(basename "${r2%_2.nonhost.fastq.gz}")"

            printf '%s\t%s\t%s\t%s\t%s\tmissing_R1\n' \
                "$species" "$project" "$sample" "$r1" "$r2" >> "$PAIRING_ISSUES"
        done < <(
            find "$species_dir" -type f -name '*_2.nonhost.fastq.gz' -print0 | LC_ALL=C sort -z
        )
    done < <(
        find "$READS_ROOT" -mindepth 1 -maxdepth 1 -type d -print0 | LC_ALL=C sort -z
    )

    printf 'species_task_id\tspecies\texpected_samples\tpairing_issues\tsource_species_dir\tscratch_species_dir\tfinal_species_dir\n' > "$SPECIES_MANIFEST"

    local task_id=0 expected issues
    while IFS= read -r species; do
        [[ -n "$species" ]] || continue

        expected="$(awk -F '\t' -v sp="$species" 'NR > 1 && $2 == sp {n++} END {print n+0}' "$SAMPLE_MANIFEST")"
        issues="$(awk -F '\t' -v sp="$species" 'NR > 1 && $1 == sp {n++} END {print n+0}' "$PAIRING_ISSUES")"

        # Empty species folders and species with pairing issues are recorded in
        # reports but excluded from the runnable Slurm array.
        if (( expected == 0 || issues > 0 )); then
            continue
        fi

        ((task_id += 1))
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$task_id" "$species" "$expected" "$issues" \
            "$READS_ROOT/$species" "$RUN_SCRATCH/$species" \
            "$ASSEMBLY_ROOT/$species" >> "$SPECIES_MANIFEST"
    done < <(
        {
            awk -F '\t' 'NR > 1 {print $2}' "$SAMPLE_MANIFEST"
            awk -F '\t' 'NR > 1 {print $1}' "$PAIRING_ISSUES"
        } | LC_ALL=C sort -u
    )

    local sample_count issue_count
    sample_count=$(( $(wc -l < "$SAMPLE_MANIFEST") - 1 ))
    issue_count=$(( $(wc -l < "$PAIRING_ISSUES") - 1 ))

    {
        printf 'SCRIPT_VERSION=%s\n' "$SCRIPT_VERSION"
        printf 'WORK_ROOT=%s\n' "$WORK_ROOT"
        printf 'PROJECT_ID=%s\n' "$PROJECT_ID"
        printf 'RUN_NAME=%s\n' "$RUN_NAME"
        printf 'READS_ROOT=%s\n' "$READS_ROOT"
        printf 'RUN_SCRATCH=%s\n' "$RUN_SCRATCH"
        printf 'OUT_ROOT=%s\n' "$OUT_ROOT"
        printf 'ASSEMBLY_ROOT=%s\n\n' "$ASSEMBLY_ROOT"
        printf 'species_tasks=%s\n' "$task_id"
        printf 'paired_samples=%s\n' "$sample_count"
        printf 'pairing_issues=%s\n' "$issue_count"
        printf 'MAX_SPECIES_CONCURRENT=%s\n' "$MAX_SPECIES_CONCURRENT"
        printf 'MEGAHIT_THREADS=%s\n' "$MEGAHIT_THREADS"
        printf 'MIN_CONTIG_LEN=%s\n' "$MIN_CONTIG_LEN"
        printf 'DELETE_NONHOST_READS=%s\n' "$DELETE_NONHOST_READS"
    } > "$MANIFEST_SUMMARY"

    cat "$MANIFEST_SUMMARY"

    if (( issue_count > 0 )); then
        log "Pairing issues detected. Affected species were excluded; inspect: $PAIRING_ISSUES"
    fi
}

run_megahit_sample() {
    local sample="$1"
    local species="$2"
    local project="$3"
    local r1="$4"
    local r2="$5"
    local scratch_assembly="$6"

    local sample_input="$RUN_SCRATCH/$species/inputs/$project/$sample"
    local sr1="$sample_input/$(basename "$r1")"
    local sr2="$sample_input/$(basename "$r2")"
    local start elapsed
    local -a extra_args=()

    start="$(date +%s)"

    if [[ "$OVERWRITE" != 1 ]] && scratch_assembly_valid "$scratch_assembly"; then
        write_scratch_done "$scratch_assembly"
        prune_scratch_assembly "$scratch_assembly"
        rm -rf "$sample_input"
        log "Reuse complete scratch assembly: $species/$project/$sample"
        return 0
    fi

    [[ -s "$r1" && -s "$r2" ]] || \
        die "Source reads missing and assembly incomplete: $species/$project/$sample"

    mkdir -p "$sample_input" "$(dirname "$scratch_assembly")"
    copy_input_atomic "$r1" "$sr1"
    copy_input_atomic "$r2" "$sr2"

    if [[ -n "$MEGAHIT_EXTRA_ARGS" ]]; then
        # MEGAHIT_EXTRA_ARGS is intentionally interpreted as shell-style words.
        # shellcheck disable=SC2206
        extra_args=( $MEGAHIT_EXTRA_ARGS )
    fi

    if [[ -d "$scratch_assembly" && ! -s "$scratch_assembly/final.contigs.fa" ]]; then
        if [[ -s "$scratch_assembly/checkpoints.txt" || -d "$scratch_assembly/intermediate_contigs" ]]; then
            log "Continue interrupted MEGAHIT assembly: $species/$project/$sample"
            if ! megahit --continue -o "$scratch_assembly"; then
                log "MEGAHIT --continue failed; restarting sample from the beginning"
                rm -rf "$scratch_assembly"
                megahit \
                    -1 "$sr1" \
                    -2 "$sr2" \
                    --num-cpu-threads "$MEGAHIT_THREADS" \
                    --min-contig-len "$MIN_CONTIG_LEN" \
                    "${extra_args[@]}" \
                    -o "$scratch_assembly"
            fi
        else
            rm -rf "$scratch_assembly"
            megahit \
                -1 "$sr1" \
                -2 "$sr2" \
                --num-cpu-threads "$MEGAHIT_THREADS" \
                --min-contig-len "$MIN_CONTIG_LEN" \
                "${extra_args[@]}" \
                -o "$scratch_assembly"
        fi
    else
        rm -rf "$scratch_assembly"
        megahit \
            -1 "$sr1" \
            -2 "$sr2" \
            --num-cpu-threads "$MEGAHIT_THREADS" \
            --min-contig-len "$MIN_CONTIG_LEN" \
            "${extra_args[@]}" \
            -o "$scratch_assembly"
    fi

    scratch_assembly_valid "$scratch_assembly" || \
        die "MEGAHIT did not produce a valid final.contigs.fa: $species/$project/$sample"

    write_scratch_done "$scratch_assembly"
    prune_scratch_assembly "$scratch_assembly"

    rm -f "$sr1" "$sr2"
    rm -rf "$sample_input"

    elapsed=$(( $(date +%s) - start ))
    log "MEGAHIT sample completed: $species/$project/$sample (${elapsed}s)"
}

species_task() {
    ensure_run_identity
    load_megahit
    need rsync
    need gzip
    need stat
    need find
    need awk

    local task="${SLURM_ARRAY_TASK_ID:-${2:-}}"
    [[ -n "$task" ]] || die "Missing species task ID"

    local line task_id species expected issues source_species scratch_species final_species
    line="$(manifest_row "$SPECIES_MANIFEST" "$task")"
    [[ -n "$line" ]] || die "Invalid species task ID: $task"

    IFS=$'\t' read -r task_id species expected issues source_species scratch_species final_species <<< "$line"
    (( issues == 0 )) || die "Species has pairing issues: $species"
    (( expected > 0 )) || die "Species has no paired samples: $species"

    mkdir -p "$scratch_species" "$SAMPLE_STATUS_DIR" "$SPECIES_STATUS_DIR"

    local status_file="$SPECIES_STATUS_DIR/$(safe_name "${task_id}_${species}").tsv"
    local start
    local files_deleted=0
    local bytes_deleted=0
    start="$(date +%s)"

    write_species_status() {
        local status="$1"
        local note="$2"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$task_id" "$species" "$expected" "$files_deleted" "$bytes_deleted" \
            "$status" "$(( $(date +%s) - start ))" "$scratch_species" "$note" \
            > "$status_file"
    }

    log "Species $species; expected samples=$expected; scratch=$scratch_species"
    megahit --version 2>&1 | head -n 2 >&2 || true

    local sample sp project r1 r2 scratch_assembly
    local final_contigs final_log final_options final_done
    local sample_status_file sample_start contigs bytes

    # Phase 1: assemble each sample in scratch, or reuse valid final/scratch data.
    while IFS=$'\t' read -r sample sp project r1 r2 scratch_assembly \
        final_contigs final_log final_options final_done; do

        [[ "$sp" == "$species" ]] || continue

        sample_status_file="$SAMPLE_STATUS_DIR/$(safe_name "${species}_${project}_${sample}").tsv"
        sample_start="$(date +%s)"

        if [[ "$OVERWRITE" != 1 ]] && final_output_valid "$final_contigs" "$final_done"; then
            contigs="$(grep -c '^>' "$final_contigs")"
            bytes="$(stat -c '%s' "$final_contigs")"
            printf '%s\t%s\t%s\t%s\t%s\t%s\tREUSED_FINAL\t0\n' \
                "$sample" "$species" "$project" "$contigs" "$bytes" \
                "$final_contigs" > "$sample_status_file"
            continue
        fi

        run_megahit_sample "$sample" "$species" "$project" "$r1" "$r2" "$scratch_assembly"

        contigs="$(grep -c '^>' "$scratch_assembly/final.contigs.fa")"
        bytes="$(stat -c '%s' "$scratch_assembly/final.contigs.fa")"
        printf '%s\t%s\t%s\t%s\t%s\t%s\tSUCCESS_SCRATCH\t%s\n' \
            "$sample" "$species" "$project" "$contigs" "$bytes" \
            "$scratch_assembly/final.contigs.fa" \
            "$(( $(date +%s) - sample_start ))" > "$sample_status_file"
    done < <(tail -n +2 "$SAMPLE_MANIFEST")

    # Phase 2: every sample must have a complete final or scratch assembly.
    local completed=0
    while IFS=$'\t' read -r sample sp project r1 r2 scratch_assembly \
        final_contigs final_log final_options final_done; do

        [[ "$sp" == "$species" ]] || continue

        if final_output_valid "$final_contigs" "$final_done"; then
            :
        elif scratch_assembly_valid "$scratch_assembly"; then
            :
        else
            write_species_status "RETAINED_INCOMPLETE_ASSEMBLY" "$project/$sample"
            die "Incomplete assembly: $species/$project/$sample"
        fi
        ((completed += 1))
    done < <(tail -n +2 "$SAMPLE_MANIFEST")

    if (( completed != expected )); then
        write_species_status "RETAINED_COMPLETION_MISMATCH" "$completed/$expected complete"
        die "Completion count mismatch for $species: $completed/$expected"
    fi

    # Phase 3: after all sample assemblies are complete, delete any remaining
    # paired nonhost FASTQ files for this species. Partial prior deletion is safe
    # during resume because all assemblies have already been validated above.
    if [[ "$DELETE_NONHOST_READS" == 1 ]]; then
        local -a delete_list=()
        local f size

        while IFS=$'\t' read -r sample sp project r1 r2 scratch_assembly \
            final_contigs final_log final_options final_done; do

            [[ "$sp" == "$species" ]] || continue
            [[ -e "$r1" ]] && delete_list+=("$r1")
            [[ -e "$r2" ]] && delete_list+=("$r2")
        done < <(tail -n +2 "$SAMPLE_MANIFEST")

        for f in "${delete_list[@]}"; do
            size="$(stat -c '%s' "$f")"
            bytes_deleted=$((bytes_deleted + size))
        done
        files_deleted="${#delete_list[@]}"

        if (( files_deleted > 0 )); then
            rm -f -- "${delete_list[@]}"
            log "Deleted $files_deleted remaining nonhost FASTQ files for $species"
        else
            log "No nonhost FASTQ files remain for $species; continuing resume"
        fi
    fi

    # Phase 4: transfer any missing scratch assemblies to flat final outputs.
    while IFS=$'\t' read -r sample sp project r1 r2 scratch_assembly \
        final_contigs final_log final_options final_done; do

        [[ "$sp" == "$species" ]] || continue

        if final_output_valid "$final_contigs" "$final_done"; then
            continue
        fi

        if ! scratch_assembly_valid "$scratch_assembly"; then
            write_species_status "TRANSFER_FAILED_SCRATCH_MISSING" "$project/$sample"
            die "Scratch assembly missing after source deletion: $species/$project/$sample"
        fi

        copy_scratch_to_flat_output \
            "$scratch_assembly" "$sample" "$final_contigs" "$final_log" \
            "$final_options" "$final_done"
    done < <(tail -n +2 "$SAMPLE_MANIFEST")

    # Phase 5: validate all final outputs, update sample reports, then remove the
    # species scratch directory.
    completed=0
    while IFS=$'\t' read -r sample sp project r1 r2 scratch_assembly \
        final_contigs final_log final_options final_done; do

        [[ "$sp" == "$species" ]] || continue

        if ! final_output_valid "$final_contigs" "$final_done"; then
            write_species_status "FINAL_VALIDATION_FAILED_SCRATCH_RETAINED" "$project/$sample"
            die "Final assembly validation failed: $species/$project/$sample"
        fi

        sample_status_file="$SAMPLE_STATUS_DIR/$(safe_name "${species}_${project}_${sample}").tsv"
        contigs="$(grep -c '^>' "$final_contigs")"
        bytes="$(stat -c '%s' "$final_contigs")"
        local old_elapsed=0
        if [[ -s "$sample_status_file" ]]; then
            old_elapsed="$(awk -F '\t' 'NR == 1 {print $8+0}' "$sample_status_file")"
        fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\tSUCCESS_FINAL\t%s\n' \
            "$sample" "$species" "$project" "$contigs" "$bytes" \
            "$final_contigs" "$old_elapsed" > "$sample_status_file"

        ((completed += 1))
    done < <(tail -n +2 "$SAMPLE_MANIFEST")

    if (( completed != expected )); then
        write_species_status "FINAL_COMPLETION_MISMATCH_SCRATCH_RETAINED" "$completed/$expected"
        die "Final completion mismatch for $species"
    fi

    rm -rf "$scratch_species"

    if [[ "$DELETE_NONHOST_READS" == 1 ]]; then
        write_species_status "SUCCESS_DELETED_NONHOST_READS" "Final: $final_species"
    else
        write_species_status "SUCCESS_RETAINED_NONHOST_READS" "Final: $final_species"
    fi

    log "Completed species $species -> $final_species"
}

collect_reports() {
    mkdir -p "$REPORT_DIR"

    printf 'sample_id\tspecies\tproject\tcontig_count\tfinal_contigs_bytes\tfinal_contigs\tstatus\telapsed_seconds\n' > "$SAMPLE_SUMMARY"
    local sample_files=("$SAMPLE_STATUS_DIR"/*.tsv)
    if (( ${#sample_files[@]} > 0 )); then
        LC_ALL=C sort "${sample_files[@]}" >> "$SAMPLE_SUMMARY"
    fi

    printf 'species_task_id\tspecies\texpected_samples\tfiles_deleted\tbytes_deleted\tstatus\telapsed_seconds\tscratch_directory\tnote\n' > "$SPECIES_SUMMARY"
    local species_files=("$SPECIES_STATUS_DIR"/*.tsv)
    if (( ${#species_files[@]} > 0 )); then
        LC_ALL=C sort -t $'\t' -k1,1n "${species_files[@]}" >> "$SPECIES_SUMMARY"
    fi

    {
        printf 'WORK_ROOT: %s\n' "$WORK_ROOT"
        printf 'RUN_NAME: %s\n' "$RUN_NAME"
        printf 'Final assemblies: %s\n' "$ASSEMBLY_ROOT"
        printf 'Scratch root: %s\n\n' "$RUN_SCRATCH"
        printf 'Species statuses:\n'
        awk -F '\t' 'NR > 1 {n[$6]++} END {for (k in n) print k "\t" n[k]}' \
            "$SPECIES_SUMMARY" | LC_ALL=C sort
    } > "$RUN_SUMMARY"

    cat "$RUN_SUMMARY"
}

export_runtime_vars() {
    export SCRIPT_VERSION WORK_ROOT PROJECT_ID PROJECT_DIR_NAME AUTO_RUN_NAME RUN_NAME
    export READS_ROOT OUT_ROOT ASSEMBLY_ROOT REPORT_DIR LOG_DIR
    export SCRATCH_ROOT RUN_SCRATCH RUN_IDENTITY_FILE
    export SBATCH_ACCOUNT SBATCH_PARTITION MEGAHIT_CPUS MEGAHIT_MEM MEGAHIT_TIME
    export MAX_SPECIES_CONCURRENT MEGAHIT_THREADS MIN_CONTIG_LEN MEGAHIT_EXTRA_ARGS
    export DELETE_NONHOST_READS VERIFY_INPUT_GZIP OVERWRITE PRUNE_INTERMEDIATE
}

submit_pipeline() {
    need sbatch

    cd "$WORK_ROOT"
    mkdir -p "$LOG_DIR" "$REPORT_DIR"

    build_manifests
    ensure_run_identity

    local species_count
    species_count=$(( $(wc -l < "$SPECIES_MANIFEST") - 1 ))
    (( species_count > 0 )) || \
        die "No runnable species. Check $MANIFEST_SUMMARY and $PAIRING_ISSUES"

    export_runtime_vars

    local submitted species_job collect_job
    submitted="$(
        sbatch --parsable \
            --export=ALL \
            --account="$SBATCH_ACCOUNT" \
            --partition="$SBATCH_PARTITION" \
            --job-name=megahit_species \
            --cpus-per-task="$MEGAHIT_CPUS" \
            --mem="$MEGAHIT_MEM" \
            --time="$MEGAHIT_TIME" \
            --array="1-${species_count}%${MAX_SPECIES_CONCURRENT}" \
            --output="$LOG_DIR/megahit_species_%A_%a.out" \
            --error="$LOG_DIR/megahit_species_%A_%a.err" \
            "$SCRIPT_PATH" species
    )"
    species_job="${submitted%%;*}"

    submitted="$(
        sbatch --parsable \
            --export=ALL \
            --account="$SBATCH_ACCOUNT" \
            --partition="$SBATCH_PARTITION" \
            --job-name=megahit_collect \
            --cpus-per-task=1 \
            --mem=2G \
            --time=00:15:00 \
            --dependency="afterany:$species_job" \
            --output="$LOG_DIR/megahit_collect_%j.out" \
            --error="$LOG_DIR/megahit_collect_%j.err" \
            "$SCRIPT_PATH" collect
    )"
    collect_job="${submitted%%;*}"

    printf '2 Slurm job IDs submitted:\n'
    printf '  species array: %s (%s species; max %s concurrently)\n' \
        "$species_job" "$species_count" "$MAX_SPECIES_CONCURRENT"
    printf '  collect:       %s\n' "$collect_job"
}

resume_pipeline() {
    need sbatch

    [[ -s "$SPECIES_MANIFEST" && -s "$SAMPLE_MANIFEST" ]] || \
        die "Existing manifests are required for resume: $REPORT_DIR"

    cd "$WORK_ROOT"
    mkdir -p "$LOG_DIR" "$REPORT_DIR"
    ensure_run_identity

    local species_count
    species_count=$(( $(wc -l < "$SPECIES_MANIFEST") - 1 ))
    (( species_count > 0 )) || die "Species manifest contains no tasks"

    export_runtime_vars

    local submitted species_job collect_job
    submitted="$(
        sbatch --parsable \
            --export=ALL \
            --account="$SBATCH_ACCOUNT" \
            --partition="$SBATCH_PARTITION" \
            --job-name=megahit_resume \
            --cpus-per-task="$MEGAHIT_CPUS" \
            --mem="$MEGAHIT_MEM" \
            --time="$MEGAHIT_TIME" \
            --array="1-${species_count}%${MAX_SPECIES_CONCURRENT}" \
            --output="$LOG_DIR/megahit_resume_%A_%a.out" \
            --error="$LOG_DIR/megahit_resume_%A_%a.err" \
            "$SCRIPT_PATH" species
    )"
    species_job="${submitted%%;*}"

    submitted="$(
        sbatch --parsable \
            --export=ALL \
            --account="$SBATCH_ACCOUNT" \
            --partition="$SBATCH_PARTITION" \
            --job-name=megahit_collect \
            --cpus-per-task=1 \
            --mem=2G \
            --time=00:15:00 \
            --dependency="afterany:$species_job" \
            --output="$LOG_DIR/megahit_collect_%j.out" \
            --error="$LOG_DIR/megahit_collect_%j.err" \
            "$SCRIPT_PATH" collect
    )"
    collect_job="${submitted%%;*}"

    printf '2 resume job IDs submitted:\n'
    printf '  species array: %s (%s species; max %s concurrently)\n' \
        "$species_job" "$species_count" "$MAX_SPECIES_CONCURRENT"
    printf '  collect:       %s\n' "$collect_job"
}

show_help() {
    cat <<EOF_HELP
Usage (run from the project directory containing 02_host_removed_reads):
  bash $(basename "$SCRIPT_PATH") manifest
  bash $(basename "$SCRIPT_PATH") submit
  bash $(basename "$SCRIPT_PATH") resume
  bash $(basename "$SCRIPT_PATH") collect

Automatic paths for the current directory:
  WORK_ROOT=$WORK_ROOT
  RUN_NAME=$RUN_NAME
  READS_ROOT=$READS_ROOT
  RUN_SCRATCH=$RUN_SCRATCH
  OUT_ROOT=$OUT_ROOT

Input layout:
  <WORK_ROOT>/02_host_removed_reads/reads/<species>/<project>/<sample>_1.nonhost.fastq.gz
  <WORK_ROOT>/02_host_removed_reads/reads/<species>/<project>/<sample>_2.nonhost.fastq.gz

Scratch layout:
  <RUN_SCRATCH>/<species>/results/<project>/<sample>/final.contigs.fa

Final flat output layout:
  <WORK_ROOT>/03_megahit_sample_assemblies/assemblies/<species>/<project>/<sample>_final.contigs.fa
  <WORK_ROOT>/03_megahit_sample_assemblies/assemblies/<species>/<project>/<sample>_megahit.log
  <WORK_ROOT>/03_megahit_sample_assemblies/assemblies/<species>/<project>/<sample>_options.json
  <WORK_ROOT>/03_megahit_sample_assemblies/assemblies/<species>/<project>/<sample>_megahit.done

Preview without submitting:
  bash $(basename "$SCRIPT_PATH") manifest

Recommended submission, maximum 3 species concurrently:
  MAX_SPECIES_CONCURRENT=3 \\
  DELETE_NONHOST_READS=1 \\
  VERIFY_INPUT_GZIP=1 \\
  bash $(basename "$SCRIPT_PATH") submit

Resume an interrupted run with the existing manifests:
  MAX_SPECIES_CONCURRENT=3 \\
  bash $(basename "$SCRIPT_PATH") resume

Useful optional settings:
  WORK_ROOT=/explicit/project/path       # overrides current directory
  RUN_NAME=unique_run_name               # overrides automatic name
  SCRATCH_ROOT=/path/to/scratch         # parent scratch directory
  MAX_SPECIES_CONCURRENT=3
  MIN_CONTIG_LEN=1000
  MEGAHIT_CPUS=32
  MEGAHIT_THREADS=32
  MEGAHIT_MEM=128G
  MEGAHIT_TIME=7-00:00:00
  SBATCH_ACCOUNT=your_account
  SBATCH_PARTITION=sapphire
  DELETE_NONHOST_READS=1
  VERIFY_INPUT_GZIP=1
  PRUNE_INTERMEDIATE=1
  OVERWRITE=0
EOF_HELP
}

case "${1:-help}" in
    manifest|prepare)
        build_manifests
        ;;
    species|run)
        species_task "$@"
        ;;
    collect)
        collect_reports
        ;;
    submit)
        submit_pipeline
        ;;
    resume)
        resume_pipeline
        ;;
    help|-h|--help)
        show_help
        ;;
    *)
        show_help
        die "Unknown mode: ${1:-}"
        ;;
esac
