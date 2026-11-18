#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:?Set PIPELINE_ROOT to the MAG and SGB results directory}"
ANN_ROOT="${ANN_ROOT:-/path/to/scratch/annotation_stack}"
ANN_RESULTS_ROOT="${ANN_RESULTS_ROOT:-$ANN_ROOT/results}"
ANN_OUT="${ANN_OUT:-$ANN_RESULTS_ROOT/sgb_annotation_all_species_3binners_v2_pan_metawrap_kma}"
DB_ROOT="${DB_ROOT:-$ANN_ROOT/databases}"
DRAM_ROOT="${DRAM_ROOT:-$ANN_ROOT/dram2}"
ANNOT_ENV="${ANNOT_ENV:-/path/to/scratch/software/anaconda3/envs/micro_annota}"
GAPMIND_ROOT="${GAPMIND_ROOT:-$ANN_ROOT/software/PaperBLAST}"
export SCRIPT_DIR PIPELINE_ROOT ANN_ROOT ANN_RESULTS_ROOT ANN_OUT DB_ROOT DRAM_ROOT ANNOT_ENV GAPMIND_ROOT

die() { printf '[ERROR] %s\n' "$*" >&2; exit 1; }
require_file() { [[ -s "$1" ]] || die "Missing or empty: $1"; }
require_dir() { [[ -d "$1" ]] || die "Missing directory: $1"; }

load_annotation_env() {
    module load Anaconda3/2024.02-1
    eval "$(conda shell.bash hook)"
    conda activate "$ANNOT_ENV"
    unset PYTHONPATH PYTHONHOME
    export PYTHONNOUSERSITE=1
    hash -r
}

mag_at_task() {
    local index="${SLURM_ARRAY_TASK_ID:-}"
    [[ "$index" =~ ^[1-9][0-9]*$ ]] || die 'SLURM_ARRAY_TASK_ID must be a positive integer'
    local offset="${ARRAY_OFFSET:-0}"
    [[ "$offset" =~ ^[0-9]+$ ]] || die 'ARRAY_OFFSET must be a non-negative integer'
    index=$((index + offset))
    awk -F '\t' -v n="$index" 'NR == n+1 { print $1; found=1 } END { if (!found) exit 1 }' "$ANN_OUT/input_manifest.tsv"
}

sgb_at_task() {
    local index="${SLURM_ARRAY_TASK_ID:-}"
    [[ "$index" =~ ^[1-9][0-9]*$ ]] || die 'SLURM_ARRAY_TASK_ID must be a positive integer'
    local offset="${ARRAY_OFFSET:-0}"
    [[ "$offset" =~ ^[0-9]+$ ]] || die 'ARRAY_OFFSET must be a non-negative integer'
    index=$((index + offset))
    awk -F '\t' -v n="$index" 'NR > 1 && !seen[$2]++ { i++; if (i == n) {print $2; found=1; exit} } END { if (!found) exit 1 }' "$ANN_OUT/input_manifest.tsv"
}
