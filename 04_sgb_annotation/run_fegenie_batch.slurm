#!/usr/bin/env bash
#SBATCH --job-name=sgb_fegenie
#SBATCH --partition=sapphire
#SBATCH --qos=publiccpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=7-00:00:00
#SBATCH --output=slurm-%x-%A_%a.out
set -Eeuo pipefail
source "${ANNOT_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}/common.sh"

FEGENIE_ROOT="${FEGENIE_ROOT:-$ANN_ROOT/software/FeGenie}"
FEGENIE_ENV="${FEGENIE_ENV:-/path/to/scratch/software/anaconda3/envs/fegenie}"
run_name="${FEGENIE_RUN_NAME:-$(basename "$ANN_OUT")}"
result_stack_root=$(dirname "$(dirname "$ANN_OUT")")
runtime="${FEGENIE_RUNTIME:-$result_stack_root/runtime/$run_name/fegenie}"
stride="${FEGENIE_BATCH_STRIDE:-500}"
total=$(awk -F '\t' 'END {print NR-1}' "$ANN_OUT/input_manifest.tsv")
task="${SLURM_ARRAY_TASK_ID:-}"

[[ "$task" =~ ^[1-9][0-9]*$ ]] || die "Invalid SLURM_ARRAY_TASK_ID"
require_file "$ANN_ROOT/software/fegenie_install/.done"
require_file "$ANN_OUT/diagnostics/fegenie_pilot/.done"
mkdir -p "$runtime/tasks" "$ANN_OUT/07_fegenie" "$ANN_OUT/diagnostics/fegenie_failures"

module load Anaconda3/2024.02-1
eval "$(conda shell.bash hook)"
conda activate "$FEGENIE_ENV"
export PATH="$FEGENIE_ROOT:$PATH"
export rscripts="$FEGENIE_ROOT/rscripts"
export iron_hmms="$FEGENIE_ROOT/hmms/iron"
export MPLBACKEND=Agg

for ((index=task; index<=total; index+=stride)); do
    mag=$(awk -F '\t' -v n="$index" 'NR==n+1 {print $1; exit}' "$ANN_OUT/input_manifest.tsv")
    [[ -n "$mag" ]] || die "No MAG at manifest index $index"
    genome="$ANN_OUT/01_orf/$mag/$mag.fa"
    final="$ANN_OUT/07_fegenie/$mag"
    require_file "$genome"
    if [[ -s "$final/.done" ]]; then
        echo "[REUSE] index=$index MAG=$mag"
        continue
    fi
    [[ ! -e "$final" ]] || die "Partial FeGenie result exists: $final"

    tmp=$(mktemp -d "$runtime/tasks/${mag}.tmp.XXXXXX")
    mkdir -p "$tmp/input" "$tmp/output"
    ln -s "$genome" "$tmp/input/$mag.fa"
    log="$tmp/FeGenie.log"
    if ! python "$FEGENIE_ROOT/FeGenie.py" \
        -bin_dir "$tmp/input" \
        -bin_ext fa \
        -out "$tmp/output" \
        -t "${SLURM_CPUS_PER_TASK:-4}" \
        --meta > "$log" 2>&1; then
        cp "$log" "$ANN_OUT/diagnostics/fegenie_failures/$mag.log"
        rm -rf "$tmp"
        die "FeGenie failed for $mag"
    fi
    if grep -Eqi 'bad file format|Traceback|cannot find the correct hmmsearch|segmentation fault' "$log"; then
        cp "$log" "$ANN_OUT/diagnostics/fegenie_failures/$mag.log"
        rm -rf "$tmp"
        die "FeGenie fatal diagnostic for $mag"
    fi
    for file in FinalSummary.csv FinalSummary-dereplicated-clustered.csv FeGenie-geneSummary-clusters.csv; do
        [[ -f "$tmp/output/$file" ]] || die "Missing FeGenie output for $mag: $file"
    done

    temporary="$final.tmp.$$"
    mkdir -p "$temporary"
    for file in \
        FinalSummary.csv \
        FinalSummary-dereplicated-clustered.csv \
        FeGenie-geneSummary.csv \
        FeGenie-geneSummary-clusters.csv \
        FeGenie-heatmap-data.csv; do
        [[ -f "$tmp/output/$file" ]] && cp "$tmp/output/$file" "$temporary/"
    done
    printf 'status=SUCCESS\nmanifest_index=%s\nmode=nucleotide_fasta_meta\n' "$index" \
        > "$temporary/.done"
    mv "$temporary" "$final"
    rm -rf "$tmp"
    echo "[SUCCESS] index=$index MAG=$mag"
done

