#!/bin/bash
#SBATCH --job-name=toga2_cafe5
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --time=5-00:00:00
#SBATCH --output=cafe5-%j.out
#SBATCH --error=cafe5-%j.err

set -euo pipefail

if [[ $# -lt 3 ]]; then
    echo "Usage: sbatch $0 <family_counts.txt> <time_tree.nwk> <results_dir>" >&2
    exit 2
fi

FAMILY_FILE=$(realpath "$1")
TREE_FILE=$(realpath "$2")
RESULTS_DIR=$(realpath -m "$3")

CAFE5_BIN=${CAFE5_BIN:-/path/to/scratch/software/biosoft/CAFE5/build/bin/cafe5}
MODEL=${MODEL:-gamma}         # base or gamma
K=${K:-3}                    # used only when MODEL=gamma
CPU=${SLURM_CPUS_PER_TASK:-6}

[[ -s "$FAMILY_FILE" ]] || { echo "Missing family file: $FAMILY_FILE" >&2; exit 1; }
[[ -s "$TREE_FILE" ]] || { echo "Missing tree file: $TREE_FILE" >&2; exit 1; }
[[ -x "$CAFE5_BIN" ]] || { echo "CAFE5 is not executable: $CAFE5_BIN" >&2; exit 1; }
[[ "$MODEL" == "base" || "$MODEL" == "gamma" ]] || {
    echo "MODEL must be base or gamma, got: $MODEL" >&2
    exit 1
}

mkdir -p "$RESULTS_DIR"
cd "$RESULTS_DIR"

if [[ "$MODEL" == "base" ]]; then
    RUN_NAME=base
    CMD=("$CAFE5_BIN" -i "$FAMILY_FILE" -t "$TREE_FILE" -p -c "$CPU" -o "$RUN_NAME")
else
    [[ "$K" =~ ^[1-9][0-9]*$ ]] || { echo "K must be a positive integer" >&2; exit 1; }
    RUN_NAME="gamma_k${K}"
    CMD=("$CAFE5_BIN" -i "$FAMILY_FILE" -t "$TREE_FILE" -p -c "$CPU" -k "$K" -o "$RUN_NAME")
fi

{
    echo "date=$(date -Is)"
    echo "host=$(hostname)"
    echo "family_file=$FAMILY_FILE"
    echo "tree_file=$TREE_FILE"
    echo "model=$MODEL"
    echo "gamma_categories=$K"
    echo "threads=$CPU"
    printf 'command='
    printf '%q ' "${CMD[@]}"
    printf '\n'
} | tee "${RUN_NAME}.command.txt"

srun "${CMD[@]}" 2>&1 | tee "${RUN_NAME}.log"

if [[ "$MODEL" == "base" ]]; then
    EXPECTED_RESULT="$RUN_NAME/Base_results.txt"
else
    EXPECTED_RESULT="$RUN_NAME/Gamma_results.txt"
fi

if [[ ! -s "$EXPECTED_RESULT" ]]; then
    echo "ERROR: CAFE5 exited without producing $EXPECTED_RESULT" >&2
    echo "Inspect: $RESULTS_DIR/${RUN_NAME}.log" >&2
    exit 1
fi

echo "CAFE5 completed: $RESULTS_DIR/$RUN_NAME"
