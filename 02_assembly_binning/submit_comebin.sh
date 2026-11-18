#!/usr/bin/env bash
# Submit CPU COMEBin for the selected per-sample manifest rows.
set -Eeuo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="${BASE:-/path/to/project}"
MANIFEST="${MANIFEST:-${BASE}/06_comebin/per_sample_manifest.tsv}"
INDEX_FILE="${INDEX_FILE:-}"
LOG_DIR="${LOG_DIR:-${BASE}/06_comebin/logs}"
MAX_CONCURRENT="${MAX_CONCURRENT:-2}"

fail() { echo "ERROR: $*" >&2; exit 1; }
command -v sbatch >/dev/null || fail 'sbatch is required'
[[ -s "$MANIFEST" && -s "$HERE/run_comebin.slurm" && -s "$HERE/prepare_comebin_mapping.slurm" ]] \
    || fail 'Missing manifest or COMEBin scripts'
[[ "$MAX_CONCURRENT" =~ ^[1-9][0-9]*$ ]] || fail 'MAX_CONCURRENT must be positive'
mkdir -p "$LOG_DIR"
if [[ -z "$INDEX_FILE" ]]; then
    INDEX_FILE="$LOG_DIR/comebin_indices_$(date +%Y%m%d_%H%M%S).txt"
    awk -F '\t' 'NR>1 {print $1}' "$MANIFEST" > "$INDEX_FILE"
fi
[[ -s "$INDEX_FILE" ]] || fail "Empty index file: $INDEX_FILE"
awk 'NF && $1 !~ /^#/ && $1 !~ /^[1-9][0-9]*$/ {exit 1}' "$INDEX_FILE" || fail 'Index file contains a non-positive index'
count=$(awk 'NF && $1 !~ /^#/ {n++} END {print n+0}' "$INDEX_FILE")
(( count > 0 )) || fail 'No COMEBin targets selected'
unique=$(awk 'NF && $1 !~ /^#/ {print $1}' "$INDEX_FILE" | sort -u | wc -l | tr -d '[:space:]')
(( unique == count )) || fail 'Index file contains duplicate array indices'
job=$(sbatch --parsable --job-name=comebin_cpu --partition=sapphire --qos=publiccpu \
    --cpus-per-task=64 --mem=128G --time=4-00:00:00 \
    --array="1-${count}%${MAX_CONCURRENT}" \
    --output="$LOG_DIR/%x_%A_%a.out" --error="$LOG_DIR/%x_%A_%a.err" \
    --export="ALL,BASE=${BASE},MANIFEST=${MANIFEST},INDEX_FILE=${INDEX_FILE}" \
    "$HERE/run_comebin.slurm")
printf 'COMEBin job: %s\nTargets: %s\nIndex file: %s\n' "${job%%;*}" "$count" "$INDEX_FILE"
