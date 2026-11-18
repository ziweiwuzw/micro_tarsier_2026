#!/usr/bin/env bash
# Prepare all host groups, then submit only the MetaWRAP/direct-bin stage.
set -Eeuo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PIPELINE_SCRIPT_DIR="$HERE"
: "${PIPELINE_ROOT:?Set PIPELINE_ROOT to a fresh candidate-bin output root}"
export SOURCE_ROOT="$PIPELINE_ROOT"
source "$HERE/../../03_mag_sgb_kma/common.sh"
command -v sbatch >/dev/null || die 'sbatch unavailable; submit on Spartan'
: "${JOINT_ROOT:?Set JOINT_ROOT to the finished group-binning root}"
[[ "$PIPELINE_ROOT" != "$JOINT_ROOT" ]] || die 'PIPELINE_ROOT must differ from JOINT_ROOT'
[[ ! -e "$PIPELINE_ROOT" ]] || die "Choose a fresh PIPELINE_ROOT: $PIPELINE_ROOT"
args=(--joint-root "$JOINT_ROOT" --output-root "$PIPELINE_ROOT" --all-groups)
[[ -n "${EXCLUDE_COMEBIN_GROUPS:-}" ]] && args+=(--exclude-comebin-groups "$EXCLUDE_COMEBIN_GROUPS")
[[ "${ALLOW_UNMARKED:-0}" == 1 ]] && args+=(--allow-unmarked)
python3 "$HERE/prepare_refinement_inputs.py" "${args[@]}"
require_file "$PILOT_GROUPS_TSV"
count=$(awk 'END {print NR-1}' "$PILOT_GROUPS_TSV")
(( count == 32 )) || die "Expected 32 host groups; found $count"
max="${MAX_CONCURRENT_REFINEMENTS:-3}"
[[ "$max" =~ ^[1-9][0-9]*$ ]] || die 'MAX_CONCURRENT_REFINEMENTS must be positive'
mkdir -p "$LOG_ROOT"
sbatch --parsable --export=ALL --array="1-${count}%${max}" \
    --output="$LOG_ROOT/%x_%A_%a.out" --error="$LOG_ROOT/%x_%A_%a.err" \
    "$HERE/run_metawrap_refinement.slurm"
