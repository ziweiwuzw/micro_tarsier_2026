#!/usr/bin/env bash
#SBATCH --partition=cascade,sapphire
#SBATCH --qos=publiccpu
#SBATCH --output=slurm-%x-%j.out
# Resource requests are set per mode by submit.sh.
set -Eeuo pipefail
HERE="${PIPELINE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "$HERE/common.sh"
require_file "$HERE/workflow.py"
mode="${1:-}"
case "$mode" in
    index|sample)
        [[ -x "$KMA_ENV/bin/kma" ]] || die "KMA executable missing: $KMA_ENV/bin/kma"
        export PATH="$KMA_ENV/bin:$PATH"
        ;;
    merge) ;;
    *) die 'Usage: kma_job.slurm index|sample|merge' ;;
esac
if [[ "$mode" == sample ]]; then
    : "${SLURM_ARRAY_TASK_ID:?Submit sample mode as a Slurm array}"
    python3 "$HERE/workflow.py" sample "$SLURM_ARRAY_TASK_ID"
else
    python3 "$HERE/workflow.py" "$mode"
fi
