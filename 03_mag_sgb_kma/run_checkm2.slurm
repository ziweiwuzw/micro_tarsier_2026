#!/usr/bin/env bash
#SBATCH --job-name=mag_checkm2
#SBATCH --partition=cascade,sapphire
#SBATCH --qos=publiccpu
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --time=14-00:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -Eeuo pipefail
trap 'rc=$?; echo "[FAILED] line=${LINENO} status=${rc}" >&2; exit "$rc"' ERR
HERE="${PIPELINE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
[[ -f "$HERE/common.sh" ]] || { echo "[ERROR] common.sh not found under script directory: $HERE" >&2; exit 1; }
# shellcheck disable=SC1091
source "$HERE/common.sh"

require_file "$CANDIDATE_MANIFEST"
require_dir "$CANDIDATE_DIR"
quality_report="${CHECKM2_OUT}/quality_report.tsv"
done_file="${CHECKM2_OUT}/.checkm2.done"
if [[ -s "$done_file" && -s "$quality_report" ]]; then
    echo "[REUSE] $quality_report"
    exit 0
fi
[[ ! -e "$CHECKM2_OUT" ]] || die "Incomplete CheckM2 output exists: $CHECKM2_OUT"

module purge
module load CheckM2/1.0.1
load_checkm2_database
command -v checkm2 >/dev/null || die "checkm2 not found"
checkm2_tmp="/tmp/magc2_${USER:-user}_${SLURM_JOB_ID:-$$}"
mkdir -p "$checkm2_tmp"
chmod 700 "$checkm2_tmp"
export TMPDIR="$checkm2_tmp" TMP="$checkm2_tmp" TEMP="$checkm2_tmp"
echo "CheckM2 temporary directory: $TMPDIR"
args=(predict --threads "${SLURM_CPUS_PER_TASK:-64}" --input "$CANDIDATE_DIR"
    --output-directory "$CHECKM2_OUT" -x fa)
[[ -n "$CHECKM2_DB" ]] && args+=(--database_path "$CHECKM2_DB")
[[ "${CHECKM2_LOWMEM:-0}" == 1 ]] && args+=(--lowmem)
checkm2 "${args[@]}"

require_file "$quality_report"
printf 'status=SUCCESS\nquality_report=%s\ncompleted_at=%s\n' \
    "$quality_report" "$(date -Is)" > "$done_file"
echo "[SUCCESS] CheckM2 genomes=$(awk 'NR>1{n++} END{print n+0}' "$quality_report")"
