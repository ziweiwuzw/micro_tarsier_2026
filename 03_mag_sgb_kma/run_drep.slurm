#!/usr/bin/env bash
#SBATCH --job-name=mag_drep
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

genome_info="$(dirname "$HQ_DIR")/drep_genome_info.csv"
require_file "$HQ_MANIFEST"
require_file "$genome_info"
require_dir "$HQ_DIR"
done_file="${DREP_OUT}/.drep.done"
if [[ -s "$done_file" && -d "${DREP_OUT}/dereplicated_genomes" ]]; then
    echo "[REUSE] $DREP_OUT"
    exit 0
fi
[[ ! -e "$DREP_OUT" ]] || die "Incomplete dRep output exists: $DREP_OUT"
mapfile -t genomes < <(find "$HQ_DIR" -maxdepth 1 -type l -name '*.fa' -print | sort)
(( ${#genomes[@]} >= 2 )) || die "dRep pilot requires at least two quality-filtered MAGs; found ${#genomes[@]}"

module purge
module load GCC/11.3.0 OpenMPI/4.1.4 dRep/3.6.2
command -v dRep >/dev/null || die "dRep not found"
command -v mash >/dev/null || die "mash not found"
command -v fastANI >/dev/null || die "fastANI not found"
dRep dereplicate "$DREP_OUT" \
    -g "${genomes[@]}" \
    -p "${SLURM_CPUS_PER_TASK:-64}" \
    --genomeInfo "$genome_info" \
    -comp "${MIN_COMPLETENESS:-90}" \
    -con "${MAX_CONTAMINATION:-5}" \
    -pa "${PRIMARY_ANI:-0.90}" \
    -sa "${SECONDARY_ANI:-0.95}" \
    -nc "${MIN_ALIGNED_FRACTION:-0.30}" \
    --S_algorithm fastANI

require_dir "${DREP_OUT}/dereplicated_genomes"
require_file "${DREP_OUT}/data_tables/Cdb.csv"
rep_count=$(find "${DREP_OUT}/dereplicated_genomes" -maxdepth 1 -type f | wc -l | tr -d '[:space:]')
(( rep_count > 0 )) || die "dRep returned no representative genomes"
printf 'status=SUCCESS\nrepresentatives=%s\nsecondary_ani=%s\nmin_aligned_fraction=%s\ncompleted_at=%s\n' \
    "$rep_count" "${SECONDARY_ANI:-0.95}" "${MIN_ALIGNED_FRACTION:-0.30}" "$(date -Is)" > "$done_file"
echo "[SUCCESS] dRep representatives=$rep_count"
