#!/usr/bin/env bash
#SBATCH --job-name=sgb_gtdbtk
#SBATCH --partition=cascade,sapphire
#SBATCH --qos=publiccpu
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --time=14-00:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

# Final taxonomy is intentionally independent of MAG recovery and SGB creation.
set -Eeuo pipefail
trap 'rc=$?; echo "[FAILED] line=${LINENO} status=${rc}" >&2; exit "$rc"' ERR
HERE="${PIPELINE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
[[ -f "$HERE/common.sh" ]] || { echo "[ERROR] common.sh not found under script directory: $HERE" >&2; exit 1; }
# shellcheck disable=SC1091
source "$HERE/common.sh"

require_file "${SGB_ROOT}/sgb_representatives.tsv"
require_dir "$SGB_REP_DIR"
require_file "${GTDBTK_REFERENCE_PATH}/metadata/metadata.txt"

out="${SGB_ROOT}/taxonomy_gtdbtk"
run_out="${out}/classify_wf"
taxonomy="${out}/sgb_taxonomy_gtdbtk.tsv"
done_file="${out}/.taxonomy.done"
if [[ -s "$done_file" && -s "$taxonomy" ]]; then
    echo "[REUSE] $taxonomy"
    exit 0
fi
[[ ! -e "$out" ]] || die "Incomplete GTDB-Tk output exists: $out"
mkdir -p "$out" "${out}/tmp"

# GTDB-Tk is a hierarchical Spartan module. Load its compiler/MPI parents
# explicitly so the batch job does not depend on the submitter's environment.
module purge
module load GCC/11.3.0 OpenMPI/4.1.4 "$GTDBTK_MODULE"
command -v gtdbtk >/dev/null || die "gtdbtk not found after loading $GTDBTK_MODULE"
# The module sets its own central path. Restore the workflow's explicit,
# version-matched reference path so inherited/module state cannot redirect it.
GTDBTK_DATA_PATH="$GTDBTK_REFERENCE_PATH"
export GTDBTK_DATA_PATH
require_file "${GTDBTK_DATA_PATH}/metadata/metadata.txt"
export TMPDIR="${out}/tmp"
export TMP="$TMPDIR" TEMP="$TMPDIR"
mkdir -p "$(dirname "$GTDBTK_MASH_DB")"
gtdbtk --version
echo "GTDBTK_DATA_PATH=$GTDBTK_DATA_PATH"
echo "GTDBTK_MASH_DB=$GTDBTK_MASH_DB"

gtdbtk classify_wf \
    --genome_dir "$SGB_REP_DIR" \
    --extension fa \
    --out_dir "$run_out" \
    --prefix gtdbtk \
    --mash_db "$GTDBTK_MASH_DB" \
    --cpus "${SLURM_CPUS_PER_TASK:-64}" \
    --pplacer_cpus "${GTDBTK_PPLACER_CPUS:-32}" \
    --tmpdir "$TMPDIR"

summary_args=()
for summary in \
    "${run_out}/gtdbtk.bac120.summary.tsv" \
    "${run_out}/gtdbtk.ar53.summary.tsv"; do
    [[ -s "$summary" ]] && summary_args+=(--summary "$summary")
done
(( ${#summary_args[@]} > 0 )) || die "GTDB-Tk produced no bacterial or archaeal summary"

python3 "$HERE/collect_gtdbtk_taxonomy.py" \
    --representatives "${SGB_ROOT}/sgb_representatives.tsv" \
    "${summary_args[@]}" \
    --output "$taxonomy"

require_file "$taxonomy"
printf 'status=SUCCESS\nmethod=GTDB-Tk-v2.4.0-R220\nreference=%s\nmash_db=%s\ntaxonomy=%s\ncompleted_at=%s\n' \
    "$GTDBTK_DATA_PATH" "$GTDBTK_MASH_DB" "$taxonomy" "$(date -Is)" > "$done_file"
echo "[SUCCESS] GTDB-Tk taxonomy=$taxonomy"
