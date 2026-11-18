#!/usr/bin/env bash
#SBATCH --job-name=sgb_dbcan
#SBATCH --partition=cascade,sapphire
#SBATCH --qos=publiccpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --output=slurm-%x-%A_%a.out
set -Eeuo pipefail
source "${ANNOT_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}/common.sh"
require_file "$ANN_OUT/input_manifest.tsv"
load_annotation_env
command -v run_dbcan >/dev/null || die 'run_dbcan missing in annotation environment'
db="$DB_ROOT/dbcan3"
for name in CAZy.dmnd dbCAN.hmm dbCAN-sub.hmm fam-substrate-mapping.tsv; do
    require_file "$db/$name"
done
mag="$(mag_at_task)"
faa="$ANN_OUT/01_orf/$mag/$mag.faa"
require_file "$faa"
final="$ANN_OUT/04_dbcan/$mag"
if [[ -s "$final/.done" ]]; then
    require_file "$final/overview.tsv"
    echo "[REUSE] $final"
    exit 0
fi
[[ ! -e "$final" ]] || die "Partial dbCAN output exists: $final"
mkdir -p "$ANN_OUT/04_dbcan"
tmp="$(mktemp -d "$ANN_OUT/04_dbcan/$mag.tmp.XXXXXX")"
trap 'rc=$?; if ((rc != 0)); then echo "[ERROR] dbCAN partial output retained: '$tmp'" >&2; fi' EXIT
run_dbcan CAZyme_annotation --mode protein --input_raw_data "$faa" \
    --output_dir "$tmp" --db_dir "$db" --threads "${SLURM_CPUS_PER_TASK:-8}" \
    --methods diamond,hmm,dbCANsub
require_file "$tmp/overview.tsv"
[[ -f "$tmp/dbCAN_hmm_results.tsv" ]] || die 'dbCAN HMM result missing'
[[ -f "$tmp/dbCANsub_hmm_results.tsv" ]] || die 'dbCAN subfamily HMM result missing'
dbcan_version="$(python3 -c 'from importlib.metadata import version; print(version("dbcan"))' 2>/dev/null || printf unknown)"
printf 'dbcan_version=%s\n' "$dbcan_version" > "$tmp/.done"
mv "$tmp" "$final"
echo "[SUCCESS] $final"
