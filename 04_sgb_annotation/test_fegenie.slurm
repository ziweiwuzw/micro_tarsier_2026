#!/usr/bin/env bash
#SBATCH --job-name=fegenie_test
#SBATCH --partition=sapphire
#SBATCH --qos=publiccpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --output=slurm-%x-%j.out
set -Eeuo pipefail
source "${ANNOT_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}/common.sh"

FEGENIE_ROOT="${FEGENIE_ROOT:-$ANN_ROOT/software/FeGenie}"
FEGENIE_ENV="${FEGENIE_ENV:-/path/to/scratch/software/anaconda3/envs/fegenie}"
run_name="${FEGENIE_RUN_NAME:-$(basename "$ANN_OUT")}"
result_stack_root=$(dirname "$(dirname "$ANN_OUT")")
runtime="${FEGENIE_RUNTIME:-$result_stack_root/runtime/$run_name/fegenie}"
pilot="$ANN_OUT/diagnostics/fegenie_pilot"

require_file "$ANN_ROOT/software/fegenie_install/.done"
require_file "$ANN_OUT/input_manifest.tsv"
if [[ -s "$pilot/.done" ]]; then
    echo "[REUSE] $pilot"
    cat "$pilot/.done"
    exit 0
fi

module load Anaconda3/2024.02-1
eval "$(conda shell.bash hook)"
conda activate "$FEGENIE_ENV"
export PATH="$FEGENIE_ROOT:$PATH"
export rscripts="$FEGENIE_ROOT/rscripts"
export iron_hmms="$FEGENIE_ROOT/hmms/iron"
export MPLBACKEND=Agg

mkdir -p "$runtime"
tmp=$(mktemp -d "$runtime/pilot.tmp.XXXXXX")
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/input" "$tmp/output"

mapfile -t pilot_mags < <(
    awk -F '\t' '
        NR == 1 {next}
        $1 ~ /^Cephalopachus_bancanus__/ && !target {print $1; target=1; next}
        !other && $1 !~ /^Cephalopachus_bancanus__/ {print $1; other=1}
        target && other {exit}
    ' "$ANN_OUT/input_manifest.tsv"
)
if (( ${#pilot_mags[@]} < 2 )); then
    mapfile -t pilot_mags < <(awk -F '\t' 'NR>1 && NR<=3 {print $1}' "$ANN_OUT/input_manifest.tsv")
fi
(( ${#pilot_mags[@]} == 2 )) || die "Unable to select two pilot MAGs"

printf 'mag_id\n' > "$tmp/pilot_mags.tsv"
for mag in "${pilot_mags[@]}"; do
    genome="$ANN_OUT/01_orf/$mag/$mag.fa"
    require_file "$genome"
    ln -s "$genome" "$tmp/input/$mag.fa"
    printf '%s\n' "$mag" >> "$tmp/pilot_mags.tsv"
done

python "$FEGENIE_ROOT/FeGenie.py" \
    -bin_dir "$tmp/input" \
    -bin_ext fa \
    -out "$tmp/output" \
    -t "${SLURM_CPUS_PER_TASK:-4}" \
    --meta > "$tmp/FeGenie-pilot.log" 2>&1

for file in FinalSummary.csv FinalSummary-dereplicated-clustered.csv FeGenie-geneSummary-clusters.csv; do
    [[ -f "$tmp/output/$file" ]] || die "Missing FeGenie pilot output: $file"
done
if grep -Eqi 'bad file format|Traceback|cannot find the correct hmmsearch|segmentation fault' \
        "$tmp/FeGenie-pilot.log"; then
    die "FeGenie pilot log contains a fatal diagnostic"
fi

temporary="$pilot.tmp.$$"
rm -rf "$temporary"
mkdir -p "$temporary"
cp "$tmp/pilot_mags.tsv" "$tmp/FeGenie-pilot.log" "$temporary/"
find "$tmp/output" -maxdepth 1 -type f \
    \( -name '*.csv' -o -name '*.faa' -o -name '*.txt' \) \
    -exec cp {} "$temporary/" \;
{
    printf 'status=SUCCESS\n'
    printf 'MAGs=%s\n' "${#pilot_mags[@]}"
    printf 'mode=nucleotide_fasta_meta\n'
    printf 'operon_context=enabled\n'
    printf 'completed_at=%s\n' "$(date -Is)"
} > "$temporary/.done"
rm -rf "$pilot"
mv "$temporary" "$pilot"
trap - EXIT
rm -rf "$tmp"

cat "$pilot/.done"
printf 'Pilot MAGs:\n'
tail -n +2 "$pilot/pilot_mags.tsv"
printf 'Output headers:\n'
for file in FinalSummary.csv FinalSummary-dereplicated-clustered.csv FeGenie-geneSummary-clusters.csv; do
    printf '%s: ' "$file"
    head -n 1 "$pilot/$file"
done

