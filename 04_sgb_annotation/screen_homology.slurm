#!/usr/bin/env bash
#SBATCH --job-name=sgb_homology
#SBATCH --partition=cascade,sapphire
#SBATCH --qos=publiccpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=2-00:00:00
#SBATCH --output=slurm-%x-%A_%a.out
set -Eeuo pipefail
source "${ANNOT_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}/common.sh"
load_annotation_env
command -v diamond >/dev/null || die 'DIAMOND missing in annotation environment'
require_file "$ANN_OUT/05_reference/.done"
mag="$(mag_at_task)"
faa="$ANN_OUT/01_orf/$mag/$mag.faa"
require_file "$faa"
final="$ANN_OUT/06_homology/$mag"
if [[ -s "$final/.done" ]]; then
    require_file "$final/best_hits.tsv"
    echo "[REUSE] $final"
    exit 0
fi
[[ ! -e "$final" ]] || die "Partial homology output: $final"
mkdir -p "$ANN_OUT/06_homology"
tmp="$(mktemp -d "$ANN_OUT/06_homology/$mag.tmp.XXXXXX")"
mkdir -p "$tmp/raw"
outfmt=(6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore)
diamond blastp --query "$faa" --db "$ANN_OUT/05_reference/vfdb_A" \
    --out "$tmp/raw/vfdb_A.tsv" --outfmt "${outfmt[@]}" --threads "${SLURM_CPUS_PER_TASK:-8}" \
    --very-sensitive --max-target-seqs 25 --id 70 --query-cover 70 --subject-cover 70 --evalue 1e-10
diamond blastp --query "$faa" --db "$ANN_OUT/05_reference/vfdb_B" \
    --out "$tmp/raw/vfdb_B.tsv" --outfmt "${outfmt[@]}" --threads "${SLURM_CPUS_PER_TASK:-8}" \
    --very-sensitive --max-target-seqs 25 --id 70 --query-cover 70 --subject-cover 70 --evalue 1e-10
diamond blastp --query "$faa" --db "$ANN_OUT/05_reference/toxinome" \
    --out "$tmp/raw/toxinome.tsv" --outfmt "${outfmt[@]}" --threads "${SLURM_CPUS_PER_TASK:-8}" \
    --sensitive --max-target-seqs 25 --id 40 --query-cover 90 --subject-cover 60 --evalue 0.001
python3 "$SCRIPT_DIR/homology.py" parse --mag "$mag" --catalog "$ANN_OUT/05_reference/catalog.sqlite" \
    --raw-dir "$tmp/raw" --out "$tmp"
require_file "$tmp/best_hits.tsv"
printf 'vfdb=id70,qcov70,scov70,e1e-10\ntoxinome=id40,qcov90,scov60,e0.001\n' > "$tmp/.done"
mv "$tmp" "$final"
echo "[SUCCESS] $final"
