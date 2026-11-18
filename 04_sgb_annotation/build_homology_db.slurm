#!/usr/bin/env bash
#SBATCH --job-name=sgb_refdb
#SBATCH --partition=cascade,sapphire
#SBATCH --qos=publiccpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00
#SBATCH --output=slurm-%x-%j.out
set -Eeuo pipefail
source "${ANNOT_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}/common.sh"
load_annotation_env
command -v diamond >/dev/null || die 'DIAMOND missing in annotation environment'
vfdb_a="$DB_ROOT/vfdb/VFDB_setA_pro.fas"
vfdb_b="$DB_ROOT/vfdb/VFDB_setB_pro.fas"
toxin_fasta="$DB_ROOT/toxinome/toxinome_genes.fasta"
toxin_csv="$DB_ROOT/toxinome/all_genes.csv"
for file in "$vfdb_a" "$vfdb_b" "$toxin_fasta" "$toxin_csv"; do require_file "$file"; done
final="$ANN_OUT/05_reference"
if [[ -s "$final/.done" ]]; then
    for name in catalog.sqlite vfdb_A.dmnd vfdb_B.dmnd toxinome.dmnd; do require_file "$final/$name"; done
    echo "[REUSE] $final"
    exit 0
fi
[[ ! -e "$final" ]] || die "Partial reference build: $final"
mkdir -p "$ANN_OUT"
tmp="$(mktemp -d "$ANN_OUT/05_reference.tmp.XXXXXX")"
python3 "$SCRIPT_DIR/homology.py" prepare --vfdb-a "$vfdb_a" --vfdb-b "$vfdb_b" \
    --toxin-fasta "$toxin_fasta" --toxin-csv "$toxin_csv" --out "$tmp"
diamond makedb --in "$vfdb_a" --db "$tmp/vfdb_A" --threads "${SLURM_CPUS_PER_TASK:-16}"
diamond makedb --in "$vfdb_b" --db "$tmp/vfdb_B" --threads "${SLURM_CPUS_PER_TASK:-16}"
diamond makedb --in "$toxin_fasta" --db "$tmp/toxinome" --threads "${SLURM_CPUS_PER_TASK:-16}"
for name in catalog.sqlite vfdb_A.dmnd vfdb_B.dmnd toxinome.dmnd; do require_file "$tmp/$name"; done
printf 'vfdb_A=%s\nvfdb_B=%s\ntoxinome_fasta=%s\ntoxinome_csv=%s\n' \
    "$vfdb_a" "$vfdb_b" "$toxin_fasta" "$toxin_csv" > "$tmp/.done"
mv "$tmp" "$final"
echo "[SUCCESS] $final"
