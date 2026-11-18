#!/usr/bin/env bash
#SBATCH --job-name=sgb_tables
#SBATCH --partition=cascade,sapphire
#SBATCH --qos=publiccpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=slurm-%x-%j.out
set -Eeuo pipefail
source "${ANNOT_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}/common.sh"
unset PYTHONPATH PYTHONHOME
export PYTHONNOUSERSITE=1
require_file "$ANN_OUT/input_manifest.tsv"
require_file "$ANN_OUT/03_dram2_all_members/.done"
python3 "$SCRIPT_DIR/export_ko.py" --out "$ANN_OUT" --dram-dir 03_dram2_all_members
python3 "$SCRIPT_DIR/export_dram_features.py" --out "$ANN_OUT"
python3 "$SCRIPT_DIR/export_cazy.py" --out "$ANN_OUT"
python3 "$SCRIPT_DIR/export_virulence.py" --out "$ANN_OUT" \
    --toxin-metadata "$DB_ROOT/toxinome/all_genes.csv"
printf 'status=SUCCESS\n' > "$ANN_OUT/12_final_tables/.done"
echo "[SUCCESS] Source-specific tables: $ANN_OUT/12_final_tables"
