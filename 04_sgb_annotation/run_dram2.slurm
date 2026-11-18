#!/usr/bin/env bash
#SBATCH --job-name=sgb_dram2
#SBATCH --partition=cascade,sapphire
#SBATCH --qos=publiccpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=14-00:00:00
#SBATCH --output=slurm-%x-%j.out
set -Eeuo pipefail
source "${ANNOT_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}/common.sh"
require_file "$ANN_OUT/input_manifest.tsv"
require_file "$DRAM_ROOT/dram.config"
require_file "$DRAM_ROOT/databases/dram_db/dram_db.hmm"
require_file "$DRAM_ROOT/databases/dram_db/dram_db_scores.tsv"
require_file "$DRAM_ROOT/databases/db_descriptions/description_db.sqlite"
for db in kofam merops metals sulfur card; do
    require_dir "$DRAM_ROOT/databases/$db"
done
dram_out="$ANN_OUT/03_dram2_all_members"
genes_dir="$ANN_OUT/02_dram_input_all_members/genes"
if [[ -s "$dram_out/.done" ]]; then
    echo "[REUSE] $dram_out"
    exit 0
fi
resume_args=()
if [[ -e "$dram_out" ]]; then
    [[ "${DRAM_RESUME:-0}" == 1 ]] || die "Partial DRAM2 output exists; inspect it and set DRAM_RESUME=1 to resume"
    resume_args=(-resume)
fi
module load Nextflow/25.10.4
module load Apptainer/1.3.3
command -v nextflow >/dev/null || die 'Nextflow missing'
command -v apptainer >/dev/null || die 'Apptainer missing'
mkdir -p "$genes_dir" "$ANN_OUT/work"
while IFS=$'\t' read -r mag sgb rep genome genome_sha256; do
    [[ "$mag" == mag_id ]] && continue
    orf="$ANN_OUT/01_orf/$mag"
    require_file "$orf/.done"
    for suffix in faa fna; do
        require_file "$orf/$mag.$suffix"
    done
    ln -sfn "$orf/$mag.faa" "$genes_dir/$mag.faa"
    ln -sfn "$orf/$mag.fna" "$genes_dir/$mag.fna"
done < "$ANN_OUT/input_manifest.tsv"
input_count=$(find "$genes_dir" -maxdepth 1 -name '*.faa' | wc -l | tr -d '[:space:]')
expected=$(awk -F '\t' 'END {print NR-1}' "$ANN_OUT/input_manifest.tsv")
[[ "$input_count" == "$expected" ]] || die "Expected $expected MAG inputs, found $input_count"
cd "$DRAM_ROOT"
nextflow run WrightonLabCSU/DRAM -r "${DRAM_REVISION:-4536cd93c9}" \
    "${resume_args[@]}" \
    -profile apptainer,slurm -c "$DRAM_ROOT/dram.config" \
    -work-dir "$ANN_OUT/work" \
    --input_genes "$genes_dir/" \
    --outdir "$dram_out" \
    --annotate --summarize --sum_ecos gut \
    --anno_dbs dram_db,kofam,merops,metals,sulfur,card
require_file "$dram_out/ANNOTATE/raw-annotations.tsv"
diff -u \
    <(awk -F '\t' 'NR > 1 {print $1}' "$ANN_OUT/input_manifest.tsv" | sort -u) \
    <(awk -F '\t' 'NR > 1 {print $2}' "$dram_out/ANNOTATE/raw-annotations.tsv" | sort -u) ||
    die 'DRAM2 input_fasta values do not cover exactly the selected MAGs'
printf 'revision=%s\nMAGs=%s\nscope=all_members\n' "${DRAM_REVISION:-4536cd93c9}" "$expected" > "$dram_out/.done"
