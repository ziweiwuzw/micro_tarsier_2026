#!/usr/bin/env bash
#SBATCH --job-name=semibin2_prepare_sp
#
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --partition=fos
#SBATCH --qos=fos
#SBATCH --time=14-00:00:00
#SBATCH --output=slurm-%x-%A_%a.out
#SBATCH --error=slurm-%x-%A_%a.err
set -Eeuo pipefail
BASE="${BASE:-/path/to/project}"
SCRATCH="${SCRATCH:-/path/to/scratch}"
MANIFEST="${MANIFEST:-$BASE/06_comebin/per_sample_manifest.tsv}"
ROOT="${ROOT:-$SCRATCH/07_semibin2/by_species}"
ENV_PREFIX="${SEMIBIN_ENV:-$SCRATCH/software/anaconda3/envs/SemiBin2}"
CONDA_SH="${CONDA_SH:-/path/to/software/Anaconda3/2024.02-1/etc/profile.d/conda.sh}"
MAP="${MAP:-$BASE/07_semibin2/species_group_map.tsv}"
mkdir -p "$ROOT" "$BASE/07_semibin2/logs"
[[ -s "$MANIFEST" ]] || { echo "Missing manifest: $MANIFEST" >&2; exit 2; }
source "$CONDA_SH"; conda activate "$ENV_PREFIX"
# Build explicit canonical species mapping; unknown names remain unchanged.
python - "$MANIFEST" "$MAP" "$ROOT/sample_inputs.tsv" <<'PY'
import csv, pathlib, sys
manifest, mapping, out = map(pathlib.Path, sys.argv[1:])
mp={}
if mapping.exists():
    with mapping.open() as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if r['original_species'].strip(): mp[r['original_species'].strip()]=r['canonical_species'].strip()
with manifest.open() as f, out.open('w') as g:
    rows=csv.DictReader(f, delimiter='\t'); w=csv.writer(g, delimiter='\t', lineterminator='\n')
    w.writerow(['array_index','sample_id','original_species','canonical_species','contigs','r1','r2'])
    for r in rows:
        sp=r['species']; canon=mp.get(sp, sp)
        for k in ('contigs','r1','r2'):
            if not pathlib.Path(r[k]).is_file(): raise SystemExit(f"Missing {k}: {r[k]}")
        w.writerow([r['array_index'],r['sample_id'],sp,canon,r['contigs'],r['r1'],r['r2']])
PY
cut -f4 "$ROOT/sample_inputs.tsv" | tail -n +2 | sort -u | while read -r group; do
  safe=$(printf '%s' "$group" | tr -c 'A-Za-z0-9._-' '_'); dir="$ROOT/$safe"; mkdir -p "$dir"
  awk -F '\t' -v g="$group" 'NR>1 && $4==g {print $5}' "$ROOT/sample_inputs.tsv" > "$dir/contigs.list"
  mapfile -t fs < "$dir/contigs.list"
  [[ -s "$dir/contigs.list" ]] || continue
  if [[ ! -s "$dir/concatenated.fa.gz" ]]; then SemiBin2 concatenate_fasta -i "${fs[@]}" -o "$dir"; fi
done
printf 'groups=%s\nroot=%s\n' "$(cut -f4 "$ROOT/sample_inputs.tsv" | tail -n +2 | sort -u | wc -l)" "$ROOT" | tee "$ROOT/prepare.done"
