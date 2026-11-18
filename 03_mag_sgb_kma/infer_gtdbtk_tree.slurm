#!/usr/bin/env bash
# Infer separate, unrooted bacterial and archaeal trees from SGB marker MSAs.
# Usage: sbatch infer_gtdbtk_tree.slurm /absolute/path/pipeline_root
# Reuses classify_wf user MSAs only when they contain every SGB representative.
#SBATCH --job-name=sgb_tree
#SBATCH --partition=cascade,sapphire
#SBATCH --qos=publiccpu
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --time=7-00:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -Eeuo pipefail
trap 'rc=$?; echo "[FAILED] line=${LINENO} exit=${rc}" >&2; exit "$rc"' ERR

root="${1:-}"
[[ "$root" == /* && -d "$root" ]] || {
    echo 'Usage: sbatch infer_gtdbtk_tree.slurm /absolute/path/pipeline_root' >&2
    exit 2
}
root="$(cd "$root" && pwd -P)"
reps="$root/sgb/sgb_representatives.tsv"
genomes="$root/sgb/representatives"
taxonomy="$root/sgb/taxonomy_gtdbtk/sgb_taxonomy_gtdbtk.tsv"
[[ -s "$reps" && -d "$genomes" && -s "$taxonomy" ]] || {
    echo '[ERROR] Missing SGB representatives or completed GTDB-Tk taxonomy' >&2
    exit 1
}
out="$root/sgb/phylogeny_gtdbtk"
if [[ -s "$out/.done" ]]; then
    echo "[REUSE] Completed SGB phylogeny: $out"
    exit 0
fi
mkdir -p "$out/tmp" "$out/msas"

# Match the GTDB-Tk version and R220 database used for SGB classification.
module purge
module load GCC/11.3.0 OpenMPI/4.1.4 GTDB-Tk/2.4.0
command -v gtdbtk >/dev/null || { echo '[ERROR] GTDB-Tk not on PATH' >&2; exit 1; }
command -v python3 >/dev/null || { echo '[ERROR] python3 not on PATH' >&2; exit 1; }
reference="${GTDBTK_REFERENCE_PATH:-/path/to/GTDBtk/release220}"
[[ -s "$reference/metadata/metadata.txt" ]] || {
    echo "[ERROR] GTDB-Tk reference missing: $reference" >&2
    exit 1
}
export GTDBTK_DATA_PATH="$reference"
export TMPDIR="$out/tmp" TMP="$out/tmp" TEMP="$out/tmp"
cpus="${SLURM_CPUS_PER_TASK:-64}"
gtdbtk --version
echo "[INFO] GTDBTK_DATA_PATH=$GTDBTK_DATA_PATH"

# Print COMPLETE only when the two user-only alignments contain every SGB
# exactly once. Also write an audit with explicit exclusions and taxonomy.
audit_msa() {
    local alignment_dir="$1" report="$2"
    python3 - "$reps" "$taxonomy" "$alignment_dir" "$report" <<'PY'
import csv
import gzip
import sys
from pathlib import Path

representatives, taxonomy_file, alignment_dir, report = map(Path, sys.argv[1:])
with representatives.open(newline='') as handle:
    expected_rows = list(csv.DictReader(handle, delimiter='\t'))
expected = [row['sgb_id'] for row in expected_rows]
if not expected or len(expected) != len(set(expected)):
    raise SystemExit('[ERROR] Empty or duplicate SGB representative IDs')
with taxonomy_file.open(newline='') as handle:
    labels = {row['sgb_id']: row.get('domain', '') for row in csv.DictReader(handle, delimiter='\t')}
if set(labels) != set(expected):
    raise SystemExit('[ERROR] Taxonomy and SGB representatives have different IDs')

membership = {}
for domain in ('bac120', 'ar53'):
    path = alignment_dir / f'gtdbtk.{domain}.user_msa.fasta.gz'
    if not path.is_file():
        continue
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            if line.startswith('>'):
                genome = line[1:].strip().split()[0]
                if genome in membership:
                    raise SystemExit(f'[ERROR] Duplicate alignment genome: {genome}')
                membership[genome] = domain
unexpected = sorted(set(membership) - set(expected))
if unexpected:
    raise SystemExit(f'[ERROR] Unexpected alignment genome IDs: {unexpected[:5]}')

report.parent.mkdir(parents=True, exist_ok=True)
with report.open('w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
    writer.writerow(('sgb_id', 'taxonomy_domain', 'marker_set', 'status'))
    for genome in expected:
        domain = membership.get(genome, '')
        writer.writerow((genome, labels[genome], domain, 'IN_MSA' if domain else 'NOT_IN_MSA'))
missing = len(expected) - len(membership)
print(f'[AUDIT] expected={len(expected)} in_MSA={len(membership)} missing={missing}', file=sys.stderr)
print('COMPLETE' if not missing else 'INCOMPLETE')
PY
}

classify_align="$root/sgb/taxonomy_gtdbtk/classify_wf/align"
if [[ "$(audit_msa "$classify_align" "$out/classify_msa_audit.tsv")" == COMPLETE ]]; then
    alignment_dir="$classify_align"
    alignment_origin='existing_classify_wf'
    echo '[REUSE] classify_wf marker alignment contains every SGB'
else
    identify_dir="$out/identify"
    align_dir="$out/align"
    if [[ ! -s "$out/.identify.done" ]]; then
        [[ ! -e "$identify_dir" ]] || {
            echo "[ERROR] Incomplete identify output exists: $identify_dir" >&2
            exit 1
        }
        gtdbtk identify \
            --genome_dir "$genomes" --extension fa \
            --out_dir "$identify_dir" --prefix gtdbtk \
            --cpus "$cpus" --force --tmpdir "$TMPDIR"
        [[ -s "$identify_dir/identify/gtdbtk.bac120.markers_summary.tsv" ||
           -s "$identify_dir/identify/gtdbtk.ar53.markers_summary.tsv" ]] || {
            echo '[ERROR] GTDB-Tk identify produced no marker summary' >&2
            exit 1
        }
        printf 'SUCCESS\n' > "$out/.identify.done"
    fi
    if [[ ! -s "$out/.align.done" ]]; then
        [[ ! -e "$align_dir" ]] || {
            echo "[ERROR] Incomplete align output exists: $align_dir" >&2
            exit 1
        }
        gtdbtk align \
            --identify_dir "$identify_dir" --out_dir "$align_dir" \
            --skip_gtdb_refs --prefix gtdbtk \
            --cpus "$cpus" --tmpdir "$TMPDIR"
        [[ -s "$align_dir/align/gtdbtk.bac120.user_msa.fasta.gz" ||
           -s "$align_dir/align/gtdbtk.ar53.user_msa.fasta.gz" ]] || {
            echo '[ERROR] GTDB-Tk align produced no user MSA' >&2
            exit 1
        }
        printf 'SUCCESS\n' > "$out/.align.done"
    fi
    alignment_dir="$align_dir/align"
    alignment_origin='new_identify_align'
fi

[[ "$(audit_msa "$alignment_dir" "$out/sgb_tree_inclusion.tsv")" == COMPLETE ]] || {
    echo '[ERROR] GTDB-Tk marker MSAs do not contain every SGB; inspect sgb_tree_inclusion.tsv' >&2
    exit 1
}
trees=0
for domain in bac120 ar53; do
    source_msa="$alignment_dir/gtdbtk.${domain}.user_msa.fasta.gz"
    [[ -s "$source_msa" ]] || continue
    msa="$out/msas/${domain}.user_msa.fasta"
    gzip -dc "$source_msa" > "$msa"
    count="$(awk '/^>/{n++} END{print n+0}' "$msa")"
    if (( count < 2 )); then
        echo "[SKIP] $domain has only $count aligned SGB; at least two are needed"
        continue
    fi
    infer_dir="$out/infer_${domain}"
    tree="$infer_dir/gtdbtk_${domain}.unrooted.tree"
    marker="$out/.infer_${domain}.done"
    if [[ ! -s "$marker" ]]; then
        [[ ! -e "$infer_dir" ]] || {
            echo "[ERROR] Incomplete infer output exists: $infer_dir" >&2
            exit 1
        }
        gtdbtk infer \
            --msa_file "$msa" --out_dir "$infer_dir" \
            --prefix "gtdbtk_${domain}" --prot_model WAG --gamma \
            --cpus "$cpus" --tmpdir "$TMPDIR"
        [[ -s "$tree" ]] || { echo "[ERROR] Missing tree: $tree" >&2; exit 1; }
        printf 'SUCCESS\n' > "$marker"
    fi
    python3 - "$msa" "$tree" <<'PY'
import re
import sys
from collections import Counter
from pathlib import Path

msa, tree = map(Path, sys.argv[1:])
expected = [line[1:].strip().split()[0] for line in msa.open() if line.startswith('>')]
observed = re.findall(r'SGB_[0-9]+', tree.read_text())
if Counter(expected) != Counter(observed):
    raise SystemExit(f'[ERROR] Tree tips differ from the marker MSA: {tree}')
print(f'[PASS] {tree.name}: {len(observed)} SGB tips')
PY
    trees=$((trees + 1))
done
(( trees > 0 )) || { echo '[ERROR] No domain had enough aligned SGBs for inference' >&2; exit 1; }
printf 'status=SUCCESS\nmethod=GTDB-Tk identify/align/infer\nreference=%s\nalignment_origin=%s\ntrees=%s\ncompleted_at=%s\n' \
    "$reference" "$alignment_origin" "$trees" "$(date -Is)" > "$out/.done"
echo "[SUCCESS] $trees unrooted SGB tree(s) and inclusion audit: $out"
