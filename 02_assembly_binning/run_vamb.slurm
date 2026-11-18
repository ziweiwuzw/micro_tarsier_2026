#!/usr/bin/env bash
# Run CPU-only VAMB using the same coordinate-sorted BAMs consumed by SemiBin2.
# Required exports: JOINT_ROOT, GROUP_INDEX, GROUPS_TSV.
#SBATCH --job-name=joint_vamb
#
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00
#SBATCH --output=slurm-%x-%A_%a.out
#SBATCH --error=slurm-%x-%A_%a.err
set -Eeuo pipefail
trap 'rc=$?; echo "[FAILED] group=${GROUP_INDEX:-?} line=${LINENO} status=${rc}" >&2; exit "$rc"' ERR

ROOT="${JOINT_ROOT:?Set JOINT_ROOT}"
GROUPS_MANIFEST="${GROUPS_TSV:-$ROOT/joint_groups.tsv}"
INDEX="${GROUP_INDEX:?Set GROUP_INDEX}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-32}}"
VAMB_ENV="${VAMB_ENV:-/path/to/home/software/anaconda3/envs/Taxvamb}"
CONDA_SH="${CONDA_SH:-/path/to/software/Anaconda3/2024.02-1/etc/profile.d/conda.sh}"
EXPECTED_VAMB_VERSION="${EXPECTED_VAMB_VERSION:-5.0.4}"
MINLEN="${MINLEN:-2000}"
MINFASTA="${MINFASTA:-200000}"
SEED="${SEED:-42}"

[[ -s "$GROUPS_MANIFEST" ]] || { echo "Missing groups manifest: $GROUPS_MANIFEST" >&2; exit 2; }
row=$(awk -F '\t' -v n="$INDEX" 'NR > 1 && $1 == n {print; count++} END {if (count != 1) exit 1}' "$GROUPS_MANIFEST") || { echo "No unique group $INDEX" >&2; exit 2; }
IFS=$'\t' read -r _ group group_slug group_dir reference _ bam_dir _ vamb_out sample_count _ <<< "$row"
[[ -s "$reference" && -d "$bam_dir" ]] || { echo "Missing reference/BAM directory for $group_slug" >&2; exit 2; }
mapfile -t bams < <(find "$bam_dir" -maxdepth 1 -type f -name '*.sorted.bam' -print | sort)
(( ${#bams[@]} == sample_count )) || { echo "Expected $sample_count BAMs, found ${#bams[@]} for $group_slug" >&2; exit 3; }
for bam in "${bams[@]}"; do
    [[ -s "$bam" && -s "$bam.bai" ]] || { echo "Missing BAM/BAI: $bam" >&2; exit 3; }
done

if type module >/dev/null 2>&1; then module purge; fi
source "$CONDA_SH"
conda activate "$VAMB_ENV"
export PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME PIP_USER || true
command -v vamb >/dev/null || { echo "vamb not found" >&2; exit 127; }
actual_version=$(python -c 'from importlib.metadata import version; print(version("vamb"))')
[[ -z "$EXPECTED_VAMB_VERSION" || "$actual_version" == "$EXPECTED_VAMB_VERSION" ]] || { echo "Expected VAMB $EXPECTED_VAMB_VERSION, found $actual_version" >&2; exit 4; }

done_file="$vamb_out/.joint_vamb.done"
if [[ -s "$done_file" && -s "$vamb_out/vae_clusters_split.tsv" && -d "$vamb_out/bins" ]]; then
    echo "[REUSE] VAMB output: $vamb_out"
    exit 0
fi
[[ ! -e "$vamb_out" ]] || { echo "Incomplete VAMB output exists; inspect/remove: $vamb_out" >&2; exit 5; }

vamb bin default --outdir "$vamb_out" --fasta "$reference" --bamdir "$bam_dir" \
    -m "$MINLEN" -p "$THREADS" -o ':' --seed "$SEED" --minfasta "$MINFASTA"

[[ -s "$vamb_out/vae_clusters_split.tsv" && -d "$vamb_out/bins" ]] || { echo "VAMB output validation failed: $vamb_out" >&2; exit 6; }
printf 'status=SUCCESS\ngroup=%s\nreference=%s\nbam_count=%s\nresult=%s\ncompleted_at=%s\n' \
    "$group" "$reference" "${#bams[@]}" "$vamb_out" "$(date -Is)" > "$done_file"
