#!/usr/bin/env bash
# Run SemiBin2 in CPU or GPU mode for one group after mapping succeeds.
# Resource placement is selected by the submitter, not hard-coded here.
# Required exports: JOINT_ROOT, GROUP_INDEX, GROUPS_TSV.
# Optional: SEMIBIN_ENGINE=cpu|gpu, THREADS, SEMIBIN_ENV.
#SBATCH --job-name=joint_semibin2
#
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=7-00:00:00
#SBATCH --output=slurm-%x-%A_%a.out
#SBATCH --error=slurm-%x-%A_%a.err
set -Eeuo pipefail
trap 'rc=$?; echo "[FAILED] group=${GROUP_INDEX:-?} line=${LINENO} status=${rc}" >&2; exit "$rc"' ERR

ROOT="${JOINT_ROOT:?Set JOINT_ROOT}"
GROUPS_MANIFEST="${GROUPS_TSV:-$ROOT/joint_groups.tsv}"
INDEX="${GROUP_INDEX:?Set GROUP_INDEX}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"
SEMIBIN_ENGINE="${SEMIBIN_ENGINE:-cpu}"
SEMIBIN_ENV="${SEMIBIN_ENV:-/path/to/scratch/software/anaconda3/envs/SemiBin2}"
CONDA_SH="${CONDA_SH:-/path/to/software/Anaconda3/2024.02-1/etc/profile.d/conda.sh}"
SEMIBIN_PRELOAD_LIBSTDCXX="${SEMIBIN_PRELOAD_LIBSTDCXX:-0}"
SEED="${SEED:-42}"

[[ "$SEMIBIN_ENGINE" == cpu || "$SEMIBIN_ENGINE" == gpu ]] || {
    echo "SEMIBIN_ENGINE must be cpu or gpu" >&2
    exit 2
}
[[ -s "$GROUPS_MANIFEST" ]] || { echo "Missing groups manifest: $GROUPS_MANIFEST" >&2; exit 2; }
row=$(awk -F '\t' -v n="$INDEX" 'NR > 1 && $1 == n {print; count++} END {if (count != 1) exit 1}' "$GROUPS_MANIFEST") || {
    echo "No unique group $INDEX" >&2
    exit 2
}
IFS=$'\t' read -r _ group group_slug group_dir reference _ bam_dir semibin_out _ sample_count _ <<< "$row"
[[ -s "$reference" && -d "$bam_dir" ]] || { echo "Missing reference/BAM directory for $group_slug" >&2; exit 2; }

if type module >/dev/null 2>&1; then module purge; fi
source "$CONDA_SH"
conda activate "$SEMIBIN_ENV"
export PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME PIP_USER || true

# Keep the Conda C++ runtime fix from the tested CPU script. It is also needed
# by optree/torch on some GPU compute nodes.
CONDA_LIBSTDCXX="$CONDA_PREFIX/lib/libstdc++.so.6"
[[ -s "$CONDA_LIBSTDCXX" ]] || { echo "Missing Conda libstdc++: $CONDA_LIBSTDCXX" >&2; exit 127; }
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
if [[ "$SEMIBIN_PRELOAD_LIBSTDCXX" == 1 ]]; then
    export LD_PRELOAD="$CONDA_LIBSTDCXX${LD_PRELOAD:+:$LD_PRELOAD}"
fi

command -v SemiBin2 >/dev/null || { echo "SemiBin2 not found" >&2; exit 127; }
command -v samtools >/dev/null || { echo "samtools not found" >&2; exit 127; }

python - "$SEMIBIN_ENGINE" <<'PY'
import sys
import optree
import torch

engine = sys.argv[1]
model = torch.nn.Linear(2, 2)
torch.optim.Adam(model.parameters(), lr=1e-3)
print("optree:", optree.__version__, flush=True)
print("torch:", torch.__version__, flush=True)
print("compiled CUDA:", torch.version.cuda, flush=True)
print("CUDA available:", torch.cuda.is_available(), flush=True)
if engine == "gpu":
    if not torch.cuda.is_available():
        raise SystemExit("CUDA unavailable for SemiBin2 --engine gpu")
    print("GPU:", torch.cuda.get_device_name(0), flush=True)
PY

mapfile -t bams < <(find "$bam_dir" -maxdepth 1 -type f -name '*.sorted.bam' -print | sort)
(( ${#bams[@]} == sample_count )) || {
    echo "Expected $sample_count BAMs, found ${#bams[@]} for $group_slug" >&2
    exit 3
}
for bam in "${bams[@]}"; do
    [[ -s "$bam" && -s "$bam.bai" ]] || { echo "Missing BAM/BAI: $bam" >&2; exit 3; }
    samtools quickcheck "$bam" >/dev/null 2>&1 || { echo "Invalid BAM: $bam" >&2; exit 3; }
done

done_file="$semibin_out/.joint_semibin2.done"
if [[ -s "$done_file" && ( -d "$semibin_out/bins" || -d "$semibin_out/output_bins" ) ]]; then
    echo "[REUSE] SemiBin2 output: $semibin_out"
    exit 0
fi
[[ ! -e "$semibin_out" ]] || {
    echo "Incomplete SemiBin2 output exists; inspect or move it before retrying: $semibin_out" >&2
    exit 4
}

SemiBin2 multi_easy_bin --self-supervised -i "$reference" -b "${bams[@]}" \
    -o "$semibin_out" --separator ':' --engine "$SEMIBIN_ENGINE" \
    --threads "$THREADS" --random-seed "$SEED"

[[ -d "$semibin_out/bins" || -d "$semibin_out/output_bins" ]] || {
    echo "SemiBin2 produced no bins directory: $semibin_out" >&2
    exit 5
}
printf 'status=SUCCESS\ngroup=%s\nengine=%s\nreference=%s\nbam_count=%s\nresult=%s\ncompleted_at=%s\n' \
    "$group" "$SEMIBIN_ENGINE" "$reference" "${#bams[@]}" "$semibin_out" "$(date -Is)" > "$done_file"
