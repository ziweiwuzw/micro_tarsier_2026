#!/usr/bin/env bash
# Generate sorted/indexed Strobealign BAMs for one already-prepared group.
# Required exports: JOINT_ROOT, GROUP_JOBS, SLURM_ARRAY_TASK_ID.
# Optional: STROBEALIGN_ENV, THREADS, SAMTOOLS_THREADS, STROBEALIGN_R,
#           STROBEALIGN_THREADS, BAM_FILTER_FLAGS, CREATE_INDEX, SORT_TMPDIR.
#SBATCH --job-name=joint_map
#
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00
#SBATCH --output=slurm-%x-%A_%a.out
#SBATCH --error=slurm-%x-%A_%a.err
set -Eeuo pipefail
trap 'rc=$?; echo "[FAILED] task=${SLURM_ARRAY_TASK_ID:-?} line=${LINENO} status=${rc}" >&2; exit "$rc"' ERR

ROOT="${JOINT_ROOT:?Set JOINT_ROOT to the SemiBin2 by-group root}"
JOBS="${GROUP_JOBS:?Set GROUP_JOBS to one group joint_mapping_jobs.tsv}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-32}}"
SAMTOOLS_THREADS="${SAMTOOLS_THREADS:-4}"
STROBEALIGN_THREADS="${STROBEALIGN_THREADS:-$((THREADS - 2 * SAMTOOLS_THREADS))}"
STROBEALIGN_R="${STROBEALIGN_R:-6}"
BAM_FILTER_FLAGS="${BAM_FILTER_FLAGS:-4}"
CREATE_INDEX="${CREATE_INDEX:-1}"
STROBEALIGN_ENV="${STROBEALIGN_ENV:-/path/to/home/software/anaconda3/envs/Taxvamb}"
SAMTOOLS_ENV="${SAMTOOLS_ENV:-/path/to/scratch/software/anaconda3/envs/SemiBin2}"
STROBEALIGN_BIN="${STROBEALIGN_BIN:-$STROBEALIGN_ENV/bin/strobealign}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-$SAMTOOLS_ENV/bin/samtools}"
TASK="${SLURM_ARRAY_TASK_ID:?Submit this script as a Slurm array}"

(( STROBEALIGN_THREADS > 0 )) || { echo "THREADS is too small for SAMTOOLS_THREADS=$SAMTOOLS_THREADS" >&2; exit 2; }

[[ -s "$JOBS" ]] || { echo "Missing mapping jobs: $JOBS" >&2; exit 2; }
row=$(awk -F '\t' -v n="$TASK" 'NR > 1 && $1 == n {print; count++} END {if (count != 1) exit 1}' "$JOBS") || { echo "No unique mapping job $TASK" >&2; exit 2; }
IFS=$'\t' read -r _ sample group _ group_dir reference r1 r2 profile bam <<< "$row"
[[ -s "$reference" && -s "$r1" && -s "$r2" ]] || { echo "Missing input for $sample" >&2; exit 2; }
[[ "$profile" =~ ^[0-9]+$ && "$profile" -gt 0 ]] || { echo "Invalid canonical_profile for $sample: $profile" >&2; exit 2; }
mkdir -p "$(dirname "$bam")"

export PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME PIP_USER || true
[[ -x "$STROBEALIGN_BIN" ]] || { echo "Strobealign executable not found: $STROBEALIGN_BIN" >&2; exit 127; }
[[ -x "$SAMTOOLS_BIN" ]] || { echo "samtools executable not found: $SAMTOOLS_BIN" >&2; exit 127; }

if [[ -s "$bam" && -s "$bam.bai" ]] && "$SAMTOOLS_BIN" quickcheck "$bam" >/dev/null 2>&1; then
    echo "[REUSE] $bam"
    exit 0
fi

reference_index="${reference}.r${profile}.sti"
index_ready() { [[ -s "$reference_index" ]]; }
if [[ "$CREATE_INDEX" == 1 ]] && ! index_ready; then
    lock="${reference_index}.lock"
    if mkdir "$lock" 2>/dev/null; then
        trap 'rmdir "$lock" 2>/dev/null || true' EXIT
        if ! index_ready; then
            "$STROBEALIGN_BIN" --create-index -t "$STROBEALIGN_THREADS" -r "$profile" "$reference"
            touch "$group_dir/.joint_strobealign_r${profile}.created"
        fi
        rmdir "$lock"
        trap - EXIT
    else
        for _ in $(seq 1 720); do
            index_ready && break
            [[ -d "$lock" ]] || break
            sleep 10
        done
        index_ready || { echo "Timed out waiting for Strobealign index: $reference_index" >&2; exit 3; }
    fi
fi

bam_part="${bam}.part.${SLURM_JOB_ID:-$$}.${TASK}"
sort_prefix="${SORT_TMPDIR:-$(dirname "$bam")/sort_tmp_${SLURM_JOB_ID:-$$}_${TASK}}"
mkdir -p "$sort_prefix"
cleanup_tmp() { rm -f -- "$bam_part" "$bam_part.bai"; rm -rf -- "$sort_prefix"; }
trap cleanup_tmp EXIT

cmd=("$STROBEALIGN_BIN" -t "$STROBEALIGN_THREADS" -r "$profile" -R "$STROBEALIGN_R")
if [[ "$CREATE_INDEX" == 1 ]] && index_ready; then
    cmd+=(--use-index)
fi
cmd+=("$reference" "$r1" "$r2")
"${cmd[@]}" |
    "$SAMTOOLS_BIN" view -@ "$SAMTOOLS_THREADS" -b -F "$BAM_FILTER_FLAGS" - |
    "$SAMTOOLS_BIN" sort -@ "$SAMTOOLS_THREADS" -T "$sort_prefix/sort" -o "$bam_part" -
"$SAMTOOLS_BIN" quickcheck -v "$bam_part"
"$SAMTOOLS_BIN" index -@ "$SAMTOOLS_THREADS" "$bam_part"
mv "$bam_part" "$bam"
mv "$bam_part.bai" "$bam.bai"
"$SAMTOOLS_BIN" quickcheck "$bam"
printf 'status=SUCCESS\nmapper=strobealign\nsample=%s\nreference=%s\nprofile=%s\nbam=%s\ncompleted_at=%s\n' \
    "$sample" "$reference" "$profile" "$bam" "$(date -Is)" > "${bam}.done"
