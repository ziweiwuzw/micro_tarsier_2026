#!/usr/bin/env bash
# Delete only the temporary BAM/BAI and Strobealign .sti files for one group.
# This job is submitted only after both SemiBin2 and VAMB succeed.
# Required exports: JOINT_ROOT, GROUP_INDEX, GROUPS_TSV.
#SBATCH --job-name=joint_cleanup
#
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=slurm-%x-%A_%a.out
#SBATCH --error=slurm-%x-%A_%a.err
set -Eeuo pipefail
trap 'rc=$?; echo "[FAILED] group=${GROUP_INDEX:-?} line=${LINENO} status=${rc}" >&2; exit "$rc"' ERR

ROOT="${JOINT_ROOT:?Set JOINT_ROOT}"
GROUPS_MANIFEST="${GROUPS_TSV:-$ROOT/joint_groups.tsv}"
INDEX="${GROUP_INDEX:?Set GROUP_INDEX}"

[[ -s "$GROUPS_MANIFEST" ]] || { echo "Missing groups manifest: $GROUPS_MANIFEST" >&2; exit 2; }
row=$(awk -F '\t' -v n="$INDEX" 'NR > 1 && $1 == n {print; count++} END {if (count != 1) exit 1}' "$GROUPS_MANIFEST") || { echo "No unique group $INDEX" >&2; exit 2; }
IFS=$'\t' read -r _ group group_slug group_dir reference _ bam_dir _ _ sample_count _ <<< "$row"
[[ -d "$bam_dir" ]] || { echo "Missing BAM directory: $bam_dir" >&2; exit 2; }

find "$bam_dir" -maxdepth 1 -type f \( -name '*.bam' -o -name '*.bam.bai' -o -name '*.bai' -o -name '*.bam.done' -o -name '*.part.*' \) -delete
find "$bam_dir" -maxdepth 1 -type d -name 'sort_tmp_*' -exec rm -rf -- {} +

# Remove only indexes created by this joint workflow.  Pre-existing indexes
# are intentionally preserved so this cleanup cannot damage another pipeline.
for marker in "$group_dir"/.joint_strobealign_r*.created; do
    [[ -e "$marker" ]] || continue
    profile="${marker##*.joint_strobealign_r}"
    profile="${profile%.created}"
    rm -f -- "${reference}.r${profile}.sti" "$marker"
done

if find "$bam_dir" -maxdepth 1 -type f \( -name '*.bam' -o -name '*.bam.bai' -o -name '*.bai' \) -print -quit | grep -q .; then
    echo "BAM cleanup validation failed: $bam_dir" >&2
    exit 3
fi
printf 'status=SUCCESS\ngroup=%s\nremoved_bam_dir=%s\ncompleted_at=%s\n' "$group" "$bam_dir" "$(date -Is)" > "$group_dir/.joint_cleanup.done"
