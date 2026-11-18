#!/usr/bin/env bash
#SBATCH --job-name=mag_refine
#SBATCH --partition=cascade,sapphire
#SBATCH --qos=publiccpu
#SBATCH --cpus-per-task=24
#SBATCH --mem=192G
#SBATCH --time=14-00:00:00
#SBATCH --output=slurm-%x-%A_%a.out
#SBATCH --error=slurm-%x-%A_%a.err

set -Eeuo pipefail
trap 'rc=$?; echo "[FAILED] array=${SLURM_ARRAY_TASK_ID:-?} line=${LINENO} status=${rc}" >&2; exit "$rc"' ERR
HERE="${PIPELINE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
[[ -n "${PIPELINE_ROOT:-}" ]] || { echo '[ERROR] PIPELINE_ROOT is required' >&2; exit 1; }
export SOURCE_ROOT="$PIPELINE_ROOT"
# shellcheck disable=SC1091
source "$HERE/../../03_mag_sgb_kma/common.sh"

require_file "$PILOT_GROUPS_TSV"
INDEX="${SLURM_ARRAY_TASK_ID:?Submit as a Slurm array}"
THREADS="${SLURM_CPUS_PER_TASK:-24}"
MEM_GB="${METWRAP_MEM_GB:-180}"
MIN_COMPLETENESS="${METWRAP_MIN_COMPLETENESS:-50}"
MAX_CONTAMINATION="${METWRAP_MAX_CONTAMINATION:-10}"
METWRAP_QUICK="${METWRAP_QUICK:-1}"
FORCE_SINGLE_BINNER_GROUPS="${FORCE_SINGLE_BINNER_GROUPS:-Papio_ursinus:semibin2}"
CHECKM_RECURSION_LIMIT="${CHECKM_RECURSION_LIMIT:-100000}"

row=$(awk -F '\t' -v n="$INDEX" 'NR > 1 && $1 == n {print; c++} END {if (c != 1) exit 1}' "$PILOT_GROUPS_TSV") \
    || die "No unique pilot manifest row for index $INDEX"
manifest_field() {
    local column="$1"
    awk -F '\t' -v n="$INDEX" -v c="$column" \
        'NR > 1 && $1 == n {print $c; exit}' "$PILOT_GROUPS_TSV"
}
group=$(manifest_field 3)
group_slug=$(manifest_field 4)
reference=$(manifest_field 5)
semibin_input=$(manifest_field 6)
vamb_input=$(manifest_field 7)
comebin_input=$(manifest_field 8)
group_output=$(manifest_field 9)
require_file "$reference"

forced_method=""
IFS=',' read -r -a forced_rules <<< "$FORCE_SINGLE_BINNER_GROUPS"
for rule in "${forced_rules[@]}"; do
    rule="${rule//[[:space:]]/}"
    [[ -n "$rule" ]] || continue
    [[ "$rule" == *:* ]] \
        || die "Invalid FORCE_SINGLE_BINNER_GROUPS entry '$rule'; expected group:method"
    rule_group="${rule%%:*}"
    rule_method="${rule#*:}"
    if [[ "$rule_group" == "$group" || "$rule_group" == "$group_slug" ]]; then
        forced_method="${rule_method,,}"
        break
    fi
done

if [[ -n "$forced_method" ]]; then
    case "$forced_method" in
        semibin2)
            [[ -n "$semibin_input" ]] || die "Forced SemiBin2 input is unavailable for $group_slug"
            vamb_input=""
            comebin_input=""
            ;;
        vamb)
            [[ -n "$vamb_input" ]] || die "Forced VAMB input is unavailable for $group_slug"
            semibin_input=""
            comebin_input=""
            ;;
        comebin)
            [[ -n "$comebin_input" ]] || die "Forced COMEBin input is unavailable for $group_slug"
            semibin_input=""
            vamb_input=""
            ;;
        *)
            die "Unsupported forced binner '$forced_method' for $group_slug"
            ;;
    esac
    echo "[FORCE] group=$group_slug single_binner=$forced_method"
fi

bin_dirs=()
bin_labels=()
if [[ -n "$semibin_input" ]]; then
    require_dir "$semibin_input"
    bin_dirs+=("$semibin_input")
    bin_labels+=("SemiBin2")
fi
if [[ -n "$vamb_input" ]]; then
    require_dir "$vamb_input"
    bin_dirs+=("$vamb_input")
    bin_labels+=("VAMB")
fi
if [[ -n "$comebin_input" ]]; then
    require_dir "$comebin_input"
    bin_dirs+=("$comebin_input")
    bin_labels+=("COMEBin")
fi

(( ${#bin_dirs[@]} > 0 )) || die "No usable binner input directories for $group_slug"

out="${group_output}/refinement"
done_file="${out}/.metawrap_refinement.done"
final_bins="${out}/metawrap_${MIN_COMPLETENESS}_${MAX_CONTAMINATION}_bins"
if [[ -s "$done_file" && -d "$final_bins" ]]; then
    echo "[REUSE] $group_slug: $final_bins"
    exit 0
fi
[[ ! -e "$out" ]] || die "Incomplete MetaWRAP output exists: $out"
mkdir -p "${group_output}/tmp"
export TMPDIR="${group_output}/tmp"

if (( ${#bin_dirs[@]} == 1 )); then
    echo "[RUN] group=$group_slug mode=direct_single_binner method=${bin_labels[0]}"
    mkdir -p "$final_bins"
    shopt -s nullglob
    direct_bins=("${bin_dirs[0]}"/*.fa)
    (( ${#direct_bins[@]} > 0 )) || die "No staged FASTA bins in ${bin_dirs[0]}"
    for source in "${direct_bins[@]}"; do
        target="$final_bins/$(basename "$source")"
        ln -s "$source" "$target"
    done
    printf 'status=SUCCESS\nmode=direct_single_binner\nmethod=%s\nbin_count=%s\ncompleted_at=%s\n' \
        "${bin_labels[0]}" "${#direct_bins[@]}" "$(date -Is)" > "$done_file"
    echo "[SUCCESS] group=$group_slug direct_bins=${#direct_bins[@]}"
    exit 0
fi

module purge
module load GCC/11.3.0 OpenMPI/4.1.4 metaWRAP/1.3-Python-2.7.18
load_checkm1_database_module
command -v metawrap >/dev/null || die "metawrap not found"
command -v checkm >/dev/null || die "CheckM1 not found"
require_file "${CHECKM_DATA_PATH}/.dmanifest"
[[ "$CHECKM_RECURSION_LIMIT" =~ ^[1-9][0-9]*$ ]] \
    || die "CHECKM_RECURSION_LIMIT must be a positive integer"
checkm_python_site="$HERE/checkm_python_site"
require_file "$checkm_python_site/sitecustomize.py"
export CHECKM_RECURSION_LIMIT
export PYTHONPATH="$checkm_python_site${PYTHONPATH:+:$PYTHONPATH}"
echo "[CHECKM1] executable=$(command -v checkm) database=$CHECKM_DATA_PATH"
echo "[CHECKM1] python_recursion_limit=$CHECKM_RECURSION_LIMIT"

args=(bin_refinement -o "$out" -t "$THREADS" -m "$MEM_GB"
    -c "$MIN_COMPLETENESS" -x "$MAX_CONTAMINATION")
input_flag=A
for bin_dir in "${bin_dirs[@]}"; do
    args+=("-${input_flag}" "$bin_dir")
    case "$input_flag" in
        A) input_flag=B ;;
        B) input_flag=C ;;
        *) die "Too many binner input directories for MetaWRAP" ;;
    esac
done
[[ "$METWRAP_QUICK" == 1 ]] && args+=(--quick)

echo "[RUN] group=$group_slug mode=metawrap_refinement methods=${bin_labels[*]}"
metawrap "${args[@]}"

[[ -d "$final_bins" ]] || die "Expected final bins directory missing: $final_bins"
bin_count=$(find "$final_bins" -maxdepth 1 -type f -name '*.fa' | wc -l | tr -d '[:space:]')
(( bin_count > 0 )) || die "MetaWRAP returned zero final bins for $group_slug"
printf 'status=SUCCESS\ngroup=%s\nreference=%s\nbin_count=%s\ncompleted_at=%s\n' \
    "$group" "$reference" "$bin_count" "$(date -Is)" > "$done_file"
echo "[SUCCESS] group=$group_slug refined_bins=$bin_count"
