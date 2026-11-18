#!/usr/bin/env bash
# Submit one CPU dependency DAG per group:
#   Strobealign BAM array -> {SemiBin2 CPU, VAMB CPU} -> cleanup
# The next group's mapping can be released by the first successful binner while
# cleanup still waits for both binners to succeed.
set -Eeuo pipefail

ROOT="${JOINT_ROOT:?Set JOINT_ROOT to the existing SemiBin2 by-group root}"
GROUPS_MANIFEST="${GROUPS_TSV:-$ROOT/joint_groups.tsv}"
HERE="$(cd "$(dirname "$0")" && pwd)"
MAP_SCRIPT="${MAP_SCRIPT:-$HERE/map_host_group_reads.slurm}"
SEMIBIN_SCRIPT="${SEMIBIN_SCRIPT:-$HERE/run_semibin2.slurm}"
VAMB_SCRIPT="${VAMB_SCRIPT:-$HERE/run_vamb.slurm}"
CLEAN_SCRIPT="${CLEAN_SCRIPT:-$HERE/cleanup_joint_bams.slurm}"

MAP_CONCURRENCY="${MAP_CONCURRENCY:-4}"
MAP_CPUS="${MAP_CPUS:-32}"
MAP_MEM="${MAP_MEM:-128G}"
MAP_TIME="${MAP_TIME:-2-00:00:00}"
SEMIBIN_CPUS="${SEMIBIN_CPUS:-32}"
SEMIBIN_MEM="${SEMIBIN_MEM:-128G}"
SEMIBIN_TIME="${SEMIBIN_TIME:-7-00:00:00}"
VAMB_CPUS="${VAMB_CPUS:-32}"
VAMB_MEM="${VAMB_MEM:-128G}"
VAMB_TIME="${VAMB_TIME:-2-00:00:00}"

MAP_ROUTE="${MAP_ROUTE:-auto}"       # auto, fos, or public
BIN_ROUTE="${BIN_ROUTE:-fos}"        # fos or public
FOS_PARTITION="${FOS_PARTITION:-fos}"
FOS_QOS="${FOS_QOS:-fos}"
FOS_CPU_LIMIT="${FOS_CPU_LIMIT:-128}"
LOG_ROOT="${LOG_ROOT:-/path/to/project/08_vamb_semibin2/logs}"
NEXT_GROUP_AFTER="${NEXT_GROUP_AFTER:-first}"  # first, semibin2, vamb, or cleanup
START_GROUP_INDEX="${START_GROUP_INDEX:-1}"
INITIAL_DEPENDENCY="${INITIAL_DEPENDENCY:-}"
MAX_NEW_JOBS="${MAX_NEW_JOBS:-0}"  # 0 means no submit-count cap

[[ -s "$GROUPS_MANIFEST" ]] || { echo "Missing groups manifest: $GROUPS_MANIFEST" >&2; exit 2; }
[[ "$START_GROUP_INDEX" =~ ^[0-9]+$ && "$START_GROUP_INDEX" -gt 0 ]] || {
    echo "START_GROUP_INDEX must be a positive integer" >&2; exit 2;
}
[[ "$MAX_NEW_JOBS" =~ ^[0-9]+$ ]] || {
    echo "MAX_NEW_JOBS must be a non-negative integer" >&2; exit 2;
}
case "$NEXT_GROUP_AFTER" in
    first|semibin2|vamb|cleanup) ;;
    *) echo "NEXT_GROUP_AFTER must be first, semibin2, vamb, or cleanup" >&2; exit 2 ;;
esac
mkdir -p "$LOG_ROOT"

set_route_args() {
    local route="$1"
    ROUTE_ARGS=()
    case "$route" in
        fos) ROUTE_ARGS=(--partition="$FOS_PARTITION" --qos="$FOS_QOS") ;;
        public) ;; # Omit partition and QoS: use cluster defaults.
        *) echo "Route must be fos or public, got: $route" >&2; exit 2 ;;
    esac
}

effective_map_route="$MAP_ROUTE"
if [[ "$MAP_ROUTE" == auto ]]; then
    # Route oversized aggregate mapping requests away from fos when auto is used.
    if (( MAP_CONCURRENCY * MAP_CPUS >= FOS_CPU_LIMIT )); then
        effective_map_route=public
    else
        effective_map_route=fos
    fi
fi
[[ "$effective_map_route" == fos || "$effective_map_route" == public ]] || {
    echo "MAP_ROUTE must be auto, fos, or public" >&2; exit 2;
}

previous_dependency="$INITIAL_DEPENDENCY"
submitted_groups=0
submitted_job_slots=0
while IFS=$'\t' read -r index group group_slug group_dir reference mapping_jobs bam_dir semibin_out vamb_out sample_count reference_contigs; do
    [[ "$index" == "array_index" ]] && continue
    [[ -n "$index" && -s "$mapping_jobs" ]] || { echo "Invalid group row: $index $group_slug" >&2; exit 2; }
    (( sample_count > 0 )) || { echo "Group has no samples: $group_slug" >&2; exit 2; }
    (( index < START_GROUP_INDEX )) && continue

    # Slurm counts every array element toward MaxSubmitJobs. Each group adds
    # sample_count mapping elements plus SemiBin2, VAMB, and cleanup jobs.
    group_job_slots=$((sample_count + 3))
    if (( MAX_NEW_JOBS > 0 && submitted_job_slots + group_job_slots > MAX_NEW_JOBS )); then
        printf '[PAUSE] submit cap reached: submitted_groups=%s submitted_jobs=%s next_group_index=%s next_group=%s next_group_jobs=%s\n' \
            "$submitted_groups" "$submitted_job_slots" "$index" "$group_slug" "$group_job_slots"
        printf '[RESUME] START_GROUP_INDEX=%s INITIAL_DEPENDENCY=%q\n' "$index" "$previous_dependency"
        break
    fi

    dependency=()
    [[ -n "$previous_dependency" ]] && dependency=(--dependency="$previous_dependency")
    export_args="ALL,JOINT_ROOT=${ROOT},GROUPS_TSV=${GROUPS_MANIFEST},GROUP_INDEX=${index},GROUP_JOBS=${mapping_jobs}"
    set_route_args "$effective_map_route"
    map_job=$(sbatch --parsable "${dependency[@]}" "${ROUTE_ARGS[@]}" \
        --cpus-per-task="$MAP_CPUS" --mem="$MAP_MEM" --time="$MAP_TIME" \
        --array="1-${sample_count}%${MAP_CONCURRENCY}" --export="$export_args" "$MAP_SCRIPT")

    common_export="ALL,JOINT_ROOT=${ROOT},GROUPS_TSV=${GROUPS_MANIFEST},GROUP_INDEX=${index},SEMIBIN_ENGINE=cpu"
    set_route_args "$BIN_ROUTE"
    semi_job=$(sbatch --parsable "${ROUTE_ARGS[@]}" --dependency="afterok:${map_job}" \
        --cpus-per-task="$SEMIBIN_CPUS" --mem="$SEMIBIN_MEM" --time="$SEMIBIN_TIME" \
        --export="$common_export" "$SEMIBIN_SCRIPT")
    vamb_job=$(sbatch --parsable "${ROUTE_ARGS[@]}" --dependency="afterok:${map_job}" \
        --cpus-per-task="$VAMB_CPUS" --mem="$VAMB_MEM" --time="$VAMB_TIME" \
        --export="$common_export" "$VAMB_SCRIPT")
    cleanup_job=$(sbatch --parsable "${ROUTE_ARGS[@]}" --dependency="afterok:${semi_job}:${vamb_job}" \
        --cpus-per-task=1 --mem=2G --time=00:30:00 \
        --export="$common_export" "$CLEAN_SCRIPT")

    case "$NEXT_GROUP_AFTER" in
        first) previous_dependency="afterok:${semi_job}?afterok:${vamb_job}" ;;
        semibin2) previous_dependency="afterok:${semi_job}" ;;
        vamb) previous_dependency="afterok:${vamb_job}" ;;
        cleanup) previous_dependency="afterok:${cleanup_job}" ;;
    esac
    submitted_groups=$((submitted_groups + 1))
    submitted_job_slots=$((submitted_job_slots + group_job_slots))
    printf 'group=%s\tsamples=%s\tmap_route=%s\tbin_route=%s\tmap=%s\tsemibin2=%s\tvamb=%s\tcleanup=%s\n' \
        "$group_slug" "$sample_count" "$effective_map_route" "$BIN_ROUTE" "$map_job" "$semi_job" "$vamb_job" "$cleanup_job"
done < "$GROUPS_MANIFEST"
