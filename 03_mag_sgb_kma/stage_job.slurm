#!/usr/bin/env bash
#SBATCH --partition=cascade,sapphire
#SBATCH --qos=publiccpu
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err
# Resource requests are set per mode by submit.sh.
set -Eeuo pipefail
HERE="${PIPELINE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "$HERE/common.sh"

case "${1:-}" in
    collect)
        require_file "$PILOT_GROUPS_TSV"
        python3 "$HERE/collect_candidates.py" \
            --groups "$PILOT_GROUPS_TSV" \
            --source-root "$SOURCE_ROOT" \
            --output-root "$(dirname "$CANDIDATE_DIR")" \
            --completeness "${METWRAP_MIN_COMPLETENESS:-50}" \
            --contamination "${METWRAP_MAX_CONTAMINATION:-10}" \
            --pan-source "${PAN_SOURCE:-semibin2}"
        python3 "$HERE/prepare_reads.py" \
            --input "$INPUT_READS_MANIFEST" --output "$READS_MANIFEST"
        ;;
    filter)
        require_file "$CANDIDATE_MANIFEST"
        require_file "${CHECKM2_OUT}/quality_report.tsv"
        python3 "$HERE/filter_hq_mags.py" \
            --candidates "$CANDIDATE_MANIFEST" \
            --checkm2 "${CHECKM2_OUT}/quality_report.tsv" \
            --output-root "$(dirname "$HQ_DIR")" \
            --min-completeness "${MIN_COMPLETENESS:-90}" \
            --max-contamination "${MAX_CONTAMINATION:-5}"
        ;;
    sgb)
        require_file "${DREP_OUT}/.drep.done"
        require_file "$HQ_MANIFEST"
        python3 "$HERE/assign_sgb_ids.py" \
            --drep "$DREP_OUT" --hq-manifest "$HQ_MANIFEST" --output-root "$SGB_ROOT"
        ;;
    reference)
        require_file "$SGB_ROOT/sgb_representatives.tsv"
        python3 "$HERE/build_sgb_reference.py" \
            --representatives "$SGB_ROOT/sgb_representatives.tsv" \
            --reference-root "$REFERENCE_ROOT"
        ;;
    source_metadata)
        require_file "$SGB_ROOT/sgb_membership.tsv"
        require_file "$SGB_ROOT/sgb_representatives.tsv"
        require_file "$READS_MANIFEST"
        python3 "$HERE/build_sgb_source_tables.py" --root "$PIPELINE_ROOT" --out "$PIPELINE_ROOT/metadata"
        require_file "$PIPELINE_ROOT/metadata/sgb_summary.tsv"
        require_file "$PIPELINE_ROOT/metadata/mag_source_metadata.tsv"
        ;;
    *) die 'Usage: stage_job.slurm collect|filter|sgb|reference|source_metadata' ;;
esac
