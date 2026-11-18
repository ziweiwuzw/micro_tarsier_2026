#!/usr/bin/env bash
#SBATCH --partition=cascade,sapphire
#SBATCH --qos=publiccpu
#SBATCH --output=slurm-%x-%j.out
# Resource requests are set per mode by the submitter.
set -Eeuo pipefail
source "${ANNOT_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}/common.sh"

case "${1:-}" in
    summarize)
        require_file "$PIPELINE_ROOT/tables/sample_x_sgb_tpm.tsv"
        python3 "$SCRIPT_DIR/annotation.py" --pipeline-root "$PIPELINE_ROOT" --out "$ANN_OUT" summarize
        ;;
    integrate)
        require_file "$PIPELINE_ROOT/tables/sample_x_sgb_tpm.tsv"
        python3 "$SCRIPT_DIR/integrate.py" --pipeline-root "$PIPELINE_ROOT" --out "$ANN_OUT"
        ;;
    fegenie_export)
        require_file "$ANN_OUT/input_manifest.tsv"
        python3 "$SCRIPT_DIR/export_fegenie.py" --out "$ANN_OUT"
        ;;
    *) die 'Usage: summary_job.slurm summarize|integrate|fegenie_export' ;;
esac
