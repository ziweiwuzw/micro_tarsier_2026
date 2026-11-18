#!/usr/bin/env bash
#SBATCH --job-name=sgb_ann_launch
#SBATCH --partition=cascade,sapphire
#SBATCH --qos=publiccpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --output=slurm-%x-%j.out
set -Eeuo pipefail
source "${ANNOT_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}/common.sh"
unset SGB_LIST

require_file "$PIPELINE_ROOT/drep/.drep.done"
require_file "$PIPELINE_ROOT/sgb/sgb_membership.tsv"
require_file "$PIPELINE_ROOT/sgb/sgb_representatives.tsv"
require_dir "$PIPELINE_ROOT/sgb/representatives"
require_dir "$PIPELINE_ROOT/hq/genomes"

# submit.sh freezes and validates every member MAG before queuing annotation jobs.
bash "$SCRIPT_DIR/submit.sh" all
