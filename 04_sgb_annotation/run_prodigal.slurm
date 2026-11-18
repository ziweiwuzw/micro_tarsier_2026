#!/usr/bin/env bash
#SBATCH --job-name=sgb_orf
#SBATCH --partition=cascade,sapphire
#SBATCH --qos=publiccpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=slurm-%x-%A_%a.out
set -Eeuo pipefail
source "${ANNOT_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}/common.sh"
require_file "$ANN_OUT/input_manifest.tsv"
module load GCCcore/11.3.0
module load prodigal/2.6.3
command -v prodigal >/dev/null || die 'Prodigal module did not provide prodigal'
mag="$(mag_at_task)"
python3 "$SCRIPT_DIR/annotation.py" --pipeline-root "$PIPELINE_ROOT" --out "$ANN_OUT" genes --mag "$mag"
