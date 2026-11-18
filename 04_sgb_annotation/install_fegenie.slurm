#!/usr/bin/env bash
#SBATCH --job-name=fegenie_install
#SBATCH --partition=sapphire
#SBATCH --qos=publiccpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=slurm-%x-%j.out
set -Eeuo pipefail
source "${ANNOT_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}/common.sh"

FEGENIE_ROOT="${FEGENIE_ROOT:-$ANN_ROOT/software/FeGenie}"
FEGENIE_ENV="${FEGENIE_ENV:-/path/to/scratch/software/anaconda3/envs/fegenie}"
INSTALL_DIR="$ANN_ROOT/software/fegenie_install"
DONE="$INSTALL_DIR/.done"

if [[ -s "$DONE" ]]; then
    echo "[REUSE] FeGenie installation"
    cat "$DONE"
    exit 0
fi

module load Anaconda3/2024.02-1
eval "$(conda shell.bash hook)"

mkdir -p "$ANN_ROOT/software" "$INSTALL_DIR"
if [[ ! -e "$FEGENIE_ROOT" ]]; then
    git clone https://github.com/Arkadiy-Garber/FeGenie.git "$FEGENIE_ROOT"
elif [[ ! -d "$FEGENIE_ROOT/.git" ]]; then
    die "Existing FeGenie path is not a Git checkout: $FEGENIE_ROOT"
fi

if [[ ! -d "$FEGENIE_ENV/conda-meta" ]]; then
    solver=conda
    command -v mamba >/dev/null && solver=mamba
    "$solver" create --prefix "$FEGENIE_ENV" \
        -c conda-forge -c bioconda --strict-channel-priority \
        python=3.10 \
        r-base r-ggplot2 r-stringi r-ggpubr r-reshape r-reshape2 \
        r-tidyverse r-argparse r-ggdendro r-pvclust \
        hmmer diamond prodigal blast metabat2 -y
fi

mkdir -p "$FEGENIE_ENV/etc/conda/activate.d"
cat > "$FEGENIE_ENV/etc/conda/activate.d/fegenie_paths.sh" <<EOF
export PATH="$FEGENIE_ROOT:\$PATH"
export rscripts="$FEGENIE_ROOT/rscripts"
export iron_hmms="$FEGENIE_ROOT/hmms/iron"
EOF

conda activate "$FEGENIE_ENV"
export PATH="$FEGENIE_ROOT:$PATH"
export rscripts="$FEGENIE_ROOT/rscripts"
export iron_hmms="$FEGENIE_ROOT/hmms/iron"

require_file "$FEGENIE_ROOT/FeGenie.py"
require_file "$iron_hmms/HMM-bitcutoffs.txt"
require_file "$iron_hmms/FeGenie-map.txt"
for command in python hmmsearch hmmstat prodigal diamond blastp; do
    command -v "$command" >/dev/null || die "Missing FeGenie dependency: $command"
done

hmm_count=$(find "$iron_hmms" -type f -name '*.hmm' | wc -l | tr -d '[:space:]')
(( hmm_count > 0 )) || die "No FeGenie HMM files found"
failed="$INSTALL_DIR/hmmstat_failed.txt"
: > "$failed"
while IFS= read -r hmm; do
    hmmstat "$hmm" >/dev/null 2>&1 || printf '%s\n' "$hmm" >> "$failed"
done < <(find "$iron_hmms" -type f -name '*.hmm' | sort)
[[ ! -s "$failed" ]] || die "FeGenie contains HMM files unreadable by this HMMER; see $failed"

python "$FEGENIE_ROOT/FeGenie.py" -h > "$INSTALL_DIR/FeGenie-help.txt"
git -C "$FEGENIE_ROOT" rev-parse HEAD > "$INSTALL_DIR/git_commit.txt"
conda list --explicit > "$INSTALL_DIR/conda-explicit.txt"
hmmsearch_help=$(hmmsearch -h 2>&1 || true)
hmmsearch_version=$(awk '/^# HMMER/ {print; exit}' <<< "$hmmsearch_help")
prodigal_help=$(prodigal -v 2>&1 || true)
prodigal_version=$(awk 'NF {print; exit}' <<< "$prodigal_help")
{
    printf 'status=SUCCESS\n'
    printf 'installed_at=%s\n' "$(date -Is)"
    printf 'git_commit=%s\n' "$(cat "$INSTALL_DIR/git_commit.txt")"
    printf 'environment=%s\n' "$FEGENIE_ENV"
    printf 'hmm_count=%s\n' "$hmm_count"
    printf 'hmmsearch=%s\n' "$hmmsearch_version"
    printf 'prodigal=%s\n' "$prodigal_version"
} > "$DONE"

cat "$DONE"

