#!/bin/bash
#SBATCH --job-name=ccmetagen
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=ccmetagen_%j.out
#SBATCH --error=ccmetagen_%j.err

set -euo pipefail

PPN=${SLURM_CPUS_PER_TASK:-2}

# ==============================
# Paths
# ==============================
PROJECT_DIR="${PROJECT_DIR:?Set PROJECT_DIR to the focal-read project directory}"
WORK_DIR=${PROJECT_DIR}/special_downloads/macrogen_fastp_ccmetagen

KMA_DIR=${WORK_DIR}/05_kma_nt
CCM_DIR=${WORK_DIR}/06_ccmetagen_nt
LOG_DIR=${WORK_DIR}/logs

SAMPLE_LIST=${WORK_DIR}/sample_list.txt

mkdir -p "${CCM_DIR}"
mkdir -p "${LOG_DIR}"

# ==============================
# Fixed conda environment paths
# Do NOT use bare python here.
# ==============================
ENV="${CCMETAGEN_ENV:?Set CCMETAGEN_ENV to the CCMetagen environment}"
PY=${ENV}/bin/python
CCM=${ENV}/bin/CCMetagen.py
CCMM=${ENV}/bin/CCMetagen_merge.py

# Avoid user-level Python packages, e.g. ~/.local/lib/python*/site-packages
export PYTHONNOUSERSITE=1
unset PYTHONPATH

# Make Krona and other env binaries visible to CCMetagen
export PATH=${ENV}/bin:${PATH}
hash -r

echo "======================================"
echo "Start CCMetagen: $(date)"
echo "Project dir: ${PROJECT_DIR}"
echo "Work dir: ${WORK_DIR}"
echo "KMA dir: ${KMA_DIR}"
echo "CCMetagen output dir: ${CCM_DIR}"
echo "Threads: ${PPN}"
echo "Environment: ${ENV}"
echo "Python: ${PY}"
echo "======================================"

echo "Checking tools:"
"${PY}" -s "${CCM}" --version || true
echo "CCMetagen path: ${CCM}"
echo "CCMetagen_merge path: ${CCMM}"

echo "Checking pandas/numpy:"
"${PY}" -s - <<'PYCODE'
import sys
import pandas as pd
import numpy as np
print("python exe:", sys.executable)
print("pandas:", pd.__version__, pd.__file__)
print("numpy:", np.__version__, np.__file__)
PYCODE

echo "Checking Krona:"
which ktImportText || true
which ktImportTaxonomy || true

# ==============================
# Check input files
# ==============================
if [[ ! -s "${SAMPLE_LIST}" ]]; then
    echo "ERROR: sample list not found: ${SAMPLE_LIST}"
    echo "Please run process_focal_reads_fastp.slurm first."
    exit 1
fi

echo "Detected samples:"
cat "${SAMPLE_LIST}"

echo "Checking KMA result files..."
while read -r sample
do
    [[ -z "${sample}" ]] && continue

    RES=${KMA_DIR}/${sample}.res
    MAPSTAT=${KMA_DIR}/${sample}.mapstat

    if [[ ! -s "${RES}" ]]; then
        echo "ERROR: Missing KMA .res file: ${RES}"
        exit 1
    fi

    if [[ ! -s "${MAPSTAT}" ]]; then
        echo "ERROR: Missing KMA .mapstat file: ${MAPSTAT}"
        echo "KMA needs to be run with -ef to generate .mapstat."
        exit 1
    fi

    echo "OK: ${sample}.res and ${sample}.mapstat found"
done < "${SAMPLE_LIST}"

# ==============================
# Run CCMetagen
# ==============================
echo "======================================"
echo "Running CCMetagen..."
echo "Using: -ef y -du fr -map MAPSTAT"
echo "Important: removing extra KMA '## command' line from mapstat before CCMetagen"
echo "======================================"

while read -r sample
do
    [[ -z "${sample}" ]] && continue

    echo "--------------------------------------"
    echo "CCMetagen sample: ${sample}"
    echo "Start: $(date)"
    echo "--------------------------------------"

    RES=${KMA_DIR}/${sample}.res
    MAPSTAT=${KMA_DIR}/${sample}.mapstat
    MAPSTAT_CLEAN=${CCM_DIR}/${sample}.ccmetagen_compatible.mapstat
    OUT_PREFIX=${CCM_DIR}/${sample}

    # CCMetagen v1.5.0 expects the mapstat table header at line 7: header=6.
    # Newer KMA mapstat has an extra "## command" line before the real header.
    # Remove only that line so "# refSequence ..." becomes line 7.
    awk '!(NR==7 && /^## command/)' "${MAPSTAT}" > "${MAPSTAT_CLEAN}"

    echo "Original mapstat first 8 lines for ${sample}:"
    head -n 8 "${MAPSTAT}" || true

    echo "CCMetagen-compatible mapstat first 8 lines for ${sample}:"
    head -n 8 "${MAPSTAT_CLEAN}" || true

    "${PY}" -s "${CCM}" \
        -i "${RES}" \
        -o "${OUT_PREFIX}" \
        -ef y \
        -du fr \
        -map "${MAPSTAT_CLEAN}" \
        > "${LOG_DIR}/${sample}.ccmetagen.log" 2>&1

    echo "Finished CCMetagen sample: ${sample}"
    echo "Output prefix: ${OUT_PREFIX}"
    echo "End: $(date)"

done < "${SAMPLE_LIST}"

echo "======================================"
echo "CCMetagen per-sample outputs:"
ls -lh "${CCM_DIR}"
echo "======================================"

# ==============================
# Merge CCMetagen results
# ==============================
echo "======================================"
echo "Merging CCMetagen tables..."
echo "======================================"

"${PY}" -s "${CCMM}" \
    -i "${CCM_DIR}" \
    -o "${CCM_DIR}/CCM_species_table" \
    -t Species \
    > "${LOG_DIR}/CCM_species_table.merge.log" 2>&1

"${PY}" -s "${CCMM}" \
    -i "${CCM_DIR}" \
    -o "${CCM_DIR}/CCM_genus_table" \
    -t Genus \
    > "${LOG_DIR}/CCM_genus_table.merge.log" 2>&1

"${PY}" -s "${CCMM}" \
    -i "${CCM_DIR}" \
    -o "${CCM_DIR}/CCM_family_table" \
    -t Family \
    > "${LOG_DIR}/CCM_family_table.merge.log" 2>&1

"${PY}" -s "${CCMM}" \
    -i "${CCM_DIR}" \
    -o "${CCM_DIR}/CCM_closest_match_table" \
    -t Closest_match \
    > "${LOG_DIR}/CCM_closest_match_table.merge.log" 2>&1

echo "======================================"
echo "CCMetagen finished: $(date)"
echo "Main output dir:"
echo "${CCM_DIR}"
echo ""
echo "Merged tables:"
ls -lh "${CCM_DIR}"/CCM_*table* || true
echo "======================================"

sstat -j "${SLURM_JOB_ID}.batch" --format=JobID,MaxVMSize,MaxRSS,AveRSS,MinCPU 2>/dev/null || true
