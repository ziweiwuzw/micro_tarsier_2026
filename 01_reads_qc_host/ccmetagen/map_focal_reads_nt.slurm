#!/bin/bash
#SBATCH --job-name=kma_nt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=7-00:00:00
#SBATCH --output=kma_nt_%j.out
#SBATCH --error=kma_nt_%j.err

set -euo pipefail

PPN=${SLURM_CPUS_PER_TASK:-48}

# ==============================
# Paths
# ==============================
PROJECT_DIR="${PROJECT_DIR:?Set PROJECT_DIR to the focal-read project directory}"
WORK_DIR=${PROJECT_DIR}/special_downloads/macrogen_fastp_ccmetagen

CLEAN_DIR=${WORK_DIR}/01_clean_reads
KMA_DIR=${WORK_DIR}/05_kma_nt
LOG_DIR=${WORK_DIR}/logs

SAMPLE_LIST=${WORK_DIR}/sample_list.txt

# KMA indexed nt database prefix
DB_PREFIX="${DB_PREFIX:?Set DB_PREFIX to the indexed NCBI nt database}"

mkdir -p ${KMA_DIR}
mkdir -p ${LOG_DIR}

# ==============================
# Conda environment
# ==============================
module load Anaconda3
source activate "${CCMETAGEN_ENV:?Set CCMETAGEN_ENV to the CCMetagen environment}"

echo "======================================"
echo "Start KMA nt mapping: $(date)"
echo "Project dir: ${PROJECT_DIR}"
echo "Work dir: ${WORK_DIR}"
echo "Clean reads dir: ${CLEAN_DIR}"
echo "KMA output dir: ${KMA_DIR}"
echo "Database prefix: ${DB_PREFIX}"
echo "Threads: ${PPN}"
echo "======================================"

echo "Checking tools:"
which kma
kma 2>&1 | head -n 5 || true

# ==============================
# Check input files
# ==============================
if [[ ! -s ${SAMPLE_LIST} ]]; then
    echo "ERROR: sample list not found: ${SAMPLE_LIST}"
    echo "Please run process_focal_reads_fastp.slurm first."
    exit 1
fi

echo "Checking KMA database files..."
for suffix in comp.b length.b name seq.b
do
    if [[ ! -s ${DB_PREFIX}.${suffix} ]]; then
        echo "ERROR: Missing database file: ${DB_PREFIX}.${suffix}"
        exit 1
    fi
done
echo "KMA database files found."

echo "Detected samples:"
cat ${SAMPLE_LIST}

echo "Checking clean paired-end files..."
while read sample
do
    R1=${CLEAN_DIR}/${sample}_1.clean.fastq.gz
    R2=${CLEAN_DIR}/${sample}_2.clean.fastq.gz

    if [[ ! -s ${R1} ]]; then
        echo "ERROR: Missing clean R1: ${R1}"
        exit 1
    fi

    if [[ ! -s ${R2} ]]; then
        echo "ERROR: Missing clean R2: ${R2}"
        exit 1
    fi

    echo "OK: ${sample}"
done < ${SAMPLE_LIST}

# ==============================
# KMA mapping
# ==============================
echo "======================================"
echo "Running KMA..."
echo "======================================"

# Vanessa's reference script used -apm p.
# If you want to follow the CCMetagen README example strictly, change APM=p to APM=f.
APM=p

while read sample
do
    echo "--------------------------------------"
    echo "KMA sample: ${sample}"
    echo "Start: $(date)"
    echo "--------------------------------------"

    R1=${CLEAN_DIR}/${sample}_1.clean.fastq.gz
    R2=${CLEAN_DIR}/${sample}_2.clean.fastq.gz
    OUT_PREFIX=${KMA_DIR}/${sample}

    kma \
        -ipe ${R1} ${R2} \
        -o ${OUT_PREFIX} \
        -t_db ${DB_PREFIX} \
        -t ${PPN} \
        -1t1 \
        -mem_mode \
        -and \
        -apm ${APM} \
        -ef \
        > ${LOG_DIR}/${sample}.kma.log 2>&1

    echo "Finished KMA sample: ${sample}"
    echo "Output prefix: ${OUT_PREFIX}"
    echo "End: $(date)"

done < ${SAMPLE_LIST}

echo "======================================"
echo "KMA finished: $(date)"
echo "KMA output files:"
ls -lh ${KMA_DIR}
echo "======================================"

# Optional SLURM memory report
sstat -j ${SLURM_JOB_ID}.batch --format=JobID,MaxVMSize,MaxRSS,AveRSS,MinCPU 2>/dev/null || true
