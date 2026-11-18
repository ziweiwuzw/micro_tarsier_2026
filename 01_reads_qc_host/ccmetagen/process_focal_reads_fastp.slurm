#!/bin/bash
#SBATCH --job-name=fastp_clean_qc
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=2-04:00:00
#SBATCH --output=fastp_clean_qc_%j.out
#SBATCH --error=fastp_clean_qc_%j.err

set -euo pipefail

PPN=${SLURM_CPUS_PER_TASK:-8}

# ==============================
# Project paths
# ==============================
PROJECT_DIR="${PROJECT_DIR:?Set PROJECT_DIR to the focal-read project directory}"
RAW_DIR=${PROJECT_DIR}/special_downloads/macrogen_raw

WORK_DIR=${PROJECT_DIR}/special_downloads/macrogen_fastp_ccmetagen
CLEAN_DIR=${WORK_DIR}/01_clean_reads
FASTP_REPORT_DIR=${WORK_DIR}/02_fastp_reports
FASTQC_DIR=${WORK_DIR}/03_fastqc_clean
MULTIQC_DIR=${WORK_DIR}/04_multiqc_clean
LOG_DIR=${WORK_DIR}/logs

mkdir -p ${CLEAN_DIR}
mkdir -p ${FASTP_REPORT_DIR}
mkdir -p ${FASTQC_DIR}
mkdir -p ${MULTIQC_DIR}
mkdir -p ${LOG_DIR}

# ==============================
# Conda environments
# ==============================
module load Anaconda3

RNA_ENV="${RNA_ENV:?Set RNA_ENV to the fastp environment}"
MULTIQC_ENV="${MULTIQC_ENV:?Set MULTIQC_ENV to the MultiQC environment}"

echo "======================================"
echo "Start fastp clean + QC: $(date)"
echo "Project dir: ${PROJECT_DIR}"
echo "Raw data dir: ${RAW_DIR}"
echo "Work dir: ${WORK_DIR}"
echo "Threads: ${PPN}"
echo "======================================"

# ==============================
# Activate rna env for fastp and FastQC
# ==============================
source activate ${RNA_ENV}

echo "Conda env:"
which python || true

echo "Checking fastp:"
which fastp
fastp --version || true

echo "Checking fastqc:"
which fastqc
fastqc --version || true

# ==============================
# Check raw files
# ==============================
cd ${RAW_DIR}

echo "======================================"
echo "Input fastq.gz files:"
echo "======================================"
ls -lh *.fastq.gz

ls *.fastq.gz | sed -E 's/_[12]\.fastq\.gz$//' | sort -u > ${WORK_DIR}/sample_list.txt

echo "======================================"
echo "Detected samples:"
echo "======================================"
cat ${WORK_DIR}/sample_list.txt

echo "======================================"
echo "Checking paired-end files..."
echo "======================================"

while read sample
do
    if [[ ! -s ${RAW_DIR}/${sample}_1.fastq.gz ]]; then
        echo "ERROR: Missing ${RAW_DIR}/${sample}_1.fastq.gz"
        exit 1
    fi

    if [[ ! -s ${RAW_DIR}/${sample}_2.fastq.gz ]]; then
        echo "ERROR: Missing ${RAW_DIR}/${sample}_2.fastq.gz"
        exit 1
    fi

    echo "OK: ${sample}_1.fastq.gz and ${sample}_2.fastq.gz found"
done < ${WORK_DIR}/sample_list.txt

# ==============================
# Optional MD5 check
# ==============================
echo "======================================"
echo "Optional MD5 check..."
echo "======================================"

if [[ -s ${RAW_DIR}/md5.txt ]]; then
    echo "md5.txt found. Running md5sum check..."
    cd ${RAW_DIR}
    md5sum -c md5.txt > ${LOG_DIR}/md5_check.txt 2>&1 || {
        echo "WARNING: MD5 check failed. Please check:"
        echo "${LOG_DIR}/md5_check.txt"
        cat ${LOG_DIR}/md5_check.txt
        exit 1
    }
    echo "MD5 check passed."
else
    echo "No md5.txt found. Skip MD5 check."
fi

# ==============================
# Run fastp
# ==============================
echo "======================================"
echo "Running fastp..."
echo "======================================"

while read sample
do
    echo "Processing sample: ${sample}"

    R1=${RAW_DIR}/${sample}_1.fastq.gz
    R2=${RAW_DIR}/${sample}_2.fastq.gz

    OUT_R1=${CLEAN_DIR}/${sample}_1.clean.fastq.gz
    OUT_R2=${CLEAN_DIR}/${sample}_2.clean.fastq.gz

    JSON=${FASTP_REPORT_DIR}/${sample}.fastp.json
    HTML=${FASTP_REPORT_DIR}/${sample}.fastp.html
    LOG=${LOG_DIR}/${sample}.fastp.log

    fastp \
        --in1 ${R1} \
        --in2 ${R2} \
        --out1 ${OUT_R1} \
        --out2 ${OUT_R2} \
        --detect_adapter_for_pe \
        --qualified_quality_phred 25 \
        --unqualified_percent_limit 40 \
        --n_base_limit 3 \
        --length_required 36 \
        --thread ${PPN} \
        --html ${HTML} \
        --json ${JSON} \
        --report_title "${sample} fastp report" \
        > ${LOG} 2>&1

    echo "Finished fastp: ${sample}"
    echo "Clean R1: ${OUT_R1}"
    echo "Clean R2: ${OUT_R2}"

done < ${WORK_DIR}/sample_list.txt

echo "======================================"
echo "Clean reads:"
echo "======================================"
ls -lh ${CLEAN_DIR}/*.clean.fastq.gz

# ==============================
# Run FastQC on clean reads
# ==============================
echo "======================================"
echo "Running FastQC on clean reads..."
echo "======================================"

fastqc \
    --noextract \
    --threads ${PPN} \
    --outdir ${FASTQC_DIR} \
    ${CLEAN_DIR}/*.clean.fastq.gz

# ==============================
# Run MultiQC
# ==============================
echo "======================================"
echo "Running MultiQC..."
echo "======================================"

source activate ${MULTIQC_ENV}

echo "Checking multiqc:"
which multiqc
multiqc --version || true

multiqc \
    ${FASTP_REPORT_DIR} \
    ${FASTQC_DIR} \
    -o ${MULTIQC_DIR} \
    -n macrogen_clean_multiqc_report.html \
    --force

echo "======================================"
echo "fastp clean + QC finished: $(date)"
echo "======================================"
echo "Sample list:"
echo "${WORK_DIR}/sample_list.txt"
echo ""
echo "Clean reads:"
echo "${CLEAN_DIR}"
echo ""
echo "fastp reports:"
echo "${FASTP_REPORT_DIR}"
echo ""
echo "FastQC output:"
echo "${FASTQC_DIR}"
echo ""
echo "MultiQC report:"
echo "${MULTIQC_DIR}/macrogen_clean_multiqc_report.html"
echo "======================================"
