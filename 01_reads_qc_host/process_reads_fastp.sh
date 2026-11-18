#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Fastp + FastQC + MultiQC by species
#
# Project:
#   Set PROJECT_DIR to the directory holding the raw reads.
#
# Output:
#   $PROJECT_DIR/01_clean_data
#
# Environments:
#   rna      : fastp + fastqc
#   multiqc  : multiqc
#
# Test one species:
#   ONLY_SPECIES=Alouatta_caraya bash shell/process_reads_fastp.sh
#
# Run all, keep raw data:
#   bash shell/process_reads_fastp.sh
#
# Run all, delete raw FASTQ after each species succeeds:
#   DELETE_RAW_AFTER_SUCCESS=1 bash shell/process_reads_fastp.sh
#
# Nohup example:
#   DELETE_RAW_AFTER_SUCCESS=1 nohup bash shell/process_reads_fastp.sh \
#     > 01_clean_data/nohup.all.log 2>&1 &
# ============================================================

PPN="${PPN:-12}"

PROJECT_DIR="${PROJECT_DIR:?Set PROJECT_DIR to the raw-read project directory}"
WORK_DIR="${PROJECT_DIR}/01_clean_data"

CLEAN_DIR="${WORK_DIR}/01_clean_reads"
FASTP_REPORT_DIR="${WORK_DIR}/02_fastp_reports"
FASTQC_DIR="${WORK_DIR}/03_fastqc_clean"
MULTIQC_DIR="${WORK_DIR}/04_multiqc_clean"
LOG_DIR="${WORK_DIR}/logs"
MANIFEST_DIR="${WORK_DIR}/manifest"

DELETE_RAW_AFTER_SUCCESS="${DELETE_RAW_AFTER_SUCCESS:-0}"
ONLY_SPECIES="${ONLY_SPECIES:-}"

mkdir -p \
    "${CLEAN_DIR}" \
    "${FASTP_REPORT_DIR}" \
    "${FASTQC_DIR}" \
    "${MULTIQC_DIR}" \
    "${LOG_DIR}" \
    "${MANIFEST_DIR}"

MAIN_LOG="${LOG_DIR}/fastp_fastqc_multiqc_by_species.$(date +%Y%m%d_%H%M%S).log"

echo "============================================================" | tee -a "${MAIN_LOG}"
echo "Start: $(date)" | tee -a "${MAIN_LOG}"
echo "PROJECT_DIR: ${PROJECT_DIR}" | tee -a "${MAIN_LOG}"
echo "WORK_DIR: ${WORK_DIR}" | tee -a "${MAIN_LOG}"
echo "Threads: ${PPN}" | tee -a "${MAIN_LOG}"
echo "DELETE_RAW_AFTER_SUCCESS: ${DELETE_RAW_AFTER_SUCCESS}" | tee -a "${MAIN_LOG}"
echo "ONLY_SPECIES: ${ONLY_SPECIES:-ALL}" | tee -a "${MAIN_LOG}"
echo "Disk usage before running:" | tee -a "${MAIN_LOG}"
df -h "${PROJECT_DIR}" | tee -a "${MAIN_LOG}"
echo "============================================================" | tee -a "${MAIN_LOG}"

cd "${PROJECT_DIR}"

# ============================================================
# Conda setup
# ============================================================

CONDA_SH="${CONDA_SH:?Set CONDA_SH to your Conda initialization script}"

if [[ -s "${CONDA_SH}" ]]; then
    source "${CONDA_SH}"
else
    echo "ERROR: Cannot find conda.sh at ${CONDA_SH}" | tee -a "${MAIN_LOG}"
    exit 1
fi

activate_rna_env() {
    conda activate rna
}

activate_multiqc_env() {
    conda activate multiqc
}

echo "Checking rna environment..." | tee -a "${MAIN_LOG}"
activate_rna_env
echo "CONDA_PREFIX=${CONDA_PREFIX}" | tee -a "${MAIN_LOG}"

which fastp | tee -a "${MAIN_LOG}"
which fastqc | tee -a "${MAIN_LOG}"

fastp --version 2>&1 | tee -a "${MAIN_LOG}" || true
fastqc --version 2>&1 | tee -a "${MAIN_LOG}" || true

echo "Checking multiqc environment..." | tee -a "${MAIN_LOG}"
activate_multiqc_env
echo "CONDA_PREFIX=${CONDA_PREFIX}" | tee -a "${MAIN_LOG}"

which multiqc | tee -a "${MAIN_LOG}"
multiqc --version 2>&1 | tee -a "${MAIN_LOG}" || true

activate_rna_env

# ============================================================
# Raw data roots
# ============================================================

RAW_ROOTS=(
    "special_downloads/DRA009275"
    "special_downloads/Zenodo_5159560"
    "special_downloads/macrogen_raw"
    "kingfisher_download_by_host"
)

ALL_FASTQ="${MANIFEST_DIR}/all_raw_fastq_paths.tsv"
ALL_SAMPLES="${MANIFEST_DIR}/all_samples.tsv"
SPECIES_LIST="${MANIFEST_DIR}/species_list.txt"

: > "${ALL_FASTQ}"

for root in "${RAW_ROOTS[@]}"; do
    if [[ -d "${root}" ]]; then
        find "${root}" -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) \
            ! -path "*/_kingfisher_current_round/*" \
            ! -path "*/_kingfisher_final_check/*" \
            ! -path "*/_kingfisher_progress/*" \
            ! -path "*/macrogen_qc/*" \
            ! -path "*/DRA009275_orginal/*" \
            ! -path "*/01_clean_data/*" \
            | sed 's/\r$//' \
            | sort >> "${ALL_FASTQ}"
    else
        echo "WARNING: raw root not found: ${root}" | tee -a "${MAIN_LOG}"
    fi
done

sort -u "${ALL_FASTQ}" -o "${ALL_FASTQ}"

TOTAL_FASTQ=$(wc -l < "${ALL_FASTQ}")
echo "Total raw FASTQ files detected: ${TOTAL_FASTQ}" | tee -a "${MAIN_LOG}"

if [[ "${TOTAL_FASTQ}" -eq 0 ]]; then
    echo "ERROR: No raw FASTQ files found. Stop." | tee -a "${MAIN_LOG}"
    exit 1
fi

# ============================================================
# Build sample manifest
# Columns:
# source species project sample_id layout r1 r2 se
#
# Important fixes:
#   1. lineterminator="\n" prevents ^M / CRLF in TSV.
#   2. Empty r1/r2/se fields are written as NA, so bash read will not lose columns.
# ============================================================

python3 - "${ALL_FASTQ}" "${ALL_SAMPLES}" <<'PY'
import os
import re
import sys
import csv
from collections import defaultdict

all_fastq, all_samples = sys.argv[1], sys.argv[2]

pair_re = re.compile(r"(.+)_([12])\.(fastq|fq)\.gz$")
plain_re = re.compile(r"(.+)\.(fastq|fq)\.gz$")

def clean_text(x):
    return x.replace("\r", "").strip()

def parse_path(path):
    path = clean_text(path)
    parts = path.split("/")
    source = "UNKNOWN"
    species = "UNKNOWN"
    project = "UNKNOWN"

    if len(parts) >= 4 and parts[0] == "special_downloads":
        source = parts[1]
        if parts[1] in {"DRA009275", "Zenodo_5159560", "macrogen_raw"}:
            species = parts[2]
            project = parts[1]

    elif len(parts) >= 4 and parts[0] == "kingfisher_download_by_host":
        source = "kingfisher_download_by_host"
        species = parts[1]
        project = parts[2]

    return source, species, project

def parse_fastq_name(path):
    path = clean_text(path)
    fn = os.path.basename(path)

    m = pair_re.match(fn)
    if m:
        return m.group(1), m.group(2)

    m = plain_re.match(fn)
    if m:
        return m.group(1), "SE"

    return fn, "UNKNOWN"

groups = defaultdict(lambda: {"r1": "", "r2": "", "se": []})

with open(all_fastq, newline="") as f:
    for line in f:
        path = clean_text(line)
        if not path:
            continue

        source, species, project = parse_path(path)
        sample_id, read_label = parse_fastq_name(path)

        key = (source, species, project, sample_id)

        if read_label == "1":
            groups[key]["r1"] = path
        elif read_label == "2":
            groups[key]["r2"] = path
        else:
            groups[key]["se"].append(path)

rows = []

for (source, species, project, sample_id), d in sorted(groups.items()):
    r1 = d["r1"]
    r2 = d["r2"]
    se_list = d["se"]

    if r1 and r2:
        rows.append({
            "source": source,
            "species": species,
            "project": project,
            "sample_id": sample_id,
            "layout": "PE",
            "r1": r1,
            "r2": r2,
            "se": "NA",
        })

    elif r1 and not r2:
        rows.append({
            "source": source,
            "species": species,
            "project": project,
            "sample_id": sample_id,
            "layout": "ORPHAN_R1",
            "r1": r1,
            "r2": "NA",
            "se": "NA",
        })

    elif r2 and not r1:
        rows.append({
            "source": source,
            "species": species,
            "project": project,
            "sample_id": sample_id,
            "layout": "ORPHAN_R2",
            "r1": "NA",
            "r2": r2,
            "se": "NA",
        })

    else:
        for se_path in se_list:
            sample_id2, _ = parse_fastq_name(se_path)
            rows.append({
                "source": source,
                "species": species,
                "project": project,
                "sample_id": sample_id2,
                "layout": "SE",
                "r1": "NA",
                "r2": "NA",
                "se": se_path,
            })

fieldnames = ["source", "species", "project", "sample_id", "layout", "r1", "r2", "se"]

with open(all_samples, "w", newline="") as out:
    writer = csv.DictWriter(
        out,
        fieldnames=fieldnames,
        delimiter="\t",
        lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
PY

cut -f2 "${ALL_SAMPLES}" | tail -n +2 | sed 's/\r$//' | sort -u > "${SPECIES_LIST}"

if [[ -n "${ONLY_SPECIES}" ]]; then
    grep -Fx "${ONLY_SPECIES}" "${SPECIES_LIST}" > "${SPECIES_LIST}.tmp" || {
        echo "ERROR: ONLY_SPECIES=${ONLY_SPECIES} not found in species list." | tee -a "${MAIN_LOG}"
        echo "Available species examples:" | tee -a "${MAIN_LOG}"
        head -30 "${SPECIES_LIST}" | tee -a "${MAIN_LOG}"
        exit 1
    }
    mv "${SPECIES_LIST}.tmp" "${SPECIES_LIST}"
fi

TOTAL_SAMPLES=$(( $(wc -l < "${ALL_SAMPLES}") - 1 ))
TOTAL_SPECIES=$(wc -l < "${SPECIES_LIST}")

echo "Total samples detected: ${TOTAL_SAMPLES}" | tee -a "${MAIN_LOG}"
echo "Species to process: ${TOTAL_SPECIES}" | tee -a "${MAIN_LOG}"
echo "Species list file: ${SPECIES_LIST}" | tee -a "${MAIN_LOG}"

awk -F'\t' 'NR==1 || $5 ~ /^ORPHAN/' "${ALL_SAMPLES}" > "${MANIFEST_DIR}/orphan_samples.tsv"

ORPHAN_N=$(( $(wc -l < "${MANIFEST_DIR}/orphan_samples.tsv") - 1 ))
if [[ "${ORPHAN_N}" -gt 0 ]]; then
    echo "WARNING: Found orphan R1/R2 samples: ${ORPHAN_N}" | tee -a "${MAIN_LOG}"
    echo "Check: ${MANIFEST_DIR}/orphan_samples.tsv" | tee -a "${MAIN_LOG}"
fi

# ============================================================
# Helper functions
# ============================================================

clean_var_cr() {
    printf "%s" "$1" | tr -d '\r'
}

print_fastp_failure() {
    local sample_label="$1"
    local sample_log="$2"

    echo "ERROR: fastp failed for ${sample_label}" | tee -a "${MAIN_LOG}"
    echo "Fastp log: ${sample_log}" | tee -a "${MAIN_LOG}"

    if [[ -s "${sample_log}" ]]; then
        echo "Last 80 lines of fastp log:" | tee -a "${MAIN_LOG}"
        tail -80 "${sample_log}" | tee -a "${MAIN_LOG}"
    else
        echo "Fastp log is empty or missing." | tee -a "${MAIN_LOG}"
    fi
}

# ============================================================
# Process species one by one
# ============================================================

while IFS= read -r species; do
    species="$(clean_var_cr "${species}")"
    [[ -n "${species}" ]] || continue

    echo "============================================================" | tee -a "${MAIN_LOG}"
    echo "Processing species: ${species}" | tee -a "${MAIN_LOG}"
    echo "Start species: $(date)" | tee -a "${MAIN_LOG}"
    echo "Disk before species:" | tee -a "${MAIN_LOG}"
    df -h "${PROJECT_DIR}" | tee -a "${MAIN_LOG}"
    echo "============================================================" | tee -a "${MAIN_LOG}"

    SPECIES_MANIFEST="${MANIFEST_DIR}/${species}.samples.tsv"
    SPECIES_SUCCESS_RAW_LIST="${MANIFEST_DIR}/${species}.success_raw_fastq_to_delete.list"
    SPECIES_CLEAN_FASTQ_LIST="${MANIFEST_DIR}/${species}.clean_fastq.list"

    awk -F'\t' -v sp="${species}" 'NR==1 || $2==sp' "${ALL_SAMPLES}" > "${SPECIES_MANIFEST}"
    : > "${SPECIES_SUCCESS_RAW_LIST}"

    SPECIES_LOG_DIR="${LOG_DIR}/${species}"
    mkdir -p "${SPECIES_LOG_DIR}"

    activate_rna_env

    # ------------------------------
    # Run fastp sample by sample
    # ------------------------------
    while IFS=$'\t' read -r source sp project sample_id layout r1 r2 se; do
        source="$(clean_var_cr "${source}")"
        sp="$(clean_var_cr "${sp}")"
        project="$(clean_var_cr "${project}")"
        sample_id="$(clean_var_cr "${sample_id}")"
        layout="$(clean_var_cr "${layout}")"
        r1="$(clean_var_cr "${r1}")"
        r2="$(clean_var_cr "${r2}")"
        se="$(clean_var_cr "${se}")"

        if [[ "${layout}" == "layout" ]]; then
            continue
        fi

        if [[ "${layout}" == ORPHAN* ]]; then
            echo "WARNING: Skip orphan sample: ${sp}/${project}/${sample_id}/${layout}" | tee -a "${MAIN_LOG}"
            continue
        fi

        OUT_CLEAN_PROJECT_DIR="${CLEAN_DIR}/${sp}/${project}"
        OUT_FASTP_PROJECT_DIR="${FASTP_REPORT_DIR}/${sp}/${project}"
        OUT_FASTQC_PROJECT_DIR="${FASTQC_DIR}/${sp}/${project}"

        mkdir -p "${OUT_CLEAN_PROJECT_DIR}" "${OUT_FASTP_PROJECT_DIR}" "${OUT_FASTQC_PROJECT_DIR}"

        SAMPLE_LOG="${SPECIES_LOG_DIR}/${project}.${sample_id}.fastp.log"

        if [[ "${layout}" == "PE" ]]; then
            OUT_R1="${OUT_CLEAN_PROJECT_DIR}/${sample_id}_1.clean.fastq.gz"
            OUT_R2="${OUT_CLEAN_PROJECT_DIR}/${sample_id}_2.clean.fastq.gz"
            JSON="${OUT_FASTP_PROJECT_DIR}/${sample_id}.fastp.json"
            HTML="${OUT_FASTP_PROJECT_DIR}/${sample_id}.fastp.html"

            if [[ "${r1}" == "NA" || "${r2}" == "NA" || ! -s "${r1}" || ! -s "${r2}" ]]; then
                echo "ERROR: PE raw input missing for ${sp}/${project}/${sample_id}" | tee -a "${MAIN_LOG}"
                echo "R1=${r1}" | tee -a "${MAIN_LOG}"
                echo "R2=${r2}" | tee -a "${MAIN_LOG}"
                exit 1
            fi

            if [[ -s "${OUT_R1}" && -s "${OUT_R2}" && -s "${JSON}" && -s "${HTML}" ]]; then
                echo "Skip existing PE clean sample: ${sp}/${project}/${sample_id}" | tee -a "${MAIN_LOG}"
            else
                echo "fastp PE: ${sp}/${project}/${sample_id}" | tee -a "${MAIN_LOG}"

                if ! fastp \
                    --in1 "${r1}" \
                    --in2 "${r2}" \
                    --out1 "${OUT_R1}" \
                    --out2 "${OUT_R2}" \
                    --detect_adapter_for_pe \
                    --qualified_quality_phred 25 \
                    --unqualified_percent_limit 40 \
                    --n_base_limit 3 \
                    --length_required 36 \
                    --thread "${PPN}" \
                    --html "${HTML}" \
                    --json "${JSON}" \
                    --report_title "${sample_id} fastp report" \
                    > "${SAMPLE_LOG}" 2>&1
                then
                    print_fastp_failure "${sp}/${project}/${sample_id} PE" "${SAMPLE_LOG}"
                    exit 1
                fi
            fi

            if [[ ! -s "${OUT_R1}" || ! -s "${OUT_R2}" || ! -s "${JSON}" || ! -s "${HTML}" ]]; then
                echo "ERROR: Missing clean PE output/report for ${sp}/${project}/${sample_id}" | tee -a "${MAIN_LOG}"
                exit 1
            fi

            printf "%s\n%s\n" "${r1}" "${r2}" >> "${SPECIES_SUCCESS_RAW_LIST}"

        elif [[ "${layout}" == "SE" ]]; then
            OUT_SE="${OUT_CLEAN_PROJECT_DIR}/${sample_id}.clean.fastq.gz"
            JSON="${OUT_FASTP_PROJECT_DIR}/${sample_id}.fastp.json"
            HTML="${OUT_FASTP_PROJECT_DIR}/${sample_id}.fastp.html"

            if [[ "${se}" == "NA" || ! -s "${se}" ]]; then
                echo "ERROR: SE raw input missing for ${sp}/${project}/${sample_id}" | tee -a "${MAIN_LOG}"
                echo "SE=${se}" | tee -a "${MAIN_LOG}"
                exit 1
            fi

            if [[ -s "${OUT_SE}" && -s "${JSON}" && -s "${HTML}" ]]; then
                echo "Skip existing SE clean sample: ${sp}/${project}/${sample_id}" | tee -a "${MAIN_LOG}"
            else
                echo "fastp SE: ${sp}/${project}/${sample_id}" | tee -a "${MAIN_LOG}"
                echo "Input SE: ${se}" | tee -a "${MAIN_LOG}"

                if ! fastp \
                    --in1 "${se}" \
                    --out1 "${OUT_SE}" \
                    --qualified_quality_phred 25 \
                    --unqualified_percent_limit 40 \
                    --n_base_limit 3 \
                    --length_required 36 \
                    --thread "${PPN}" \
                    --html "${HTML}" \
                    --json "${JSON}" \
                    --report_title "${sample_id} fastp report" \
                    > "${SAMPLE_LOG}" 2>&1
                then
                    print_fastp_failure "${sp}/${project}/${sample_id} SE" "${SAMPLE_LOG}"
                    exit 1
                fi
            fi

            if [[ ! -s "${OUT_SE}" || ! -s "${JSON}" || ! -s "${HTML}" ]]; then
                echo "ERROR: Missing clean SE output/report for ${sp}/${project}/${sample_id}" | tee -a "${MAIN_LOG}"
                exit 1
            fi

            printf "%s\n" "${se}" >> "${SPECIES_SUCCESS_RAW_LIST}"

        else
            echo "ERROR: Unknown layout for ${sp}/${project}/${sample_id}: ${layout}" | tee -a "${MAIN_LOG}"
            exit 1
        fi

    done < "${SPECIES_MANIFEST}"

    sort -u "${SPECIES_SUCCESS_RAW_LIST}" -o "${SPECIES_SUCCESS_RAW_LIST}"

    # ------------------------------
    # FastQC for this species
    # ------------------------------
    activate_rna_env

    echo "Running FastQC for species: ${species}" | tee -a "${MAIN_LOG}"

    find "${CLEAN_DIR}/${species}" -type f -name "*.clean.fastq.gz" | sort > "${SPECIES_CLEAN_FASTQ_LIST}"

    CLEAN_N=$(wc -l < "${SPECIES_CLEAN_FASTQ_LIST}")
    echo "Clean FASTQ files for ${species}: ${CLEAN_N}" | tee -a "${MAIN_LOG}"

    if [[ "${CLEAN_N}" -eq 0 ]]; then
        echo "ERROR: No clean FASTQ found for ${species}" | tee -a "${MAIN_LOG}"
        exit 1
    fi

    mkdir -p "${FASTQC_DIR}/${species}"

    if ! xargs -r -n 20 fastqc \
        --noextract \
        --threads "${PPN}" \
        --outdir "${FASTQC_DIR}/${species}" \
        < "${SPECIES_CLEAN_FASTQ_LIST}"
    then
        echo "ERROR: FastQC failed for species: ${species}" | tee -a "${MAIN_LOG}"
        exit 1
    fi

    # ------------------------------
    # MultiQC for this species
    # ------------------------------
    activate_multiqc_env

    echo "Running MultiQC for species: ${species}" | tee -a "${MAIN_LOG}"

    mkdir -p "${MULTIQC_DIR}/${species}"

    if ! multiqc \
        "${FASTP_REPORT_DIR}/${species}" \
        "${FASTQC_DIR}/${species}" \
        -o "${MULTIQC_DIR}/${species}" \
        -n "${species}_clean_multiqc_report.html" \
        --force \
        > "${LOG_DIR}/${species}.multiqc.log" 2>&1
    then
        echo "ERROR: MultiQC failed for species: ${species}" | tee -a "${MAIN_LOG}"
        echo "MultiQC log: ${LOG_DIR}/${species}.multiqc.log" | tee -a "${MAIN_LOG}"
        tail -80 "${LOG_DIR}/${species}.multiqc.log" | tee -a "${MAIN_LOG}" || true
        exit 1
    fi

    if [[ ! -s "${MULTIQC_DIR}/${species}/${species}_clean_multiqc_report.html" ]]; then
        echo "ERROR: MultiQC report missing for ${species}" | tee -a "${MAIN_LOG}"
        exit 1
    fi

    echo "Species finished successfully: ${species}" | tee -a "${MAIN_LOG}"
    echo "Clean reads: ${CLEAN_DIR}/${species}" | tee -a "${MAIN_LOG}"
    echo "Fastp reports: ${FASTP_REPORT_DIR}/${species}" | tee -a "${MAIN_LOG}"
    echo "FastQC: ${FASTQC_DIR}/${species}" | tee -a "${MAIN_LOG}"
    echo "MultiQC: ${MULTIQC_DIR}/${species}/${species}_clean_multiqc_report.html" | tee -a "${MAIN_LOG}"

    # ------------------------------
    # Optional delete raw FASTQ for this species
    # Only delete files that were successfully processed.
    # Orphan files are not deleted.
    # ------------------------------
    if [[ "${DELETE_RAW_AFTER_SUCCESS}" == "1" ]]; then
        RAW_DELETE_N=$(wc -l < "${SPECIES_SUCCESS_RAW_LIST}")
        echo "DELETE_RAW_AFTER_SUCCESS=1" | tee -a "${MAIN_LOG}"
        echo "Raw FASTQ files to delete for ${species}: ${RAW_DELETE_N}" | tee -a "${MAIN_LOG}"
        echo "Delete list: ${SPECIES_SUCCESS_RAW_LIST}" | tee -a "${MAIN_LOG}"

        while IFS= read -r rawfq; do
            rawfq="$(clean_var_cr "${rawfq}")"
            if [[ -n "${rawfq}" && "${rawfq}" != "NA" && -f "${rawfq}" ]]; then
                rm -f "${rawfq}"
            fi
        done < "${SPECIES_SUCCESS_RAW_LIST}"

        echo "Deleted successfully processed raw FASTQ for species: ${species}" | tee -a "${MAIN_LOG}"
    else
        echo "Raw FASTQ kept for species: ${species}" | tee -a "${MAIN_LOG}"
    fi

    echo "Disk after species:" | tee -a "${MAIN_LOG}"
    df -h "${PROJECT_DIR}" | tee -a "${MAIN_LOG}"
    echo "End species: $(date)" | tee -a "${MAIN_LOG}"

done < "${SPECIES_LIST}"

# ============================================================
# Final global MultiQC
# ============================================================

activate_multiqc_env

echo "============================================================" | tee -a "${MAIN_LOG}"
echo "Running final global MultiQC..." | tee -a "${MAIN_LOG}"
echo "============================================================" | tee -a "${MAIN_LOG}"

mkdir -p "${MULTIQC_DIR}/ALL_SPECIES"

if ! multiqc \
    "${FASTP_REPORT_DIR}" \
    "${FASTQC_DIR}" \
    -o "${MULTIQC_DIR}/ALL_SPECIES" \
    -n "all_species_clean_multiqc_report.html" \
    --force \
    > "${LOG_DIR}/ALL_SPECIES.multiqc.log" 2>&1
then
    echo "ERROR: Final global MultiQC failed." | tee -a "${MAIN_LOG}"
    echo "MultiQC log: ${LOG_DIR}/ALL_SPECIES.multiqc.log" | tee -a "${MAIN_LOG}"
    tail -80 "${LOG_DIR}/ALL_SPECIES.multiqc.log" | tee -a "${MAIN_LOG}" || true
    exit 1
fi

echo "============================================================" | tee -a "${MAIN_LOG}"
echo "Finished all: $(date)" | tee -a "${MAIN_LOG}"
echo "Outputs:" | tee -a "${MAIN_LOG}"
echo "  Work dir: ${WORK_DIR}" | tee -a "${MAIN_LOG}"
echo "  Clean reads: ${CLEAN_DIR}" | tee -a "${MAIN_LOG}"
echo "  Fastp reports: ${FASTP_REPORT_DIR}" | tee -a "${MAIN_LOG}"
echo "  FastQC reports: ${FASTQC_DIR}" | tee -a "${MAIN_LOG}"
echo "  MultiQC reports: ${MULTIQC_DIR}" | tee -a "${MAIN_LOG}"
echo "  Manifest: ${MANIFEST_DIR}" | tee -a "${MAIN_LOG}"
echo "  Main log: ${MAIN_LOG}" | tee -a "${MAIN_LOG}"
echo "Disk after all:" | tee -a "${MAIN_LOG}"
df -h "${PROJECT_DIR}" | tee -a "${MAIN_LOG}"
echo "============================================================" | tee -a "${MAIN_LOG}"
