#!/usr/bin/env bash
# CPU COMEBin using same-species Strobealign BAMs.
# Reuses existing BAMs and maps missing BAMs with prepare_comebin_mapping.slurm.
# Submit as an array over INDEX_FILE positions.
#SBATCH --job-name=comebin_cpu
#
#SBATCH --partition=sapphire
#SBATCH --qos=publiccpu
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=4-00:00:00
#SBATCH --output=slurm-%x-%A_%a.out
#SBATCH --error=slurm-%x-%A_%a.err

set -Eeuo pipefail

BASE="${BASE:-/path/to/project}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CPU_RUNNER="${CPU_RUNNER:-${HERE}/prepare_comebin_mapping.slurm}"
MANIFEST="${MANIFEST:-${BASE}/06_comebin/per_sample_manifest.tsv}"
INDEX_FILE="${INDEX_FILE:?Set INDEX_FILE to a one-column array-index list}"
OUT_ROOT="${OUT_ROOT:-/path/to/scratch/06_comebin/multisample_by_species_min2000}"
COMEBIN_ENV="${COMEBIN_ENV:-/path/to/home/software/anaconda3/envs/comebin_env}"
CONDA_BASE="${CONDA_BASE:-/path/to/software/Anaconda3/2024.02-1}"
THREADS="${SLURM_CPUS_PER_TASK:-64}"
LEIDEN_WORKERS="${LEIDEN_WORKERS:-8}"
N_VIEWS="${N_VIEWS:-6}"
CLEAN_MAPPING_AFTER_SUCCESS="${CLEAN_MAPPING_AFTER_SUCCESS:-1}"

fail() { echo "ERROR: $*" >&2; exit 1; }
[[ -s "${INDEX_FILE}" ]] || fail "missing INDEX_FILE: ${INDEX_FILE}"
[[ -s "${MANIFEST}" ]] || fail "missing MANIFEST: ${MANIFEST}"
[[ -s "${CPU_RUNNER}" ]] || fail "missing CPU_RUNNER: ${CPU_RUNNER}"
[[ "${SLURM_ARRAY_TASK_ID:-}" =~ ^[0-9]+$ ]] || fail "submit as a Slurm array"

ARRAY_INDEX="$(awk -v n="${SLURM_ARRAY_TASK_ID}" 'NF && $1 !~ /^#/ {i++; if (i == n) {print $1; exit}}' "${INDEX_FILE}")"
[[ "${ARRAY_INDEX}" =~ ^[0-9]+$ ]] || fail "no index for array task ${SLURM_ARRAY_TASK_ID}"
ROW="$(awk -F '\t' -v i="${ARRAY_INDEX}" 'NR > 1 && $1 == i {print; n++} END {if (n != 1) exit 1}' "${MANIFEST}")" \
    || fail "array_index ${ARRAY_INDEX} is not unique in manifest"
IFS=$'\t' read -r _ TARGET_SAMPLE SPECIES TARGET_PROJECT SOURCE_CONTIGS _ <<< "${ROW}"
SAMPLE_ROOT="${OUT_ROOT}/${SPECIES}/${TARGET_PROJECT}/${TARGET_SAMPLE}"
INPUT_DIR="${SAMPLE_ROOT}/input"
BAM_DIR="${SAMPLE_ROOT}/bam"
RUN_DIR="${SAMPLE_ROOT}/run"
LOCAL_CONTIGS="${INPUT_DIR}/${TARGET_SAMPLE}.min2000.contigs.fa"
DONE_FILE="${SAMPLE_ROOT}/comebin.done"
FAILED_FILE="${SAMPLE_ROOT}/comebin.failed"
mkdir -p "${INPUT_DIR}" "${BAM_DIR}" "${BASE}/06_comebin/logs"

if [[ -s "${DONE_FILE}" && -s "${RUN_DIR}/comebin_res/comebin_res.tsv" && -d "${RUN_DIR}/comebin_res/comebin_res_bins" ]]; then
    echo "Already complete: ${SAMPLE_ROOT}"
    exit 0
fi

on_error() {
    status=$?
    {
        echo "status=FAILED"
        echo "time=$(date -Is)"
        echo "array_index=${ARRAY_INDEX}"
        echo "target_sample=${TARGET_SAMPLE}"
        echo "species=${SPECIES}"
        echo "project=${TARGET_PROJECT}"
        echo "device=cpu"
        echo "exit_status=${status}"
        echo "command=${BASH_COMMAND:-unknown}"
    } > "${FAILED_FILE}.part" 2>/dev/null || true
    mv "${FAILED_FILE}.part" "${FAILED_FILE}" 2>/dev/null || true
    exit "${status}"
}
trap on_error ERR

echo "CPU COMEBin target=${TARGET_SAMPLE} species=${SPECIES} index=${ARRAY_INDEX}"
echo "Resources: threads=${THREADS} leiden_workers=${LEIDEN_WORKERS}"

# Reuse existing Strobealign BAMs; map any missing same-species BAMs.
SLURM_ARRAY_TASK_ID="${ARRAY_INDEX}" MAPPER=strobealign \
    bash "${CPU_RUNNER}"

[[ -s "${LOCAL_CONTIGS}" ]] || fail "missing filtered contigs: ${LOCAL_CONTIGS}"
[[ -s "${SAMPLE_ROOT}/mapping.done" ]] || fail "missing mapping.done: ${SAMPLE_ROOT}"
grep -qx 'status=PASS' "${SAMPLE_ROOT}/mapping.done" || fail "mapping.done is not PASS"
FIRST_BAM="$(find "${BAM_DIR}" -maxdepth 1 -type f -name '*.bam' -size +0c -print -quit)"
[[ -n "${FIRST_BAM}" ]] || fail "no BAMs found in ${BAM_DIR}"

source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${COMEBIN_ENV}"
export PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME PIP_USER || true
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
command -v run_comebin.sh >/dev/null 2>&1 || fail "run_comebin.sh not found in ${COMEBIN_ENV}"

ARGS=(-a "${LOCAL_CONTIGS}" -p "${BAM_DIR}" -o "${RUN_DIR}" -n "${N_VIEWS}" -d cpu -t "${THREADS}" -s 42 -w "${LEIDEN_WORKERS}")
echo "Command: run_comebin.sh ${ARGS[*]}"
run_comebin.sh "${ARGS[@]}"

RESULT_TSV="${RUN_DIR}/comebin_res/comebin_res.tsv"
RESULT_BINS="${RUN_DIR}/comebin_res/comebin_res_bins"
[[ -s "${RESULT_TSV}" ]] || fail "missing result TSV: ${RESULT_TSV}"
[[ -d "${RESULT_BINS}" ]] || fail "missing bins directory: ${RESULT_BINS}"
BIN_COUNT="$(find "${RESULT_BINS}" -maxdepth 1 -type f -name '*.fa' | wc -l)"
(( BIN_COUNT > 0 )) || fail "no FASTA bins produced"

{
    echo "status=PASS"
    echo "completed=$(date -Is)"
    echo "array_index=${ARRAY_INDEX}"
    echo "target_sample=${TARGET_SAMPLE}"
    echo "species=${SPECIES}"
    echo "target_project=${TARGET_PROJECT}"
    echo "mapper=strobealign"
    echo "device=cpu"
    echo "bin_count=${BIN_COUNT}"
    echo "result_tsv=${RESULT_TSV}"
    echo "result_bins=${RESULT_BINS}"
} > "${DONE_FILE}.part"
mv "${DONE_FILE}.part" "${DONE_FILE}"
rm -f "${FAILED_FILE}"

if [[ "${CLEAN_MAPPING_AFTER_SUCCESS}" == 1 ]]; then
    rm -rf -- "${INPUT_DIR}" "${BAM_DIR}"
fi
echo "CPU COMEBin complete: target=${TARGET_SAMPLE} bins=${BIN_COUNT}"
