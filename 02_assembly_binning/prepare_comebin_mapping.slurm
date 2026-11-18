#!/usr/bin/env bash
#SBATCH --job-name=comebin_map
#
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=7-00:00:00
#SBATCH --output=slurm-%x-%A_%a.out
#SBATCH --error=slurm-%x-%A_%a.err

set -Eeuo pipefail

BASE="${BASE:-/path/to/project}"
SCRATCH_BASE="${SCRATCH_BASE:-/path/to/scratch}"
COMEBIN_ENV="${COMEBIN_ENV:-/path/to/home/software/anaconda3/envs/comebin_env}"
# These two executable locations come from the supplied, already-used
# map_host_group_reads.slurm.
STROBEALIGN_ENV="${STROBEALIGN_ENV:-/path/to/home/software/anaconda3/envs/Taxvamb}"
SAMTOOLS_ENV="${SAMTOOLS_ENV:-/path/to/scratch/software/anaconda3/envs/SemiBin2}"
STROBEALIGN_BIN="${STROBEALIGN_BIN:-${STROBEALIGN_ENV}/bin/strobealign}"
STROBEALIGN_SAMTOOLS_BIN="${STROBEALIGN_SAMTOOLS_BIN:-${SAMTOOLS_ENV}/bin/samtools}"
COMEBIN_SAMTOOLS_BIN="${COMEBIN_SAMTOOLS_BIN:-${COMEBIN_ENV}/bin/samtools}"
BWA_BIN="${BWA_BIN:-${COMEBIN_ENV}/bin/bwa}"
PYTHON_BIN="${PYTHON_BIN:-${COMEBIN_ENV}/bin/python}"
MANIFEST="${MANIFEST:-${BASE}/06_comebin/per_sample_manifest.tsv}"
MIN_CONTIG_LENGTH="${MIN_CONTIG_LENGTH:-2000}"
OUT_ROOT="${OUT_ROOT:-${SCRATCH_BASE}/06_comebin/multisample_by_species_min${MIN_CONTIG_LENGTH}}"

# The final COMEBin run used Strobealign. BWA remains available for older BAMs.
MAPPER="${MAPPER:-strobealign}"
FILTER_UNMAPPED="${FILTER_UNMAPPED:-0}"
THREADS="${SLURM_CPUS_PER_TASK:-32}"
SORT_THREADS="${SORT_THREADS:-8}"
MAP_THREADS="${MAP_THREADS:-$((THREADS - SORT_THREADS))}"
SORT_MEMORY_PER_THREAD="${SORT_MEMORY_PER_THREAD:-4G}"
STROBEALIGN_R="${STROBEALIGN_R:-6}"
CREATE_STROBEALIGN_INDEX="${CREATE_STROBEALIGN_INDEX:-1}"

if [[ -z "${SAMTOOLS_BIN:-}" ]]; then
    if [[ "${MAPPER}" == "strobealign" ]]; then
        SAMTOOLS_BIN="${STROBEALIGN_SAMTOOLS_BIN}"
    else
        SAMTOOLS_BIN="${COMEBIN_SAMTOOLS_BIN}"
    fi
fi

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

[[ "${SLURM_ARRAY_TASK_ID:-}" =~ ^[0-9]+$ ]] \
    || fail "SLURM_ARRAY_TASK_ID must be a positive manifest array_index"
TASK_ID="${SLURM_ARRAY_TASK_ID}"
[[ "${TASK_ID}" -ge 1 ]] || fail "SLURM_ARRAY_TASK_ID must be >=1"
[[ "${MAPPER}" == "bwa" || "${MAPPER}" == "strobealign" ]] \
    || fail "MAPPER must be bwa or strobealign"
[[ "${FILTER_UNMAPPED}" == "0" || "${FILTER_UNMAPPED}" == "1" ]] \
    || fail "FILTER_UNMAPPED must be 0 or 1"
[[ "${CREATE_STROBEALIGN_INDEX}" == "0" || "${CREATE_STROBEALIGN_INDEX}" == "1" ]] \
    || fail "CREATE_STROBEALIGN_INDEX must be 0 or 1"
[[ "${STROBEALIGN_R}" =~ ^[0-9]+$ ]] || fail "STROBEALIGN_R must be a non-negative integer"
for item in THREADS SORT_THREADS MAP_THREADS; do
    value="${!item}"
    [[ "${value}" =~ ^[0-9]+$ && "${value}" -ge 1 ]] \
        || fail "${item} must be a positive integer"
done
[[ -s "${MANIFEST}" ]] || fail "manifest missing or empty: ${MANIFEST}"

EXPECTED_HEADER=$'array_index\tsample_id\tspecies\tproject\tcontigs\tr1\tr2\tassembly_contig_count\tassembly_bytes'
[[ "$(head -n 1 "${MANIFEST}")" == "${EXPECTED_HEADER}" ]] \
    || fail "unexpected manifest header"
ROW="$(awk -F '\t' -v i="${TASK_ID}" 'NR>1 && $1==i {print; n++} END {if(n!=1) exit 1}' "${MANIFEST}")" \
    || fail "array_index ${TASK_ID} does not occur exactly once"
IFS=$'\t' read -r ARRAY_INDEX TARGET_SAMPLE SPECIES TARGET_PROJECT SOURCE_CONTIGS TARGET_R1 TARGET_R2 SOURCE_CONTIG_COUNT SOURCE_CONTIG_BYTES <<< "${ROW}"

for item in TARGET_SAMPLE SPECIES TARGET_PROJECT; do
    value="${!item}"
    [[ -n "${value}" && "${value}" =~ ^[A-Za-z0-9._-]+$ ]] \
        || fail "unsafe ${item}: ${value}"
done
[[ -s "${SOURCE_CONTIGS}" ]] || fail "assembly missing or empty: ${SOURCE_CONTIGS}"
[[ "$(stat -c '%s' "${SOURCE_CONTIGS}")" -eq "${SOURCE_CONTIG_BYTES}" ]] \
    || fail "assembly size differs from manifest: ${SOURCE_CONTIGS}"

SAMPLE_ROOT="${OUT_ROOT}/${SPECIES}/${TARGET_PROJECT}/${TARGET_SAMPLE}"
INPUT_DIR="${SAMPLE_ROOT}/input"
BAM_DIR="${SAMPLE_ROOT}/bam"
RUN_DIR="${SAMPLE_ROOT}/run"
LOCAL_CONTIGS="${INPUT_DIR}/${TARGET_SAMPLE}.min${MIN_CONTIG_LENGTH}.contigs.fa"
FILTER_INFO="${INPUT_DIR}/filter.info"
GROUP_READS="${SAMPLE_ROOT}/same_species_reads.tsv"
MAPPING_DONE="${SAMPLE_ROOT}/mapping.done"
MAPPING_FAILED="${SAMPLE_ROOT}/mapping.failed"
MAPPING_RUNNING="${SAMPLE_ROOT}/mapping.inprogress"
COMEBIN_DONE="${SAMPLE_ROOT}/comebin.done"

mkdir -p "${INPUT_DIR}" "${BAM_DIR}" "${BASE}/06_comebin/logs"

if [[ -s "${COMEBIN_DONE}" && -s "${RUN_DIR}/comebin_res/comebin_res.tsv" && \
      -d "${RUN_DIR}/comebin_res/comebin_res_bins" ]]; then
    echo "COMEBin output already complete; mapping is unnecessary: ${SAMPLE_ROOT}"
    exit 0
fi

on_error() {
    status=$?
    trap - ERR INT TERM
    {
        echo "status=FAILED"
        echo "time=$(date -Is)"
        echo "exit_status=${status}"
        echo "array_index=${TASK_ID}"
        echo "target_sample=${TARGET_SAMPLE}"
        echo "species=${SPECIES}"
        echo "mapper=${MAPPER}"
        echo "command=${BASH_COMMAND:-unknown}"
    } > "${MAPPING_FAILED}.part" 2>/dev/null || true
    mv "${MAPPING_FAILED}.part" "${MAPPING_FAILED}" 2>/dev/null || true
    rm -f -- "${MAPPING_RUNNING}" 2>/dev/null || true
    exit "${status}"
}
trap on_error ERR INT TERM

on_exit() {
    status=$?
    if [[ "${status}" -ne 0 && -e "${MAPPING_RUNNING}" ]]; then
        if [[ ! -s "${MAPPING_FAILED}" ]]; then
            {
                echo "status=FAILED"
                echo "time=$(date -Is)"
                echo "exit_status=${status}"
                echo "array_index=${TASK_ID}"
                echo "target_sample=${TARGET_SAMPLE}"
                echo "species=${SPECIES}"
                echo "mapper=${MAPPER}"
                echo "command=explicit_exit_or_scheduler_signal"
            } > "${MAPPING_FAILED}.part" 2>/dev/null || true
            mv "${MAPPING_FAILED}.part" "${MAPPING_FAILED}" 2>/dev/null || true
        fi
        rm -f -- "${MAPPING_RUNNING}" 2>/dev/null || true
    fi
}
trap on_exit EXIT

{
    echo "status=RUNNING"
    echo "started=$(date -Is)"
    echo "job_id=${SLURM_JOB_ID:-manual}"
    echo "array_index=${TASK_ID}"
    echo "target_sample=${TARGET_SAMPLE}"
    echo "species=${SPECIES}"
    echo "mapper=${MAPPER}"
} > "${MAPPING_RUNNING}.part"
mv "${MAPPING_RUNNING}.part" "${MAPPING_RUNNING}"

echo "CPU mapping target=${TARGET_SAMPLE} species=${SPECIES} mapper=${MAPPER}"
echo "Threads: mapper=${MAP_THREADS}, samtools_sort=${SORT_THREADS}"

module --force purge >/dev/null 2>&1 || true
export PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME PIP_USER || true
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

[[ -x "${SAMTOOLS_BIN}" ]] || fail "samtools executable not found: ${SAMTOOLS_BIN}"
[[ -x "${PYTHON_BIN}" ]] || fail "Python executable not found: ${PYTHON_BIN}"
if [[ "${MAPPER}" == "strobealign" ]]; then
    [[ -x "${STROBEALIGN_BIN}" ]] || fail "Strobealign executable not found: ${STROBEALIGN_BIN}"
    MAPPER_BIN="${STROBEALIGN_BIN}"
else
    [[ -x "${BWA_BIN}" ]] || fail "BWA executable not found: ${BWA_BIN}"
    MAPPER_BIN="${BWA_BIN}"
fi
MAPPER_VERSION="$(${MAPPER_BIN} --version 2>&1 | awk 'NR==1{first=$0} END{print first}' || true)"
SAMTOOLS_VERSION="$(${SAMTOOLS_BIN} --version | awk 'NR==1{first=$0} END{print first}')"

detect_strobealign_profile() {
    "${PYTHON_BIN}" - "$1" <<'PY'
import gzip
import sys
from pathlib import Path

path = Path(sys.argv[1])
opener = gzip.open if path.suffix == ".gz" else open
lengths = []
with opener(path, "rt", encoding="utf-8", errors="strict") as handle:
    while len(lengths) < 500:
        header = handle.readline()
        if not header:
            break
        sequence = handle.readline().strip()
        plus = handle.readline()
        quality = handle.readline().strip()
        if not (header.startswith("@") and plus.startswith("+") and len(sequence) == len(quality)):
            raise SystemExit(f"Malformed FASTQ record in {path}")
        lengths.append(len(sequence))
if not lengths:
    raise SystemExit(f"No FASTQ reads found in {path}")
mean_length = sum(lengths) / len(lengths)
canonical = (50, 75, 100, 125, 150, 250, 400)
print(min(canonical, key=lambda value: (abs(value - mean_length), value)))
PY
}

ensure_strobealign_index() {
    local profile="$1"
    local index_file="${LOCAL_CONTIGS}.r${profile}.sti"
    [[ "${CREATE_STROBEALIGN_INDEX}" == 1 ]] || return 0
    [[ -s "${index_file}" ]] && return 0
    local lock="${index_file}.lock"
    if mkdir "${lock}" 2>/dev/null; then
        if [[ ! -s "${index_file}" ]]; then
            if ! "${STROBEALIGN_BIN}" --create-index -t "${MAP_THREADS}" -r "${profile}" "${LOCAL_CONTIGS}"; then
                rmdir "${lock}" 2>/dev/null || true
                fail "Strobealign index creation failed for profile r${profile}"
            fi
        fi
        rmdir "${lock}"
    else
        local attempt
        for attempt in $(seq 1 720); do
            [[ -s "${index_file}" ]] && break
            [[ -d "${lock}" ]] || break
            sleep 10
        done
        [[ -s "${index_file}" ]] || fail "timed out waiting for Strobealign index: ${index_file}"
    fi
}

GROUP_READS_PART="${GROUP_READS}.part.${SLURM_JOB_ID:-$$}"
awk -F '\t' -v species="${SPECIES}" 'BEGIN{OFS="\t"} NR>1 && $3==species {print $2,$4,$6,$7}' \
    "${MANIFEST}" > "${GROUP_READS_PART}"
mv "${GROUP_READS_PART}" "${GROUP_READS}"
GROUP_SAMPLE_COUNT="$(wc -l < "${GROUP_READS}")"
[[ "${GROUP_SAMPLE_COUNT}" -ge 1 ]] || fail "species group has no reads"

GROUP_READ_BYTES=0
while IFS=$'\t' read -r SOURCE_SAMPLE SOURCE_PROJECT R1 R2; do
    [[ -s "${R1}" ]] || fail "R1 missing or empty: ${R1}"
    [[ -s "${R2}" ]] || fail "R2 missing or empty: ${R2}"
    GROUP_READ_BYTES=$((GROUP_READ_BYTES + $(stat -c '%s' "${R1}") + $(stat -c '%s' "${R2}")))
done < "${GROUP_READS}"

FILTER_START="$(date +%s)"
if [[ ! -s "${LOCAL_CONTIGS}" || ! -s "${FILTER_INFO}" ]]; then
    LOCAL_PART="${LOCAL_CONTIGS}.part.${SLURM_JOB_ID:-$$}"
    INFO_PART="${FILTER_INFO}.part.${SLURM_JOB_ID:-$$}"
    "${PYTHON_BIN}" - "${SOURCE_CONTIGS}" "${LOCAL_PART}" "${MIN_CONTIG_LENGTH}" > "${INFO_PART}" <<'PY'
import sys
from pathlib import Path

source, output, minimum = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3])
seen = set()
kept = total_bp = 0
header = None
parts = []

def publish(handle):
    global kept, total_bp
    if header is None:
        return
    sequence = "".join(parts)
    if len(sequence) < minimum:
        return
    identifier = header[1:].split()[0]
    if not identifier or identifier in seen:
        raise RuntimeError(f"Invalid or duplicate FASTA identifier: {identifier!r}")
    seen.add(identifier)
    handle.write(f">{identifier}\n")
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start:start + 80] + "\n")
    kept += 1
    total_bp += len(sequence)

with source.open(encoding="utf-8") as src, output.open("w", encoding="utf-8", newline="\n") as dst:
    for number, raw in enumerate(src, 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            publish(dst)
            header, parts = line, []
        else:
            if header is None:
                raise RuntimeError(f"Sequence before first header at line {number}")
            parts.append("".join(line.split()).upper())
    publish(dst)
if kept == 0:
    raise RuntimeError(f"No contigs >= {minimum} bp")
print(f"minimum_contig_length={minimum}")
print(f"filtered_contig_count={kept}")
print(f"filtered_total_bp={total_bp}")
PY
    mv "${LOCAL_PART}" "${LOCAL_CONTIGS}"
    mv "${INFO_PART}" "${FILTER_INFO}"
fi
grep -qx "minimum_contig_length=${MIN_CONTIG_LENGTH}" "${FILTER_INFO}" \
    || fail "filtered assembly threshold mismatch"

FILTERED_CONTIG_COUNT="$(awk -F= '$1=="filtered_contig_count"{print $2}' "${FILTER_INFO}")"
FILTERED_TOTAL_BP="$(awk -F= '$1=="filtered_total_bp"{print $2}' "${FILTER_INFO}")"
REFERENCE_SHA256="$(sha256sum "${LOCAL_CONTIGS}" | awk '{print $1}')"
"${SAMTOOLS_BIN}" faidx "${LOCAL_CONTIGS}"
REFERENCE_DICTIONARY_SHA256="$(awk 'BEGIN{OFS="\t"}{print $1,$2}' "${LOCAL_CONTIGS}.fai" | sha256sum | awk '{print $1}')"
FILTER_SECONDS="$(($(date +%s) - FILTER_START))"

if [[ "${MAPPER}" == "bwa" ]]; then
    INDEX_PREFIX="${INPUT_DIR}/${TARGET_SAMPLE}.bwa_index"
    INDEX_COMPLETE=1
    for suffix in amb ann bwt pac sa; do
        [[ -s "${INDEX_PREFIX}.${suffix}" ]] || INDEX_COMPLETE=0
    done
    if [[ "${INDEX_COMPLETE}" != 1 ]]; then
        "${BWA_BIN}" index -p "${INDEX_PREFIX}.part" "${LOCAL_CONTIGS}"
        for suffix in amb ann bwt pac sa; do
            mv "${INDEX_PREFIX}.part.${suffix}" "${INDEX_PREFIX}.${suffix}"
        done
    fi
fi

MAPPING_START="$(date +%s)"
MAPPINGS_CREATED=0
MAPPINGS_REUSED=0
BAM_BYTES_TOTAL=0

while IFS=$'\t' read -r SOURCE_SAMPLE SOURCE_PROJECT R1 R2; do
    BAM="${BAM_DIR}/${SOURCE_SAMPLE}.bam"
    BAM_DONE="${BAM}.done"
    R1_BYTES="$(stat -c '%s' "${R1}")"
    R2_BYTES="$(stat -c '%s' "${R2}")"
    STROBEALIGN_PROFILE=""
    if [[ "${MAPPER}" == "strobealign" ]]; then
        STROBEALIGN_PROFILE="$(detect_strobealign_profile "${R1}")"
        ensure_strobealign_index "${STROBEALIGN_PROFILE}"
    fi
    REUSE=0
    LEGACY_BWA_PROVENANCE=0
    RECORDED_MAPPER="$(awk -F= '$1=="mapper"{print substr($0,index($0,"=")+1); exit}' "${BAM_DONE}" 2>/dev/null || true)"
    if [[ "${MAPPER}" == "bwa" && -z "${RECORDED_MAPPER}" ]]; then
        # BAM markers made by the supplied historical runner predate the
        # mapper= field, but that runner is unambiguously BWA-only.
        LEGACY_BWA_PROVENANCE=1
    fi
    PROFILE_COMPATIBLE=1
    if [[ "${MAPPER}" == "strobealign" ]] && \
       { ! grep -qx "strobealign_profile=${STROBEALIGN_PROFILE}" "${BAM_DONE}" 2>/dev/null || \
         ! grep -qx "strobealign_R=${STROBEALIGN_R}" "${BAM_DONE}" 2>/dev/null; }; then
            PROFILE_COMPATIBLE=0
    fi
    if [[ -s "${BAM}" && -s "${BAM_DONE}" ]] && "${SAMTOOLS_BIN}" quickcheck "${BAM}" >/dev/null 2>&1 && \
       grep -Fqx "reference_sha256=${REFERENCE_SHA256}" "${BAM_DONE}" && \
       grep -Fqx "r1=${R1}" "${BAM_DONE}" && grep -Fqx "r2=${R2}" "${BAM_DONE}" && \
       grep -qx "r1_bytes=${R1_BYTES}" "${BAM_DONE}" && grep -qx "r2_bytes=${R2_BYTES}" "${BAM_DONE}" && \
       grep -qx "filter_unmapped=${FILTER_UNMAPPED}" "${BAM_DONE}" && \
       [[ "${PROFILE_COMPATIBLE}" == 1 ]] && \
       { grep -qx "mapper=${MAPPER}" "${BAM_DONE}" || [[ "${LEGACY_BWA_PROVENANCE}" == 1 ]]; }; then
        BAM_DICTIONARY_SHA256="$("${SAMTOOLS_BIN}" view -H "${BAM}" | awk -F '\t' 'BEGIN{OFS="\t"}$1=="@SQ"{sn=ln="";for(i=2;i<=NF;i++){if($i~/^SN:/){sn=substr($i,4)}else if($i~/^LN:/){ln=substr($i,4)}};print sn,ln}' | sha256sum | awk '{print $1}')"
        [[ "${BAM_DICTIONARY_SHA256}" == "${REFERENCE_DICTIONARY_SHA256}" ]] && REUSE=1
    fi

    if [[ "${REUSE}" == 1 ]]; then
        echo "[reuse ${MAPPER}] ${SOURCE_SAMPLE}.bam"
        if [[ "${LEGACY_BWA_PROVENANCE}" == 1 ]]; then
            {
                cat "${BAM_DONE}"
                echo "mapper=bwa"
                echo "mapper_version=${MAPPER_VERSION}"
                echo "samtools_version=${SAMTOOLS_VERSION}"
                echo "reference_dictionary_sha256=${REFERENCE_DICTIONARY_SHA256}"
                echo "provenance_upgraded=$(date -Is)"
            } > "${BAM_DONE}.part"
            mv "${BAM_DONE}.part" "${BAM_DONE}"
        fi
        MAPPINGS_REUSED=$((MAPPINGS_REUSED + 1))
    else
        echo "[map ${MAPPER}] ${SOURCE_SAMPLE} -> ${TARGET_SAMPLE}"
        BAM_PART="${BAM}.part"
        TMP_PREFIX="${BAM_DIR}/sort_${SOURCE_SAMPLE}_${SLURM_JOB_ID:-$$}"
        rm -f -- "${BAM_PART}" "${BAM_PART}.bai"

        if [[ "${MAPPER}" == "bwa" ]]; then
            if [[ "${FILTER_UNMAPPED}" == 1 ]]; then
                "${BWA_BIN}" mem -v 1 -t "${MAP_THREADS}" "${INDEX_PREFIX}" "${R1}" "${R2}" \
                    | "${SAMTOOLS_BIN}" view -u -F 4 - \
                    | "${SAMTOOLS_BIN}" sort -@ "${SORT_THREADS}" -m "${SORT_MEMORY_PER_THREAD}" -T "${TMP_PREFIX}" -O BAM -o "${BAM_PART}"
            else
                "${BWA_BIN}" mem -v 1 -t "${MAP_THREADS}" "${INDEX_PREFIX}" "${R1}" "${R2}" \
                    | "${SAMTOOLS_BIN}" sort -@ "${SORT_THREADS}" -m "${SORT_MEMORY_PER_THREAD}" -T "${TMP_PREFIX}" -O BAM -o "${BAM_PART}"
            fi
        else
            STROBEALIGN_ARGS=(-t "${MAP_THREADS}" -r "${STROBEALIGN_PROFILE}" -R "${STROBEALIGN_R}")
            [[ "${CREATE_STROBEALIGN_INDEX}" == 1 ]] && STROBEALIGN_ARGS+=(--use-index)
            if [[ "${FILTER_UNMAPPED}" == 1 ]]; then
                "${STROBEALIGN_BIN}" "${STROBEALIGN_ARGS[@]}" "${LOCAL_CONTIGS}" "${R1}" "${R2}" \
                    | "${SAMTOOLS_BIN}" view -u -F 4 - \
                    | "${SAMTOOLS_BIN}" sort -@ "${SORT_THREADS}" -m "${SORT_MEMORY_PER_THREAD}" -T "${TMP_PREFIX}" -O BAM -o "${BAM_PART}"
            else
                "${STROBEALIGN_BIN}" "${STROBEALIGN_ARGS[@]}" "${LOCAL_CONTIGS}" "${R1}" "${R2}" \
                    | "${SAMTOOLS_BIN}" sort -@ "${SORT_THREADS}" -m "${SORT_MEMORY_PER_THREAD}" -T "${TMP_PREFIX}" -O BAM -o "${BAM_PART}"
            fi
        fi

        "${SAMTOOLS_BIN}" quickcheck -v "${BAM_PART}"
        # Do not exit awk early here. With pipefail enabled, early exit closes
        # the pipe and makes a successful `samtools view -H` report SIGPIPE
        # (status 141), which previously caused valid sorted BAMs to fail.
        SORT_ORDER="$("${SAMTOOLS_BIN}" view -H "${BAM_PART}" | awk -F '\t' '
            $1=="@HD" {
                for(i=2;i<=NF;i++) {
                    if($i~/^SO:/) {
                        sort_order=substr($i,4)
                        break
                    }
                }
            }
            END {print sort_order}
        ')"
        [[ "${SORT_ORDER}" == "coordinate" ]] || fail "BAM is not coordinate sorted: ${SOURCE_SAMPLE}"
        BAM_DICTIONARY_SHA256="$("${SAMTOOLS_BIN}" view -H "${BAM_PART}" | awk -F '\t' 'BEGIN{OFS="\t"}$1=="@SQ"{sn=ln="";for(i=2;i<=NF;i++){if($i~/^SN:/){sn=substr($i,4)}else if($i~/^LN:/){ln=substr($i,4)}};print sn,ln}' | sha256sum | awk '{print $1}')"
        [[ "${BAM_DICTIONARY_SHA256}" == "${REFERENCE_DICTIONARY_SHA256}" ]] \
            || fail "BAM/reference dictionary mismatch: ${SOURCE_SAMPLE}"
        mv "${BAM_PART}" "${BAM}"
        "${SAMTOOLS_BIN}" index -@ "${SORT_THREADS}" "${BAM}"
        {
            echo "status=PASS"
            echo "completed=$(date -Is)"
            echo "target_sample=${TARGET_SAMPLE}"
            echo "source_sample=${SOURCE_SAMPLE}"
            echo "reference_sha256=${REFERENCE_SHA256}"
            echo "reference_dictionary_sha256=${REFERENCE_DICTIONARY_SHA256}"
            echo "r1=${R1}"
            echo "r2=${R2}"
            echo "r1_bytes=${R1_BYTES}"
            echo "r2_bytes=${R2_BYTES}"
            echo "filter_unmapped=${FILTER_UNMAPPED}"
            echo "mapper=${MAPPER}"
            echo "mapper_version=${MAPPER_VERSION}"
            [[ -n "${STROBEALIGN_PROFILE}" ]] && echo "strobealign_profile=${STROBEALIGN_PROFILE}"
            [[ "${MAPPER}" == "strobealign" ]] && echo "strobealign_R=${STROBEALIGN_R}"
            echo "samtools_version=${SAMTOOLS_VERSION}"
        } > "${BAM_DONE}.part"
        mv "${BAM_DONE}.part" "${BAM_DONE}"
        MAPPINGS_CREATED=$((MAPPINGS_CREATED + 1))
    fi
    [[ -s "${BAM}.bai" ]] || "${SAMTOOLS_BIN}" index -@ "${SORT_THREADS}" "${BAM}"
    BAM_BYTES_TOTAL=$((BAM_BYTES_TOTAL + $(stat -c '%s' "${BAM}")))
done < "${GROUP_READS}"

BAM_COUNT="$(find "${BAM_DIR}" -maxdepth 1 -type f -name '*.bam' | wc -l)"
[[ "${BAM_COUNT}" -eq "${GROUP_SAMPLE_COUNT}" ]] \
    || fail "expected ${GROUP_SAMPLE_COUNT} BAMs, found ${BAM_COUNT}"
MAPPING_SECONDS="$(($(date +%s) - MAPPING_START))"

{
    echo "status=PASS"
    echo "completed=$(date -Is)"
    echo "array_index=${TASK_ID}"
    echo "target_sample=${TARGET_SAMPLE}"
    echo "species=${SPECIES}"
    echo "target_project=${TARGET_PROJECT}"
    echo "minimum_contig_length=${MIN_CONTIG_LENGTH}"
    echo "filtered_contig_count=${FILTERED_CONTIG_COUNT}"
    echo "filtered_total_bp=${FILTERED_TOTAL_BP}"
    echo "reference_sha256=${REFERENCE_SHA256}"
    echo "reference_dictionary_sha256=${REFERENCE_DICTIONARY_SHA256}"
    echo "coverage_sample_count=${GROUP_SAMPLE_COUNT}"
    echo "group_read_bytes=${GROUP_READ_BYTES}"
    echo "bam_count=${BAM_COUNT}"
    echo "bam_bytes=${BAM_BYTES_TOTAL}"
    echo "mappings_created=${MAPPINGS_CREATED}"
    echo "mappings_reused=${MAPPINGS_REUSED}"
    echo "mapper=${MAPPER}"
    echo "mapper_version=${MAPPER_VERSION}"
    echo "filter_unmapped=${FILTER_UNMAPPED}"
    echo "filter_seconds=${FILTER_SECONDS}"
    echo "mapping_seconds=${MAPPING_SECONDS}"
} > "${MAPPING_DONE}.part"
mv "${MAPPING_DONE}.part" "${MAPPING_DONE}"
rm -f -- "${MAPPING_FAILED}" "${MAPPING_RUNNING}"

echo "Mapping complete: target=${TARGET_SAMPLE} BAMs=${BAM_COUNT} mapper=${MAPPER}"
echo "Mapping marker: ${MAPPING_DONE}"
