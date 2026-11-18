#!/usr/bin/env bash
#SBATCH --job-name=nonhost_stats
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --output=nonhost_stats_%j.out
#SBATCH --error=nonhost_stats_%j.err

set -Eeuo pipefail
shopt -s nullglob


# ============================================================
# 功能
#
# 1. 使用 SeqKit 统计去宿主后的双端 FASTQ：
#
#    02_host_removed_reads/reads/<species>/<project>/
#
# 2. 调用同目录中的：
#
#    build_nonhost_summary.py
#
# 3. 生成：
#
#    00_nonhost_file_seqkit_stats.tsv
#    01_nonhost_sample_summary.tsv
#    02_nonhost_species_summary.tsv
#    03_nonhost_species_top5_average.tsv
#    04_nonhost_problem_samples.tsv
#    05_unmatched_species_metadata.tsv
#    06_unmatched_species_ddi.tsv
#    07_skipped_nonhost_fastq_names.tsv
#
# 4. Top 5 使用不同项目之间的平衡抽样。
#
#
# 使用方法
# ------------------------------------------------------------
#
# 正常提交：
#
#   bash shell/summarize_nonhost.slurm
#
# 强制重新统计所有 FASTQ：
#
#   FORCE_SEQKIT=1 bash shell/summarize_nonhost.slurm
#
# 将 Top 5 改为 Top 10：
#
#   TOP_N=10 bash shell/summarize_nonhost.slurm
#
# ============================================================


# ============================================================
# 脚本位置
# ============================================================

SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"


# ============================================================
# 输入与输出
#
# 因为设置了：
#
#   Run from the directory containing 02_host_removed_reads.
#
# 所以下面的 02_host_removed_reads 是相对于该目录的。
# ============================================================

READS_ROOT="02_host_removed_reads/reads"

OUTDIR="02_host_removed_reads/summary_seqkit"

SEQKIT_STATS="${OUTDIR}/00_nonhost_file_seqkit_stats.tsv"

PYTHON_SCRIPT="$(dirname "$SCRIPT_PATH")/build_nonhost_summary.py"


# ============================================================
# 元数据
# ============================================================

PROJECT_META_ROOT="${PROJECT_META_ROOT:?Set PROJECT_META_ROOT to the metadata directory}"

SPECIES_META="${PROJECT_META_ROOT}/ENA_gut_mirco_source_Primates.csv"

RUN_META="${PROJECT_META_ROOT}/all_bioproject_runs_metadata_sort.csv"

MACHADO="${PROJECT_META_ROOT}/Machado_Diet_Breadth_TableS1.csv"


# ============================================================
# 软件
#
# SeqKit 使用 works 环境中的可执行文件。
# Python 使用 module load Anaconda3 提供的 Python。
#
# 不激活 works 环境，防止 NumPy、pandas、pyarrow 混用。
# ============================================================

SEQKIT="${PROJECT_META_ROOT}/software/anaconda3/envs/works/bin/seqkit"


# ============================================================
# 参数
# ============================================================

TOP_N="${TOP_N:-5}"

# FORCE_SEQKIT=0：
#   如果已有完整 SeqKit 结果，并且 FASTQ 没有变化，则复用。
#
# FORCE_SEQKIT=1：
#   强制重新统计全部 FASTQ。
FORCE_SEQKIT="${FORCE_SEQKIT:-0}"


# ============================================================
# 日志函数
# ============================================================

log() {
    printf '[%s] %s\n' \
        "$(date '+%F %T')" \
        "$*"
}


die() {
    log "ERROR: $*" >&2
    exit 1
}


# ============================================================
# 登录节点：自动提交 Slurm
#
# 当运行：
#
#   bash shell/summarize_nonhost.slurm
#
# 脚本只调用 sbatch，不在登录节点执行 SeqKit。
# ============================================================

if [[ -z "${SLURM_JOB_ID:-}" ]]; then

    command -v sbatch >/dev/null 2>&1 || {
        die "sbatch command not found"
    }

    SUBMISSION="$(
        sbatch \
            --parsable \
            --export="ALL,FORCE_SEQKIT=${FORCE_SEQKIT},TOP_N=${TOP_N}" \
            "${SCRIPT_PATH}"
    )"

    JOB_ID="${SUBMISSION%%;*}"

    echo
    echo "Slurm job submitted: ${JOB_ID}"
    echo
    echo "查看任务："
    echo "  squeue -j ${JOB_ID}"
    echo
    echo "查看标准输出："
    echo "  tail -f nonhost_stats_${JOB_ID}.out"
    echo
    echo "查看错误日志："
    echo "  cat nonhost_stats_${JOB_ID}.err"
    echo
    echo "取消任务："
    echo "  scancel ${JOB_ID}"
    echo

    exit 0
fi


# ============================================================
# 以下部分只在 Slurm 计算节点运行
# ============================================================

log "============================================================"
log "Nonhost FASTQ summary started"
log "SLURM_JOB_ID=${SLURM_JOB_ID}"
log "SLURM_JOB_NODELIST=${SLURM_JOB_NODELIST:-NA}"
log "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-NA}"
log "Working directory=$(pwd)"
log "============================================================"


THREADS="${SLURM_CPUS_PER_TASK:-12}"


# ============================================================
# 加载 Anaconda3
# ============================================================

if ! command -v module >/dev/null 2>&1; then

    [[ -f /etc/profile.d/modules.sh ]] && \
        source /etc/profile.d/modules.sh

    [[ -f /usr/share/Modules/init/bash ]] && \
        source /usr/share/Modules/init/bash

fi


command -v module >/dev/null 2>&1 || {
    die "module command is unavailable"
}


log "Loading Anaconda3 module"

module load Anaconda3


# ============================================================
# 清除登录节点中可能继承的 Python 用户包设置
# ============================================================

export PYTHONNOUSERSITE=1

unset PYTHONPATH || true


# EasyBuild 的 Anaconda3 模块通常会提供 EBROOTANACONDA3。
if [[ -n "${EBROOTANACONDA3:-}" ]]; then

    PYTHON="${EBROOTANACONDA3}/bin/python3"

else

    PYTHON="$(command -v python3 || true)"

fi


# ============================================================
# 检查软件
# ============================================================

[[ -x "${SEQKIT}" ]] || {
    die "SeqKit executable not found: ${SEQKIT}"
}


[[ -n "${PYTHON}" && -x "${PYTHON}" ]] || {
    die "Anaconda3 Python executable not found"
}


log "SeqKit executable=${SEQKIT}"
log "Python executable=${PYTHON}"


"${SEQKIT}" version

"${PYTHON}" --version


log "Checking NumPy and pandas"


"${PYTHON}" - <<'PY_CHECK'
import sys
import numpy
import pandas

print("Python executable:", sys.executable)
print("NumPy version:", numpy.__version__)
print("NumPy location:", numpy.__file__)
print("pandas version:", pandas.__version__)
print("pandas location:", pandas.__file__)
PY_CHECK


# ============================================================
# 检查输入
# ============================================================

[[ -d "${READS_ROOT}" ]] || {
    die "Reads directory not found: ${READS_ROOT}"
}


[[ -s "${PYTHON_SCRIPT}" ]] || {
    die "Python script not found or empty: ${PYTHON_SCRIPT}"
}


[[ -s "${SPECIES_META}" ]] || {
    die "Species metadata not found or empty: ${SPECIES_META}"
}


[[ -s "${RUN_META}" ]] || {
    die "Run metadata not found or empty: ${RUN_META}"
}


[[ -s "${MACHADO}" ]] || {
    die "Machado DDI table not found or empty: ${MACHADO}"
}


mkdir -p "${OUTDIR}"


# ============================================================
# 搜索所有 nonhost FASTQ
# ============================================================

log "Searching for nonhost FASTQ files under ${READS_ROOT}"


mapfile -d '' FASTQ_FILES < <(
    find "${READS_ROOT}" \
        -type f \
        \( \
            -name '*.nonhost.fastq.gz' \
            -o -name '*.nonhost.fq.gz' \
            -o -name '*.nonhost.fastq' \
            -o -name '*.nonhost.fq' \
        \) \
        -print0 |
    LC_ALL=C sort -z
)


FASTQ_COUNT="${#FASTQ_FILES[@]}"


if ((FASTQ_COUNT == 0)); then
    die "No nonhost FASTQ files found under ${READS_ROOT}"
fi


log "FASTQ files found: ${FASTQ_COUNT}"


# ============================================================
# 判断是否需要重新运行 SeqKit
# ============================================================

NEED_SEQKIT=0
REBUILD_REASON=""


if [[ "${FORCE_SEQKIT}" == "1" ]]; then

    NEED_SEQKIT=1
    REBUILD_REASON="FORCE_SEQKIT=1"

elif [[ ! -s "${SEQKIT_STATS}" ]]; then

    NEED_SEQKIT=1
    REBUILD_REASON="SeqKit statistics table does not exist"

else

    EXISTING_ROWS="$(
        awk '
            NR > 1 {
                n++
            }

            END {
                print n + 0
            }
        ' "${SEQKIT_STATS}"
    )"


    if [[ "${EXISTING_ROWS}" -ne "${FASTQ_COUNT}" ]]; then

        NEED_SEQKIT=1

        REBUILD_REASON="$(
            printf \
                'FASTQ count changed: current=%s, existing_table=%s' \
                "${FASTQ_COUNT}" \
                "${EXISTING_ROWS}"
        )"

    elif find "${READS_ROOT}" \
        -type f \
        \( \
            -name '*.nonhost.fastq.gz' \
            -o -name '*.nonhost.fq.gz' \
            -o -name '*.nonhost.fastq' \
            -o -name '*.nonhost.fq' \
        \) \
        -newer "${SEQKIT_STATS}" \
        -print \
        -quit |
        grep -q .
    then

        NEED_SEQKIT=1
        REBUILD_REASON="At least one FASTQ is newer than the SeqKit table"

    fi

fi


# ============================================================
# SeqKit 统计
# ============================================================

if [[ "${NEED_SEQKIT}" == "1" ]]; then

    log "SeqKit table will be rebuilt"
    log "Reason: ${REBUILD_REASON}"
    log "Threads: ${THREADS}"


    TMP_STATS="${SEQKIT_STATS}.tmp.${SLURM_JOB_ID}"


    cleanup_tmp() {
        rm -f "${TMP_STATS:-}"
    }


    trap cleanup_tmp EXIT INT TERM


    rm -f "${TMP_STATS}"


    "${SEQKIT}" stats \
        -T \
        -j "${THREADS}" \
        "${FASTQ_FILES[@]}" \
        > "${TMP_STATS}"


    [[ -s "${TMP_STATS}" ]] || {
        die "SeqKit produced an empty output table"
    }


    OUTPUT_ROWS="$(
        awk '
            NR > 1 {
                n++
            }

            END {
                print n + 0
            }
        ' "${TMP_STATS}"
    )"


    if [[ "${OUTPUT_ROWS}" -ne "${FASTQ_COUNT}" ]]; then

        die \
            "SeqKit row count mismatch: expected=${FASTQ_COUNT}, observed=${OUTPUT_ROWS}"

    fi


    mv -f \
        "${TMP_STATS}" \
        "${SEQKIT_STATS}"


    trap - EXIT INT TERM


    log "SeqKit statistics written: ${SEQKIT_STATS}"
    log "SeqKit rows: ${OUTPUT_ROWS}"

else

    log "Reusing existing SeqKit table: ${SEQKIT_STATS}"

fi


# ============================================================
# 运行 Python 汇总
# ============================================================

log "Running Python summary script"
log "Balanced Top N=${TOP_N}"


"${PYTHON}" "${PYTHON_SCRIPT}" \
    --reads-root "${READS_ROOT}" \
    --seqkit-stats "${SEQKIT_STATS}" \
    --species-meta "${SPECIES_META}" \
    --run-meta "${RUN_META}" \
    --machado "${MACHADO}" \
    --outdir "${OUTDIR}" \
    --top-n "${TOP_N}"


# ============================================================
# 检查主要输出
# ============================================================

EXPECTED_OUTPUTS=(
    "${OUTDIR}/01_nonhost_sample_summary.tsv"
    "${OUTDIR}/02_nonhost_species_summary.tsv"
    "${OUTDIR}/03_nonhost_species_top5_average.tsv"
    "${OUTDIR}/04_nonhost_problem_samples.tsv"
    "${OUTDIR}/05_unmatched_species_metadata.tsv"
    "${OUTDIR}/06_unmatched_species_ddi.tsv"
    "${OUTDIR}/07_skipped_nonhost_fastq_names.tsv"
)


for OUTPUT_FILE in "${EXPECTED_OUTPUTS[@]}"; do

    [[ -f "${OUTPUT_FILE}" ]] || {
        die "Expected output was not generated: ${OUTPUT_FILE}"
    }

done


# ============================================================
# 输出简单统计
# ============================================================

SAMPLE_COUNT="$(
    awk '
        NR > 1 {
            n++
        }

        END {
            print n + 0
        }
    ' "${OUTDIR}/01_nonhost_sample_summary.tsv"
)"


PROBLEM_COUNT="$(
    awk '
        NR > 1 {
            n++
        }

        END {
            print n + 0
        }
    ' "${OUTDIR}/04_nonhost_problem_samples.tsv"
)"


SPECIES_COUNT="$(
    awk '
        NR > 1 {
            n++
        }

        END {
            print n + 0
        }
    ' "${OUTDIR}/02_nonhost_species_summary.tsv"
)"


TOP_SPECIES_COUNT="$(
    awk '
        NR > 1 {
            n++
        }

        END {
            print n + 0
        }
    ' "${OUTDIR}/03_nonhost_species_top5_average.tsv"
)"


# ============================================================
# 完成
# ============================================================

log "============================================================"
log "Nonhost summary completed successfully"
log "FASTQ files=${FASTQ_COUNT}"
log "Samples=${SAMPLE_COUNT}"
log "Species=${SPECIES_COUNT}"
log "Top-N species=${TOP_SPECIES_COUNT}"
log "Problem samples=${PROBLEM_COUNT}"
log "Results directory=${OUTDIR}"
log
log "Main outputs:"
log "  ${OUTDIR}/01_nonhost_sample_summary.tsv"
log "  ${OUTDIR}/02_nonhost_species_summary.tsv"
log "  ${OUTDIR}/03_nonhost_species_top5_average.tsv"
log "  ${OUTDIR}/04_nonhost_problem_samples.tsv"
log "============================================================"
