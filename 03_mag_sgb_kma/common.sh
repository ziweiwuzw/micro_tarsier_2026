#!/usr/bin/env bash
# Shared configuration for the MetaWRAP/SemiBin2 to KMA workflow.

set -Eeuo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/path/to/project}"
SCRATCH_ROOT="${SCRATCH_ROOT:-/path/to/scratch}"
JOINT_ROOT="${JOINT_ROOT:-${SCRATCH_ROOT}/07_semibin2/by_species}"
SOURCE_ROOT="${SOURCE_ROOT:-${SCRATCH_ROOT}/08_mag_pipeline/all_species_3binners_v2}"
PIPELINE_ROOT="${PIPELINE_ROOT:-${SCRATCH_ROOT}/08_mag_pipeline/all_species_3binners_v2_kma}"
LOG_ROOT="${LOG_ROOT:-${PROJECT_ROOT}/08_vamb_semibin2/logs/mag_all_species_3binners_v2_kma}"
INPUT_READS_MANIFEST="${INPUT_READS_MANIFEST:-/path/to/metadata/04_taxvamb_abundance/manifest.tsv}"
READS_MANIFEST="${READS_MANIFEST:-${PIPELINE_ROOT}/reads.tsv}"
REFERENCE_ROOT="${REFERENCE_ROOT:-${PIPELINE_ROOT}/sgb_reference}"
KMA_ROOT="${KMA_ROOT:-${PIPELINE_ROOT}/kma}"
KMA_ENV="${KMA_ENV:-/path/to/project/software/envs/kma_microbecensus}"

PILOT_GROUPS_TSV="${PILOT_GROUPS_TSV:-${SOURCE_ROOT}/pilot_groups.tsv}"
CANDIDATE_DIR="${CANDIDATE_DIR:-${PIPELINE_ROOT}/candidates/genomes}"
CANDIDATE_MANIFEST="${CANDIDATE_MANIFEST:-${PIPELINE_ROOT}/candidates/candidate_manifest.tsv}"
CHECKM2_OUT="${CHECKM2_OUT:-${PIPELINE_ROOT}/checkm2}"
GUNC_OUT="${GUNC_OUT:-${PIPELINE_ROOT}/gunc}"
HQ_DIR="${HQ_DIR:-${PIPELINE_ROOT}/hq/genomes}"
HQ_MANIFEST="${HQ_MANIFEST:-${PIPELINE_ROOT}/hq/hq_manifest.tsv}"
DREP_OUT="${DREP_OUT:-${PIPELINE_ROOT}/drep}"
SGB_ROOT="${SGB_ROOT:-${PIPELINE_ROOT}/sgb}"
SGB_REP_DIR="${SGB_REP_DIR:-${SGB_ROOT}/representatives}"

GTDB_MMSEQS_ROOT="${GTDB_MMSEQS_ROOT:-${PROJECT_ROOT}/00_databases/mmseqs}"
GTDB_MMSEQS_DB="${GTDB_MMSEQS_DB:-${GTDB_MMSEQS_ROOT}/GTDB}"
TAXVAMB_ENV="${TAXVAMB_ENV:-/path/to/home/software/anaconda3/envs/Taxvamb}"
CONDA_SH="${CONDA_SH:-/path/to/software/Anaconda3/2024.02-1/etc/profile.d/conda.sh}"

GUNC_ENV="${GUNC_ENV:-${SCRATCH_ROOT}/software/anaconda3/envs/gunc}"
GUNC_DB_ROOT="${GUNC_DB_ROOT:-${PROJECT_ROOT}/00_databases/gunc/gtdb214}"
GUNC_DB="${GUNC_DB:-}"
CHECKM1_DB_MODULE="${CHECKM1_DB_MODULE:-CheckM-Database/2015_01_16}"
CHECKM1_DB_ROOT="${CHECKM1_DB_ROOT:-/path/to/software/CheckM-Database/2015_01_16}"
CHECKM2_DB_ROOT="${CHECKM2_DB_ROOT:-${PROJECT_ROOT}/00_databases/checkm2/1.0.1}"
CHECKM2_DB="${CHECKM2_DB:-}"
GTDBTK_MODULE="${GTDBTK_MODULE:-GTDB-Tk/2.4.0}"
GTDBTK_REFERENCE_PATH="${GTDBTK_REFERENCE_PATH:-/path/to/GTDBtk/release220}"
# Do not inherit a mismatched R232 path from an interactive Conda environment.
GTDBTK_DATA_PATH="$GTDBTK_REFERENCE_PATH"
GTDBTK_MASH_DB="${GTDBTK_MASH_DB:-${PIPELINE_ROOT}/gtdbtk_cache/gtdb_r220.msh}"

export PROJECT_ROOT SCRATCH_ROOT JOINT_ROOT SOURCE_ROOT PIPELINE_ROOT LOG_ROOT
export INPUT_READS_MANIFEST READS_MANIFEST REFERENCE_ROOT KMA_ROOT KMA_ENV
export PILOT_GROUPS_TSV CANDIDATE_DIR CANDIDATE_MANIFEST CHECKM2_OUT GUNC_OUT
export HQ_DIR HQ_MANIFEST DREP_OUT SGB_ROOT SGB_REP_DIR
export GTDB_MMSEQS_ROOT GTDB_MMSEQS_DB TAXVAMB_ENV CONDA_SH
export GUNC_ENV GUNC_DB_ROOT GUNC_DB CHECKM1_DB_MODULE CHECKM1_DB_ROOT CHECKM2_DB_ROOT CHECKM2_DB
export GTDBTK_MODULE GTDBTK_REFERENCE_PATH GTDBTK_DATA_PATH GTDBTK_MASH_DB

die() {
    echo "[ERROR] $*" >&2
    exit 1
}

require_file() {
    [[ -s "$1" ]] || die "Missing or empty file: $1"
}

require_dir() {
    [[ -d "$1" ]] || die "Missing directory: $1"
}

data_rows() {
    awk 'NR > 1 {n++} END {print n + 0}' "$1"
}

load_checkm1_database_module() {
    # metaWRAP/1.3 already loads its compatible CheckM 1.0.18/Python-2 stack.
    # Only add the site reference-data module; do not replace its CheckM binary.
    module load "$CHECKM1_DB_MODULE"
    CHECKM_DATA_PATH="$CHECKM1_DB_ROOT"
    require_dir "$CHECKM_DATA_PATH"
    require_file "${CHECKM_DATA_PATH}/.dmanifest"
    command -v checkm >/dev/null 2>&1 || die "CheckM loaded by MetaWRAP is unavailable"
    export CHECKM_DATA_PATH
}

load_checkm2_database() {
    if [[ -z "$CHECKM2_DB" ]]; then
        mapfile -t checkm2_dbs < <(
            find "$CHECKM2_DB_ROOT" -maxdepth 3 -type f -name 'uniref100.KO.1.dmnd' -print 2>/dev/null | sort
        )
        (( ${#checkm2_dbs[@]} == 1 )) \
            || die "Expected one CheckM2 database under ${CHECKM2_DB_ROOT}; found ${#checkm2_dbs[@]}. Set CHECKM2_DB to the installed database."
        CHECKM2_DB="${checkm2_dbs[0]}"
    fi

    require_file "$CHECKM2_DB"
    CHECKM2DB="$CHECKM2_DB"
    export CHECKM2_DB CHECKM2DB
}

load_gunc_environment() {
    if [[ -n "${GUNC_MODULE:-}" ]]; then
        module load "$GUNC_MODULE"
    elif [[ -x "${GUNC_ENV}/bin/gunc" ]]; then
        require_file "$CONDA_SH"
        # shellcheck disable=SC1090
        source "$CONDA_SH"
        conda activate "$GUNC_ENV"
        export PYTHONNOUSERSITE=1
        unset PYTHONPATH PYTHONHOME PIP_USER || true
    else
        die "GUNC is unavailable. Set GUNC_MODULE or install the environment at ${GUNC_ENV}."
    fi

    command -v gunc >/dev/null 2>&1 || die "gunc is not on PATH"

    if [[ -z "$GUNC_DB" ]]; then
        mapfile -t gunc_dbs < <(find "$GUNC_DB_ROOT" -maxdepth 2 -type f -name '*.dmnd' -print 2>/dev/null | sort)
        (( ${#gunc_dbs[@]} == 1 )) || die "Expected one GUNC .dmnd database under ${GUNC_DB_ROOT}; found ${#gunc_dbs[@]}. Set GUNC_DB explicitly."
        GUNC_DB="${gunc_dbs[0]}"
        export GUNC_DB
    fi
    require_file "$GUNC_DB"
}
