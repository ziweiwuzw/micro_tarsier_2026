#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Summarize fastp results by species / project / merged species
#
# Input:
#   01_clean_data/manifest/all_samples.tsv
#   01_clean_data/02_fastp_reports/*/*/*.fastp.json
#
# Output:
#   01_clean_data/05_summary_tables/
#
# Merged species rule:
#   All trinomial names are merged to binomial level:
#   Genus_species_subspecies -> Genus_species
#
# Usage:
#   bash shell/summarize_fastp.sh
# ============================================================

PROJECT_DIR="${PROJECT_DIR:?Set PROJECT_DIR to the raw-read project directory}"
WORK_DIR="${PROJECT_DIR}/01_clean_data"

MANIFEST="${WORK_DIR}/manifest/all_samples.tsv"
FASTP_REPORT_DIR="${WORK_DIR}/02_fastp_reports"
OUTDIR="${WORK_DIR}/05_summary_tables"

mkdir -p "${OUTDIR}"

SAMPLE_SUMMARY="${OUTDIR}/sample_fastp_summary.tsv"
SPECIES_ORIGINAL_SUMMARY="${OUTDIR}/species_summary_original_level.tsv"
SPECIES_MERGED_SUMMARY="${OUTDIR}/species_summary_merged_level.tsv"
SPECIES_PROJECT_SUMMARY="${OUTDIR}/species_project_summary_original_level.tsv"
MISSING_JSON="${OUTDIR}/missing_fastp_json_samples.tsv"
LOG="${OUTDIR}/summarize_fastp_results.log"

echo "============================================================" | tee "${LOG}"
echo "Start summary: $(date)" | tee -a "${LOG}"
echo "PROJECT_DIR=${PROJECT_DIR}" | tee -a "${LOG}"
echo "MANIFEST=${MANIFEST}" | tee -a "${LOG}"
echo "FASTP_REPORT_DIR=${FASTP_REPORT_DIR}" | tee -a "${LOG}"
echo "OUTDIR=${OUTDIR}" | tee -a "${LOG}"
echo "============================================================" | tee -a "${LOG}"

if [[ ! -s "${MANIFEST}" ]]; then
    echo "ERROR: manifest not found or empty: ${MANIFEST}" | tee -a "${LOG}"
    exit 1
fi

python3 - "${MANIFEST}" "${FASTP_REPORT_DIR}" "${SAMPLE_SUMMARY}" "${SPECIES_ORIGINAL_SUMMARY}" "${SPECIES_MERGED_SUMMARY}" "${SPECIES_PROJECT_SUMMARY}" "${MISSING_JSON}" <<'PY'
import csv
import json
import os
import sys
from collections import defaultdict

manifest, fastp_report_dir, sample_summary, species_original_summary, species_merged_summary, species_project_summary, missing_json = sys.argv[1:]

def clean(x):
    if x is None:
        return ""
    return str(x).replace("\r", "").strip()

def merged_species(sp):
    """
    Merge all trinomial names to binomial species level.

    Examples:
      Gorilla_beringei_beringei -> Gorilla_beringei
      Gorilla_gorilla_gorilla -> Gorilla_gorilla
      Homo_sapiens_sapiens -> Homo_sapiens
      Pan_troglodytes_verus -> Pan_troglodytes
      Pan_troglodytes_troglodytes -> Pan_troglodytes
      Pan_troglodytes_schweinfurthii -> Pan_troglodytes
    """
    sp = clean(sp)
    parts = sp.split("_")

    if len(parts) >= 3:
        return "_".join(parts[:2])

    return sp

def read_fastp_json(path):
    with open(path) as f:
        data = json.load(f)

    before = data.get("summary", {}).get("before_filtering", {})
    after = data.get("summary", {}).get("after_filtering", {})

    return {
        "raw_reads": int(before.get("total_reads", 0) or 0),
        "raw_bases": int(before.get("total_bases", 0) or 0),
        "clean_reads": int(after.get("total_reads", 0) or 0),
        "clean_bases": int(after.get("total_bases", 0) or 0),
        "raw_q20_bases": int(before.get("q20_bases", 0) or 0),
        "raw_q30_bases": int(before.get("q30_bases", 0) or 0),
        "clean_q20_bases": int(after.get("q20_bases", 0) or 0),
        "clean_q30_bases": int(after.get("q30_bases", 0) or 0),
    }

sample_rows = []
missing_rows = []

with open(manifest, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        source = clean(row.get("source"))
        species_original = clean(row.get("species"))
        species_merged = merged_species(species_original)
        project = clean(row.get("project"))
        sample_id = clean(row.get("sample_id"))
        layout = clean(row.get("layout"))

        if not species_original or not project or not sample_id:
            continue

        if layout.startswith("ORPHAN"):
            continue

        json_path = os.path.join(
            fastp_report_dir,
            species_original,
            project,
            f"{sample_id}.fastp.json"
        )

        found_json = os.path.exists(json_path) and os.path.getsize(json_path) > 0

        if found_json:
            stats = read_fastp_json(json_path)
        else:
            stats = {
                "raw_reads": 0,
                "raw_bases": 0,
                "clean_reads": 0,
                "clean_bases": 0,
                "raw_q20_bases": 0,
                "raw_q30_bases": 0,
                "clean_q20_bases": 0,
                "clean_q30_bases": 0,
            }
            missing_rows.append({
                "source": source,
                "species_original": species_original,
                "species_merged": species_merged,
                "project": project,
                "sample_id": sample_id,
                "layout": layout,
                "expected_json": json_path,
            })

        if layout == "PE":
            raw_read_pairs_or_se_reads = stats["raw_reads"] // 2
            clean_read_pairs_or_se_reads = stats["clean_reads"] // 2
            pe_samples = 1
            se_samples = 0
        elif layout == "SE":
            raw_read_pairs_or_se_reads = stats["raw_reads"]
            clean_read_pairs_or_se_reads = stats["clean_reads"]
            pe_samples = 0
            se_samples = 1
        else:
            raw_read_pairs_or_se_reads = 0
            clean_read_pairs_or_se_reads = 0
            pe_samples = 0
            se_samples = 0

        sample_rows.append({
            "source": source,
            "species_original": species_original,
            "species_merged": species_merged,
            "project": project,
            "sample_id": sample_id,
            "layout": layout,
            "found_fastp_json": "YES" if found_json else "NO",
            "raw_reads_total": stats["raw_reads"],
            "raw_bases_total": stats["raw_bases"],
            "clean_reads_total": stats["clean_reads"],
            "clean_bases_total": stats["clean_bases"],
            "raw_read_pairs_or_se_reads": raw_read_pairs_or_se_reads,
            "clean_read_pairs_or_se_reads": clean_read_pairs_or_se_reads,
            "raw_q20_bases": stats["raw_q20_bases"],
            "raw_q30_bases": stats["raw_q30_bases"],
            "clean_q20_bases": stats["clean_q20_bases"],
            "clean_q30_bases": stats["clean_q30_bases"],
            "fastp_json": json_path,
            "pe_samples": pe_samples,
            "se_samples": se_samples,
        })

sample_fields = [
    "source",
    "species_original",
    "species_merged",
    "project",
    "sample_id",
    "layout",
    "found_fastp_json",
    "raw_reads_total",
    "raw_bases_total",
    "clean_reads_total",
    "clean_bases_total",
    "raw_read_pairs_or_se_reads",
    "clean_read_pairs_or_se_reads",
    "raw_q20_bases",
    "raw_q30_bases",
    "clean_q20_bases",
    "clean_q30_bases",
    "fastp_json",
]

with open(sample_summary, "w", newline="") as out:
    writer = csv.DictWriter(
        out,
        fieldnames=sample_fields,
        delimiter="\t",
        lineterminator="\n"
    )
    writer.writeheader()
    for r in sample_rows:
        writer.writerow({k: r[k] for k in sample_fields})

missing_fields = [
    "source",
    "species_original",
    "species_merged",
    "project",
    "sample_id",
    "layout",
    "expected_json",
]

with open(missing_json, "w", newline="") as out:
    writer = csv.DictWriter(
        out,
        fieldnames=missing_fields,
        delimiter="\t",
        lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(missing_rows)

def write_species_summary(out_path, species_key):
    group = defaultdict(lambda: {
        "project_set": set(),
        "samples": 0,
        "pe_samples": 0,
        "se_samples": 0,
        "missing_fastp_json_samples": 0,
        "raw_reads_total": 0,
        "raw_bases_total": 0,
        "clean_reads_total": 0,
        "clean_bases_total": 0,
        "raw_read_pairs_or_se_reads": 0,
        "clean_read_pairs_or_se_reads": 0,
        "raw_q20_bases": 0,
        "raw_q30_bases": 0,
        "clean_q20_bases": 0,
        "clean_q30_bases": 0,
    })

    for r in sample_rows:
        sp = r[species_key]
        g = group[sp]

        g["project_set"].add(r["project"])
        g["samples"] += 1
        g["pe_samples"] += int(r["pe_samples"])
        g["se_samples"] += int(r["se_samples"])

        if r["found_fastp_json"] != "YES":
            g["missing_fastp_json_samples"] += 1

        for k in [
            "raw_reads_total",
            "raw_bases_total",
            "clean_reads_total",
            "clean_bases_total",
            "raw_read_pairs_or_se_reads",
            "clean_read_pairs_or_se_reads",
            "raw_q20_bases",
            "raw_q30_bases",
            "clean_q20_bases",
            "clean_q30_bases",
        ]:
            g[k] += int(r[k])

    rows = []

    for sp, g in group.items():
        projects = sorted(g["project_set"])

        rows.append({
            "species": sp,
            "dataset_count": len(projects),
            "datasets": ";".join(projects),
            "samples": g["samples"],
            "pe_samples": g["pe_samples"],
            "se_samples": g["se_samples"],
            "missing_fastp_json_samples": g["missing_fastp_json_samples"],
            "raw_reads_total": g["raw_reads_total"],
            "raw_bases_total": g["raw_bases_total"],
            "clean_reads_total": g["clean_reads_total"],
            "clean_bases_total": g["clean_bases_total"],
            "raw_read_pairs_or_se_reads": g["raw_read_pairs_or_se_reads"],
            "clean_read_pairs_or_se_reads": g["clean_read_pairs_or_se_reads"],
            "raw_q20_bases": g["raw_q20_bases"],
            "raw_q30_bases": g["raw_q30_bases"],
            "clean_q20_bases": g["clean_q20_bases"],
            "clean_q30_bases": g["clean_q30_bases"],
        })

    rows.sort(key=lambda x: (-x["samples"], x["species"]))

    fields = [
        "species",
        "dataset_count",
        "datasets",
        "samples",
        "pe_samples",
        "se_samples",
        "missing_fastp_json_samples",
        "raw_reads_total",
        "raw_bases_total",
        "clean_reads_total",
        "clean_bases_total",
        "raw_read_pairs_or_se_reads",
        "clean_read_pairs_or_se_reads",
        "raw_q20_bases",
        "raw_q30_bases",
        "clean_q20_bases",
        "clean_q30_bases",
    ]

    with open(out_path, "w", newline="") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

def write_species_project_summary(out_path):
    group = defaultdict(lambda: {
        "samples": 0,
        "pe_samples": 0,
        "se_samples": 0,
        "missing_fastp_json_samples": 0,
        "raw_reads_total": 0,
        "raw_bases_total": 0,
        "clean_reads_total": 0,
        "clean_bases_total": 0,
        "raw_read_pairs_or_se_reads": 0,
        "clean_read_pairs_or_se_reads": 0,
    })

    for r in sample_rows:
        key = (
            r["species_original"],
            r["species_merged"],
            r["project"],
            r["source"],
        )
        g = group[key]

        g["samples"] += 1
        g["pe_samples"] += int(r["pe_samples"])
        g["se_samples"] += int(r["se_samples"])

        if r["found_fastp_json"] != "YES":
            g["missing_fastp_json_samples"] += 1

        for k in [
            "raw_reads_total",
            "raw_bases_total",
            "clean_reads_total",
            "clean_bases_total",
            "raw_read_pairs_or_se_reads",
            "clean_read_pairs_or_se_reads",
        ]:
            g[k] += int(r[k])

    rows = []

    for (species_original, species_merged, project, source), g in group.items():
        rows.append({
            "species_original": species_original,
            "species_merged": species_merged,
            "project": project,
            "source": source,
            **g
        })

    rows.sort(key=lambda x: (
        x["species_merged"],
        x["species_original"],
        x["project"],
        x["source"],
    ))

    fields = [
        "species_original",
        "species_merged",
        "project",
        "source",
        "samples",
        "pe_samples",
        "se_samples",
        "missing_fastp_json_samples",
        "raw_reads_total",
        "raw_bases_total",
        "clean_reads_total",
        "clean_bases_total",
        "raw_read_pairs_or_se_reads",
        "clean_read_pairs_or_se_reads",
    ]

    with open(out_path, "w", newline="") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

write_species_summary(species_original_summary, "species_original")
write_species_summary(species_merged_summary, "species_merged")
write_species_project_summary(species_project_summary)

print(f"Sample rows: {len(sample_rows)}")
print(f"Missing fastp JSON samples: {len(missing_rows)}")
print(f"Written: {sample_summary}")
print(f"Written: {species_original_summary}")
print(f"Written: {species_merged_summary}")
print(f"Written: {species_project_summary}")
print(f"Written: {missing_json}")
PY

echo "============================================================" | tee -a "${LOG}"
echo "Finished summary: $(date)" | tee -a "${LOG}"
echo "Outputs:" | tee -a "${LOG}"
echo "  ${SAMPLE_SUMMARY}" | tee -a "${LOG}"
echo "  ${SPECIES_ORIGINAL_SUMMARY}" | tee -a "${LOG}"
echo "  ${SPECIES_MERGED_SUMMARY}" | tee -a "${LOG}"
echo "  ${SPECIES_PROJECT_SUMMARY}" | tee -a "${LOG}"
echo "  ${MISSING_JSON}" | tee -a "${LOG}"
echo "============================================================" | tee -a "${LOG}"
