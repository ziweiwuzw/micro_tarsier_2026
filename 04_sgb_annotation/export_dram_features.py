#!/usr/bin/env python3
"""Summarize primary DRAM2 annotations across every SGB member MAG."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


SOURCES = {
    "card": "card_id",
    "dram_db": "dram_db_id",
    "kofam": "kofam_id",
    "merops": "merops_id",
    "merops_family": "merops_family",
    "metals": "metals_id",
    "sulfur": "sulfur_id",
}
GENE_FIELDS = ("sgb_id", "mag_id", "orf_id", "source", "feature_id")
SGB_FIELDS = ("sgb_id", "source", "feature_id", "members_with_feature", "members_total",
              "prevalence", "orf_count_across_members", "mean_orfs_per_member")


def write_tsv(path: Path, fields: tuple[str, ...], records: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def export(out: Path) -> None:
    raw = out / "03_dram2_all_members/ANNOTATE/raw-annotations.tsv"
    if not (out / "03_dram2_all_members/.done").is_file() or not raw.is_file():
        raise ValueError("All-member DRAM2 result is incomplete")
    with (out / "input_manifest.tsv").open(newline="", encoding="utf-8") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    membership = {row["mag_id"]: row["sgb_id"] for row in manifest}
    if not membership or len(membership) != len(manifest):
        raise ValueError("Empty or duplicate MAG in manifest")
    members: dict[str, set[str]] = defaultdict(set)
    for mag, sgb in membership.items():
        members[sgb].add(mag)

    genes = []
    by_mag_feature: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    seen_mags = set()
    seen_orfs = set()
    with raw.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"input_fasta", "query_id", *SOURCES.values()} <= set(reader.fieldnames or ()):
            raise ValueError(f"Unexpected DRAM2 raw-annotations header: {raw}")
        for row in reader:
            mag, orf = row["input_fasta"], row["query_id"]
            if mag not in membership or not orf.startswith(mag + "_contig_") or orf in seen_orfs:
                raise ValueError(f"Foreign or duplicate DRAM2 ORF: {mag} {orf}")
            seen_mags.add(mag)
            seen_orfs.add(orf)
            for source, column in SOURCES.items():
                feature = row[column].strip()
                if not feature:
                    continue
                if ";" in feature:
                    raise ValueError(f"Non-unique primary DRAM2 assignment: {orf} {column}={feature}")
                genes.append({"sgb_id": membership[mag], "mag_id": mag, "orf_id": orf,
                              "source": source, "feature_id": feature})
                by_mag_feature[mag, source, feature].add(orf)
    if seen_mags != set(membership):
        raise ValueError(f"DRAM2 MAG coverage differs from manifest: {seen_mags ^ set(membership)}")

    by_sgb_feature: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    for mag, source, feature in by_mag_feature:
        by_sgb_feature[membership[mag], source, feature].add(mag)
    summaries = []
    for (sgb, source, feature), positive in sorted(by_sgb_feature.items()):
        count = sum(len(by_mag_feature[mag, source, feature]) for mag in positive)
        summaries.append({"sgb_id": sgb, "source": source, "feature_id": feature,
                          "members_with_feature": len(positive), "members_total": len(members[sgb]),
                          "prevalence": f"{len(positive) / len(members[sgb]):.8g}",
                          "orf_count_across_members": count,
                          "mean_orfs_per_member": f"{count / len(members[sgb]):.8g}"})
    destination = out / "12_final_tables/dram2"
    write_tsv(destination / "gene_x_dram2_features_long.tsv", GENE_FIELDS, genes)
    write_tsv(destination / "sgb_x_dram2_features_long.tsv", SGB_FIELDS, summaries)
    print(f"DRAM2 primary features: MAGs={len(membership)} ORFs={len(seen_orfs)} "
          f"gene-source rows={len(genes)} SGB-feature rows={len(summaries)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    export(parser.parse_args().out)
