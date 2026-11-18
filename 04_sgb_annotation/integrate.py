#!/usr/bin/env python3
"""Integrate separate MAG evidence into SGB prevalence and weighted sample tables."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

from annotation import check_abundance_columns, rows, write_table


HOMOLOGY_TRAITS = {
    ("vfdb_A", "virulence_candidate"): "VFDB::SetA_candidate",
    ("vfdb_B", "virulence_candidate"): "VFDB::SetB_candidate",
    ("toxinome", "toxin"): "Toxinome::toxin_candidate",
    ("toxinome", "immunity"): "Toxinome::immunity_candidate",
    ("toxinome", "ambiguous"): "Toxinome::ambiguous",
}


def read_optional(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pipeline-root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    manifest = rows(args.out / "input_manifest.tsv")
    mag_to_sgb = {r["mag_id"]: r["sgb_id"] for r in manifest}
    sgb_order = list(dict.fromkeys(r["sgb_id"] for r in manifest))
    by_sgb: dict[str, list[str]] = defaultdict(list)
    for mag, sgb in mag_to_sgb.items():
        by_sgb[sgb].append(mag)
    trait_mag: dict[str, set[str]] = defaultdict(set)
    evidence = []
    for mag, sgb in mag_to_sgb.items():
        homo = args.out / "06_homology" / mag
        if not (homo / ".done").is_file():
            raise ValueError(f"Incomplete homology result: {mag}")
        best_hits = read_optional(homo / "best_hits.tsv")
        for hit in best_hits:
            if hit["mag_id"] != mag or not hit["orf_id"].startswith(mag + "_contig_"):
                raise ValueError(f"Foreign ORF in {homo}: {hit}")
            trait = HOMOLOGY_TRAITS.get((hit["source"], hit["class"]))
            if trait is None:
                if hit["source"] == "toxinome" and hit["class"] == "unclassified":
                    continue
                raise ValueError(f"Unexpected homology class: {hit['source']} {hit['class']}")
            trait_mag[trait].add(mag)
            evidence.append({"mag_id": mag, "sgb_id": sgb, "trait_id": trait,
                             "orf_id": hit["orf_id"], "source": hit["source"],
                             "reference_id": hit["reference_id"], "bitscore": hit["bitscore"]})
    cazy = rows(args.out / "10_trait_matrix/sgb_x_cazy_family_long.tsv")
    prevalence: dict[tuple[str, str], float] = {}
    trait_ids = set(HOMOLOGY_TRAITS.values())
    for row in cazy:
        sgb, family = row["sgb_id"], row["family"]
        trait = f"CAZy::{family}"
        prevalence[sgb, trait] = float(row["prevalence"])
        trait_ids.add(trait)
    trait_ids = sorted(trait_ids)
    matrix = []
    for sgb in sgb_order:
        members = by_sgb[sgb]
        row = {"sgb_id": sgb}
        for trait in trait_ids:
            value = prevalence.get((sgb, trait))
            if value is None:
                value = len(trait_mag[trait].intersection(members)) / len(members)
            row[trait] = f"{value:.8g}"
            prevalence[sgb, trait] = value
        matrix.append(row)
    out = args.out / "10_trait_matrix"
    write_table(out / "homology_best_hit_evidence.tsv", ["mag_id", "sgb_id", "trait_id", "orf_id", "source",
                                                     "reference_id", "bitscore"], evidence)
    write_table(out / "sgb_x_trait_prevalence.tsv", ["sgb_id", *trait_ids], matrix)
    abundance = args.pipeline_root / "tables/sample_x_sgb_tpm.tsv"
    with abundance.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        check_abundance_columns(reader.fieldnames, sgb_order, args.out)
        weighted = []
        seen = set()
        for sample in reader:
            sample_id = sample["sample_id"]
            if not sample_id or sample_id in seen:
                raise ValueError(f"Invalid sample ID: {sample_id}")
            seen.add(sample_id)
            tpm = {sgb: float(sample[sgb]) for sgb in sgb_order}
            if any(not math.isfinite(v) or v < 0 for v in tpm.values()):
                raise ValueError(f"Invalid Salmon TPM for {sample_id}")
            weighted.append({"sample_id": sample_id, **{
                trait: f"{sum(tpm[sgb] * prevalence[sgb, trait] for sgb in sgb_order):.8g}"
                for trait in trait_ids}})
    if not weighted:
        raise ValueError("No sample rows in Salmon TPM matrix")
    write_table(args.out / "11_abundance_weighted/sample_x_trait_weighted_tpm.tsv",
                ["sample_id", *trait_ids], weighted)
    print(f"SGBs={len(sgb_order)} traits={len(trait_ids)} samples={len(weighted)}")


if __name__ == "__main__":
    main()
