#!/usr/bin/env python3
"""Export accepted KOfam hits and SGB-level KO trait matrices."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

from trait_io import dict_writer


FOCUS_KOS = ("K01183", "K01233", "K01222", "K03478", "K01207", "K01443")
GENE_FIELDS = ("sgb_id", "mag_id", "orf_id", "kofam_id", "kofam_ids", "kofam_EC",
               "kofam_domain_bitScore", "kofam_full_bitScore")
LONG_FIELDS = ("sgb_id", "ko_id", "members_with_ko", "members_total", "prevalence",
               "orf_count_across_members", "mean_orfs_per_member")


def export(out: Path, dram_dir: str = "03_dram2_all_members") -> None:
    dram = out / dram_dir
    raw = dram / "ANNOTATE/raw-annotations.tsv"
    if not (dram / ".done").is_file() or not raw.is_file():
        raise ValueError("KOfam DRAM2 result is incomplete")

    with (out / "input_manifest.tsv").open(newline="", encoding="utf-8") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    membership = {row["mag_id"]: row["sgb_id"] for row in manifest}
    if not membership or len(membership) != len(manifest):
        raise ValueError("Empty or duplicate MAG in manifest")
    members: dict[str, set[str]] = defaultdict(set)
    for mag, sgb in membership.items():
        members[sgb].add(mag)

    dest = out / "12_final_tables/ko"
    member_counts: Counter[tuple[str, str]] = Counter()
    orf_counts: Counter[tuple[str, str]] = Counter()
    ko_ids = set(FOCUS_KOS)
    observed_mags = set()
    current_mag = None
    current_orfs: set[str] = set()
    current_ko_counts: Counter[str] = Counter()
    total_orfs = 0
    hit_orfs = 0

    def flush_mag() -> None:
        if current_mag is None:
            return
        sgb = membership[current_mag]
        for ko, count in current_ko_counts.items():
            member_counts[sgb, ko] += 1
            orf_counts[sgb, ko] += count

    with dict_writer(dest / "gene_x_ko_hits.tsv.gz", GENE_FIELDS) as writer:
        with raw.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"query_id", "input_fasta", "kofam_id", "kofam_ids", "kofam_EC",
                        "kofam_domain_bitScore", "kofam_full_bitScore"}
            if not required <= set(reader.fieldnames or ()):
                raise ValueError("Unexpected DRAM2 raw-annotations header")
            for row in reader:
                mag = row["input_fasta"]
                orf = row["query_id"]
                if mag not in membership or not orf.startswith(mag + "_contig_"):
                    raise ValueError(f"Foreign DRAM2 ORF: {mag} {orf}")
                if mag != current_mag:
                    flush_mag()
                    if mag in observed_mags:
                        raise ValueError(f"Non-contiguous DRAM2 MAG rows: {mag}")
                    observed_mags.add(mag)
                    current_mag = mag
                    current_orfs = set()
                    current_ko_counts = Counter()
                if orf in current_orfs:
                    raise ValueError(f"Duplicate DRAM2 ORF: {mag} {orf}")
                current_orfs.add(orf)
                total_orfs += 1
                ko = row["kofam_id"].strip()
                if not ko:
                    continue
                if not re.fullmatch(r"K\d{5}", ko):
                    raise ValueError(f"Invalid primary KO: {orf} {ko}")
                candidates = set(re.findall(r"K\d{5}", row["kofam_ids"]))
                if ko not in candidates:
                    raise ValueError(f"Primary KO missing from kofam_ids: {orf} {ko}")
                writer.writerow({"sgb_id": membership[mag], "mag_id": mag, "orf_id": orf,
                                 **{field: row[field] for field in GENE_FIELDS[3:]}})
                current_ko_counts[ko] += 1
                ko_ids.add(ko)
                hit_orfs += 1
        flush_mag()

    if observed_mags != set(membership):
        raise ValueError(f"DRAM2 MAG coverage differs from manifest: "
                         f"{observed_mags ^ set(membership)}")

    ko_order = sorted(ko_ids)
    with (
        dict_writer(dest / "sgb_x_ko.tsv", ["sgb_id", *ko_order]) as prevalence_writer,
        dict_writer(dest / "sgb_x_ko_mean_copy.tsv", ["sgb_id", *ko_order]) as copy_writer,
    ):
        for sgb in sorted(members):
            total = len(members[sgb])
            prevalence = {"sgb_id": sgb}
            mean_copy = {"sgb_id": sgb}
            for ko in ko_order:
                prevalence[ko] = f"{member_counts[sgb, ko] / total:.8g}"
                mean_copy[ko] = f"{orf_counts[sgb, ko] / total:.8g}"
            prevalence_writer.writerow(prevalence)
            copy_writer.writerow(mean_copy)

    with dict_writer(dest / "sgb_x_ko_long_nonzero.tsv", LONG_FIELDS) as writer:
        for sgb, ko in sorted(member_counts):
            total = len(members[sgb])
            positive = member_counts[sgb, ko]
            count = orf_counts[sgb, ko]
            writer.writerow({"sgb_id": sgb, "ko_id": ko,
                             "members_with_ko": positive, "members_total": total,
                             "prevalence": f"{positive / total:.8g}",
                             "orf_count_across_members": count,
                             "mean_orfs_per_member": f"{count / total:.8g}"})

    with dict_writer(dest / "sgb_x_selected_ko_long.tsv", LONG_FIELDS) as writer:
        for sgb in sorted(members):
            total = len(members[sgb])
            for ko in FOCUS_KOS:
                positive = member_counts[sgb, ko]
                count = orf_counts[sgb, ko]
                writer.writerow({"sgb_id": sgb, "ko_id": ko,
                                 "members_with_ko": positive, "members_total": total,
                                 "prevalence": f"{positive / total:.8g}",
                                 "orf_count_across_members": count,
                                 "mean_orfs_per_member": f"{count / total:.8g}"})

    print(f"KOfam traits: MAGs={len(membership)} SGBs={len(members)} "
          f"ORFs={total_orfs} KO_hits={hit_orfs} KOs={len(ko_order)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--dram-dir", default="03_dram2_all_members")
    args = parser.parse_args()
    export(args.out, args.dram_dir)
