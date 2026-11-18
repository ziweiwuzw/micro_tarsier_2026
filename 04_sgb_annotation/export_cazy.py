#!/usr/bin/env python3
"""Export dbCAN evidence and SGB-level CAZyme trait matrices."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

from trait_io import dict_writer


FAMILY = re.compile(r"(?<![A-Za-z0-9])(?:GH|GT|PL|CE|AA|CBM)\d+")
FOCUS_FAMILIES = ("GH18", "GH19", "GH37", "GH65")
METHODS = ("dbCAN_hmm", "dbCAN_sub", "DIAMOND")
GENE_FIELDS = ("sgb_id", "mag_id", "orf_id", "family", "dbCAN_hmm", "dbCAN_sub",
               "DIAMOND", "methods_supporting_family", "method_count_for_family",
               "selected", "overview_method_count", "recommend_result", "EC", "substrate")
LONG_FIELDS = ("sgb_id", "family", "members_with_family", "members_total", "prevalence",
               "orf_count_across_members", "mean_orfs_per_member")


def parse_families(value: str) -> set[str]:
    return set(FAMILY.findall(value))


def write_matrices(destination: Path, stem: str, families: list[str], members,
                   member_counts, orf_counts) -> None:
    fields = ["sgb_id", *families]
    with (
        dict_writer(destination / f"sgb_x_{stem}_prevalence.tsv", fields) as prev_writer,
        dict_writer(destination / f"sgb_x_{stem}_mean_copy.tsv", fields) as copy_writer,
    ):
        for sgb in sorted(members):
            total = len(members[sgb])
            prevalence = {"sgb_id": sgb}
            mean_copy = {"sgb_id": sgb}
            for family in families:
                prevalence[family] = f"{member_counts[sgb, family] / total:.8g}"
                mean_copy[family] = f"{orf_counts[sgb, family] / total:.8g}"
            prev_writer.writerow(prevalence)
            copy_writer.writerow(mean_copy)


def export(out: Path) -> None:
    with (out / "input_manifest.tsv").open(newline="", encoding="utf-8") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    membership = {row["mag_id"]: row["sgb_id"] for row in manifest}
    if not membership or len(membership) != len(manifest):
        raise ValueError("Empty or duplicate MAG in manifest")
    members: dict[str, set[str]] = defaultdict(set)
    for mag, sgb in membership.items():
        members[sgb].add(mag)

    all_member_counts: Counter[tuple[str, str]] = Counter()
    all_orf_counts: Counter[tuple[str, str]] = Counter()
    selected_member_counts: Counter[tuple[str, str]] = Counter()
    selected_orf_counts: Counter[tuple[str, str]] = Counter()
    families = set(FOCUS_FAMILIES)
    evidence_rows = 0
    destination = out / "12_final_tables/dbcan"

    with dict_writer(destination / "gene_x_cazy_evidence.tsv.gz", GENE_FIELDS) as writer:
        for mag, sgb in membership.items():
            folder = out / "04_dbcan" / mag
            if not (folder / ".done").is_file():
                raise ValueError(f"Incomplete dbCAN result: {mag}")
            path = folder / "overview.tsv"
            all_counts: Counter[str] = Counter()
            selected_counts: Counter[str] = Counter()
            with path.open(newline="", encoding="utf-8") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                required = {"Gene ID", "EC#", *METHODS, "#ofTools", "Recommend Results",
                            "Substrate"}
                if not required <= set(reader.fieldnames or ()):
                    raise ValueError(f"Unexpected dbCAN overview header: {path}")
                seen_orfs = set()
                for row in reader:
                    orf = row["Gene ID"].strip()
                    if not orf.startswith(mag + "_contig_") or orf in seen_orfs:
                        raise ValueError(f"Foreign or duplicate dbCAN ORF: {path}: {orf}")
                    seen_orfs.add(orf)
                    calls = {method: parse_families(row[method]) for method in METHODS}
                    for family in sorted(set().union(*calls.values())):
                        support = [method for method in METHODS if family in calls[method]]
                        selected = len(support) >= 2
                        families.add(family)
                        all_counts[family] += 1
                        if selected:
                            selected_counts[family] += 1
                        writer.writerow({"sgb_id": sgb, "mag_id": mag, "orf_id": orf,
                                         "family": family,
                                         **{method: row[method] for method in METHODS},
                                         "methods_supporting_family": ";".join(support),
                                         "method_count_for_family": len(support),
                                         "selected": int(selected),
                                         "overview_method_count": row["#ofTools"],
                                         "recommend_result": row["Recommend Results"],
                                         "EC": row["EC#"], "substrate": row["Substrate"]})
                        evidence_rows += 1
            for family, count in all_counts.items():
                all_member_counts[sgb, family] += 1
                all_orf_counts[sgb, family] += count
            for family, count in selected_counts.items():
                selected_member_counts[sgb, family] += 1
                selected_orf_counts[sgb, family] += count

    family_order = sorted(families)
    write_matrices(destination, "allcazy", family_order, members,
                   all_member_counts, all_orf_counts)
    write_matrices(destination, "selected_cazy", family_order, members,
                   selected_member_counts, selected_orf_counts)

    with dict_writer(destination / "sgb_x_selected_cazy_long_nonzero.tsv", LONG_FIELDS) as writer:
        for sgb, family in sorted(selected_member_counts):
            total = len(members[sgb])
            positive = selected_member_counts[sgb, family]
            count = selected_orf_counts[sgb, family]
            writer.writerow({"sgb_id": sgb, "family": family,
                             "members_with_family": positive, "members_total": total,
                             "prevalence": f"{positive / total:.8g}",
                             "orf_count_across_members": count,
                             "mean_orfs_per_member": f"{count / total:.8g}"})

    with dict_writer(destination / "sgb_x_selected_cazy_focus_long.tsv", LONG_FIELDS) as writer:
        for sgb in sorted(members):
            total = len(members[sgb])
            for family in FOCUS_FAMILIES:
                positive = selected_member_counts[sgb, family]
                count = selected_orf_counts[sgb, family]
                writer.writerow({"sgb_id": sgb, "family": family,
                                 "members_with_family": positive, "members_total": total,
                                 "prevalence": f"{positive / total:.8g}",
                                 "orf_count_across_members": count,
                                 "mean_orfs_per_member": f"{count / total:.8g}"})

    print(f"dbCAN traits: MAGs={len(membership)} SGBs={len(members)} "
          f"evidence_rows={evidence_rows} families={len(family_order)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    export(parser.parse_args().out)
