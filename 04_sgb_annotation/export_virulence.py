#!/usr/bin/env python3
"""Export separate VFDB and Toxinome traits plus an ORF-deduplicated summary."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

from trait_io import dict_writer


HIT_FIELDS = ("mag_id", "orf_id", "source", "reference_id", "class", "description",
              "reference_source", "pident", "query_coverage", "subject_coverage",
              "evalue", "bitscore")
TOXINOME_CLASSES = ("toxin", "immunity", "ambiguous", "unclassified")
TOXIN_METADATA_FIELDS = ("toxinome_cluster", "toxinome_protein_toxicity",
                         "toxinome_reference_genome", "toxinome_reference_scaffold")
LONG_FIELDS = ("sgb_id", "feature", "members_with_feature", "members_total", "prevalence",
               "orf_count_across_members", "mean_orfs_per_member")
SUMMARY_FIELDS = (
    "sgb_id", "members_total", "members_with_vfdb_setA", "vfdb_setA_prevalence",
    "vfdb_setA_orfs_across_members", "mean_vfdb_setA_orfs_per_member",
    "members_with_toxinome_toxin", "toxinome_toxin_prevalence",
    "toxinome_toxin_orfs_across_members", "mean_toxinome_toxin_orfs_per_member",
    "members_with_both_sources", "members_with_shared_orf", "shared_orfs_across_members",
    "members_with_union", "union_prevalence", "union_orfs_across_members",
    "mean_union_orfs_per_member",
)


def header_id(description: str, prefix: str) -> str:
    matches = set(re.findall(rf"\b{prefix}\d+\b", description))
    return next(iter(matches)) if len(matches) == 1 else ""


def hit_rows(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not set(HIT_FIELDS) <= set(reader.fieldnames or ()):
            raise ValueError(f"Unexpected homology header: {path}")
        yield from reader


def toxin_metadata(path: Path, wanted: set[str]) -> dict[str, dict[str, str]]:
    found = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"proteinID", "protein_toxicity", "genome", "seqID", "clus"}
        if not required <= set(reader.fieldnames or ()):
            raise ValueError(f"Unexpected Toxinome metadata header: {path}")
        for row in reader:
            protein = row["proteinID"].strip()
            if protein not in wanted:
                continue
            if protein in found:
                raise ValueError(f"Duplicate Toxinome reference ID: {protein}")
            cluster = re.search(r"\bCluster\s+(\d+)\b", row["clus"], flags=re.IGNORECASE)
            found[protein] = {
                "toxinome_cluster": cluster.group(1) if cluster else "",
                "toxinome_protein_toxicity": row["protein_toxicity"].strip(),
                "toxinome_reference_genome": row["genome"].strip(),
                "toxinome_reference_scaffold": row["seqID"].strip(),
            }
    if set(found) != wanted:
        raise ValueError(f"Toxinome hits absent from metadata: {len(wanted - set(found))}")
    return found


def export(out: Path, toxin_csv: Path) -> None:
    with (out / "input_manifest.tsv").open(newline="", encoding="utf-8") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    membership = {row["mag_id"]: row["sgb_id"] for row in manifest}
    if not membership or len(membership) != len(manifest):
        raise ValueError("Empty or duplicate MAG in manifest")
    members: dict[str, set[str]] = defaultdict(set)
    for mag, sgb in membership.items():
        members[sgb].add(mag)

    wanted_toxin_refs = set()
    for mag in membership:
        folder = out / "06_homology" / mag
        if not (folder / ".done").is_file():
            raise ValueError(f"Incomplete homology result: {mag}")
        for row in hit_rows(folder / "toxinome_best_hits.tsv"):
            wanted_toxin_refs.add(row["reference_id"])
    toxin_details = toxin_metadata(toxin_csv, wanted_toxin_refs)

    vf_members: Counter[tuple[str, str]] = Counter()
    vf_orfs: Counter[tuple[str, str]] = Counter()
    toxin_members: Counter[tuple[str, str]] = Counter()
    toxin_orfs: Counter[tuple[str, str]] = Counter()
    summary: dict[str, Counter[str]] = defaultdict(Counter)
    vf_hit_count = 0
    toxin_hit_count = 0
    destination = out / "12_final_tables"

    with (
        dict_writer(destination / "vfdb/gene_x_vfdb_setA.tsv.gz",
                    ["sgb_id", *HIT_FIELDS, "vf_id", "vfc_id"]) as vf_writer,
        dict_writer(destination / "toxinome/gene_x_toxinome.tsv.gz",
                    ["sgb_id", *HIT_FIELDS, *TOXIN_METADATA_FIELDS]) as toxin_writer,
    ):
        for mag, sgb in membership.items():
            folder = out / "06_homology" / mag
            vf_by_category: dict[str, set[str]] = defaultdict(set)
            toxin_by_class: dict[str, set[str]] = defaultdict(set)
            seen_vf = set()
            seen_toxin = set()

            for row in hit_rows(folder / "vfdb_A_best_hits.tsv"):
                orf = row["orf_id"]
                if (row["mag_id"] != mag or row["source"] != "vfdb_A" or
                        row["class"] != "virulence_candidate" or
                        not orf.startswith(mag + "_contig_") or orf in seen_vf):
                    raise ValueError(f"Foreign or duplicate VFDB hit: {mag}: {orf}")
                seen_vf.add(orf)
                vfc_id = header_id(row["description"], "VFC") or "unclassified"
                vf_by_category[vfc_id].add(orf)
                vf_writer.writerow({"sgb_id": sgb, **{field: row[field] for field in HIT_FIELDS},
                                    "vf_id": header_id(row["description"], "VF"),
                                    "vfc_id": vfc_id})
                vf_hit_count += 1

            for row in hit_rows(folder / "toxinome_best_hits.tsv"):
                orf = row["orf_id"]
                category = row["class"]
                if (row["mag_id"] != mag or row["source"] != "toxinome" or
                        category not in TOXINOME_CLASSES or
                        not orf.startswith(mag + "_contig_") or orf in seen_toxin):
                    raise ValueError(f"Foreign or duplicate Toxinome hit: {mag}: {orf}")
                seen_toxin.add(orf)
                toxin_by_class[category].add(orf)
                toxin_writer.writerow({"sgb_id": sgb,
                                       **{field: row[field] for field in HIT_FIELDS},
                                       **toxin_details[row["reference_id"]]})
                toxin_hit_count += 1

            for category, orfs in vf_by_category.items():
                vf_members[sgb, category] += 1
                vf_orfs[sgb, category] += len(orfs)
            for category, orfs in toxin_by_class.items():
                toxin_members[sgb, category] += 1
                toxin_orfs[sgb, category] += len(orfs)

            vf_set = set().union(*vf_by_category.values()) if vf_by_category else set()
            toxin_set = toxin_by_class.get("toxin", set())
            shared = vf_set & toxin_set
            union = vf_set | toxin_set
            stats = summary[sgb]
            stats["vf_members"] += bool(vf_set)
            stats["vf_orfs"] += len(vf_set)
            stats["toxin_members"] += bool(toxin_set)
            stats["toxin_orfs"] += len(toxin_set)
            stats["both_members"] += bool(vf_set and toxin_set)
            stats["shared_members"] += bool(shared)
            stats["shared_orfs"] += len(shared)
            stats["union_members"] += bool(union)
            stats["union_orfs"] += len(union)

    with dict_writer(destination / "vfdb/sgb_x_vfdb_category_long.tsv", LONG_FIELDS) as writer:
        for sgb, category in sorted(vf_members):
            total = len(members[sgb])
            positive = vf_members[sgb, category]
            count = vf_orfs[sgb, category]
            writer.writerow({"sgb_id": sgb, "feature": category,
                             "members_with_feature": positive, "members_total": total,
                             "prevalence": f"{positive / total:.8g}",
                             "orf_count_across_members": count,
                             "mean_orfs_per_member": f"{count / total:.8g}"})

    with dict_writer(destination / "toxinome/sgb_x_toxinome_class_long.tsv",
                     LONG_FIELDS) as writer:
        for sgb in sorted(members):
            total = len(members[sgb])
            for category in TOXINOME_CLASSES:
                positive = toxin_members[sgb, category]
                count = toxin_orfs[sgb, category]
                writer.writerow({"sgb_id": sgb, "feature": category,
                                 "members_with_feature": positive, "members_total": total,
                                 "prevalence": f"{positive / total:.8g}",
                                 "orf_count_across_members": count,
                                 "mean_orfs_per_member": f"{count / total:.8g}"})

    with dict_writer(destination / "virulence/sgb_x_virulence_summary.tsv",
                     SUMMARY_FIELDS) as writer:
        for sgb in sorted(members):
            total = len(members[sgb])
            stats = summary[sgb]
            writer.writerow({
                "sgb_id": sgb, "members_total": total,
                "members_with_vfdb_setA": stats["vf_members"],
                "vfdb_setA_prevalence": f"{stats['vf_members'] / total:.8g}",
                "vfdb_setA_orfs_across_members": stats["vf_orfs"],
                "mean_vfdb_setA_orfs_per_member": f"{stats['vf_orfs'] / total:.8g}",
                "members_with_toxinome_toxin": stats["toxin_members"],
                "toxinome_toxin_prevalence": f"{stats['toxin_members'] / total:.8g}",
                "toxinome_toxin_orfs_across_members": stats["toxin_orfs"],
                "mean_toxinome_toxin_orfs_per_member": f"{stats['toxin_orfs'] / total:.8g}",
                "members_with_both_sources": stats["both_members"],
                "members_with_shared_orf": stats["shared_members"],
                "shared_orfs_across_members": stats["shared_orfs"],
                "members_with_union": stats["union_members"],
                "union_prevalence": f"{stats['union_members'] / total:.8g}",
                "union_orfs_across_members": stats["union_orfs"],
                "mean_union_orfs_per_member": f"{stats['union_orfs'] / total:.8g}",
            })

    print(f"Virulence traits: MAGs={len(membership)} SGBs={len(members)} "
          f"VFDB_hits={vf_hit_count} Toxinome_hits={toxin_hit_count}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--toxin-metadata", type=Path, required=True)
    args = parser.parse_args()
    export(args.out, args.toxin_metadata)
