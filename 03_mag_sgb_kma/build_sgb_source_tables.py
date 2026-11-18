#!/usr/bin/env python3
"""Build SGB and MAG assembly-source tables from the completed MAG pipeline.

The sample IDs here describe where contigs were assembled, not the samples in
which KMA later detects an SGB.
"""

from __future__ import annotations

import argparse
import csv
import os
from collections import Counter, defaultdict
from pathlib import Path


def rows(path: Path, required: set[str], delimiter: str = "\t") -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        result = list(reader)
    if not result:
        raise ValueError(f"{path}: no data rows")
    return result


def source_counts(fasta: Path, sample_project: dict[str, str]) -> Counter[str]:
    if not fasta.is_file():
        raise FileNotFoundError(fasta)
    counts: Counter[str] = Counter()
    with fasta.open("rb", buffering=1024 * 1024) as handle:
        for line in handle:
            if not line.startswith(b">"):
                continue
            token = line[1:].split(maxsplit=1)[0].decode("ascii")
            prefix, separator, _ = token.partition(":")
            if not separator:
                raise ValueError(f"{fasta}: contig lacks sample prefix: {token}")
            if prefix in sample_project:
                sample = prefix
            elif prefix.endswith("_final") and prefix[:-6] in sample_project:
                sample = prefix[:-6]
            else:
                raise ValueError(
                    f"{fasta}: sample prefix {prefix!r} absent from reads.tsv"
                )
            counts[sample] += 1
    if not counts:
        raise ValueError(f"{fasta}: no contig headers")
    return counts


def write_tsv(path: Path, fieldnames: list[str], data: list[dict[str, str]]) -> None:
    temporary = path.with_name(path.name + f".part.{os.getpid()}")
    try:
        with temporary.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(data)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    root = args.root

    sample_rows = rows(root / "reads.tsv", {"sample_id", "project", "group"})
    sample_project: dict[str, str] = {}
    sample_group: dict[str, str] = {}
    for row in sample_rows:
        sample = row["sample_id"].strip()
        if not sample or sample in sample_project:
            raise ValueError(f"Duplicate or empty sample_id in reads.tsv: {sample!r}")
        sample_project[sample] = row["project"].strip() or "NA"
        sample_group[sample] = row["group"].strip() or "NA"

    membership = rows(
        root / "sgb/sgb_membership.tsv",
        {"sgb_id", "mag_id", "is_representative"},
    )
    representatives = rows(
        root / "sgb/sgb_representatives.tsv",
        {"sgb_id", "representative_mag", "completeness", "contamination"},
    )
    hq_rows = rows(
        root / "hq/hq_manifest.tsv",
        {"mag_id", "selected_hq", "group", "completeness", "contamination"},
    )
    hq = {row["mag_id"]: row for row in hq_rows}
    if len(hq) != len(hq_rows):
        raise ValueError("Duplicate mag_id in hq_manifest.tsv")
    reps = {row["sgb_id"]: row for row in representatives}
    if len(reps) != len(representatives):
        raise ValueError("Duplicate sgb_id in sgb_representatives.tsv")

    members_by_sgb: dict[str, list[str]] = defaultdict(list)
    flagged_reps: dict[str, list[str]] = defaultdict(list)
    mag_seen: set[str] = set()
    mag_rows: list[dict[str, str]] = []
    sources_by_mag: dict[str, set[str]] = {}
    alias_by_mag: dict[str, str] = {}
    aliases_seen: set[str] = set()
    for row in membership:
        sgb, mag = row["sgb_id"], row["mag_id"]
        if mag in mag_seen:
            raise ValueError(f"MAG occurs more than once in membership: {mag}")
        mag_seen.add(mag)
        if sgb not in reps:
            raise ValueError(f"Membership SGB absent from representatives: {sgb}")
        meta = hq.get(mag)
        if meta is None or meta["selected_hq"].strip().casefold() != "true":
            raise ValueError(f"Member MAG missing from selected HQ set: {mag}")
        if row["is_representative"].strip().casefold() == "true":
            flagged_reps[sgb].append(mag)
        elif row["is_representative"].strip().casefold() != "false":
            raise ValueError(f"Invalid is_representative for {mag}")

        fasta = root / "hq/genomes" / f"{mag}.fa"
        counts = source_counts(fasta, sample_project)
        samples = sorted(counts)
        projects = sorted({sample_project[sample] for sample in samples})
        alias = "NA"
        if len(samples) == 1:
            prefix, marker, number = mag.rpartition("__MAG_")
            if not prefix or not marker or not number.isdigit():
                raise ValueError(f"Unexpected MAG ID format: {mag}")
            alias = f"{samples[0]}_MAG_{number}"
            if alias in aliases_seen:
                raise ValueError(f"Sample-first MAG alias is not unique: {alias}")
            aliases_seen.add(alias)
        alias_by_mag[mag] = alias
        sources_by_mag[mag] = set(samples)
        members_by_sgb[sgb].append(mag)
        mag_rows.append({
            "sgb_id": sgb,
            "mag_id": mag,
            "sample_first_mag_alias": alias,
            "is_representative": row["is_representative"],
            "group": meta["group"],
            "completeness_percent": meta["completeness"],
            "contamination_percent": meta["contamination"],
            "assembly_sample_count": str(len(samples)),
            "assembly_sample_ids": ";".join(samples),
            "assembly_sample_groups": ";".join(sorted({sample_group[s] for s in samples})),
            "assembly_project_ids": ";".join(projects),
            "assembly_sample_project_pairs": ";".join(
                f"{s}:{sample_project[s]}" for s in samples
            ),
            "contigs_by_assembly_sample": ";".join(f"{s}:{counts[s]}" for s in samples),
            "contig_count": str(sum(counts.values())),
        })

    if set(members_by_sgb) != set(reps):
        raise ValueError("SGB sets differ between membership and representative tables")
    sgb_rows: list[dict[str, str]] = []
    for sgb in sorted(reps):
        rep = reps[sgb]
        representative_mag = rep["representative_mag"]
        if flagged_reps[sgb] != [representative_mag]:
            raise ValueError(f"SGB {sgb} has inconsistent representative flags")
        meta = hq[representative_mag]
        for key in ("completeness", "contamination"):
            if abs(float(rep[key]) - float(meta[key])) > 1e-5:
                raise ValueError(f"SGB {sgb} representative {key} differs from HQ manifest")
        rep_samples = sources_by_mag[representative_mag]
        member_samples = set().union(*(sources_by_mag[m] for m in members_by_sgb[sgb]))
        sgb_rows.append({
            "sgb_id": sgb,
            "mag_member_count": str(len(members_by_sgb[sgb])),
            "representative_mag": representative_mag,
            "representative_sample_first_alias": alias_by_mag[representative_mag],
            "representative_group": meta["group"],
            "member_groups": ";".join(sorted({hq[m]["group"] for m in members_by_sgb[sgb]})),
            "representative_completeness_percent": rep["completeness"],
            "representative_contamination_percent": rep["contamination"],
            "representative_assembly_sample_ids": ";".join(sorted(rep_samples)),
            "representative_assembly_sample_groups": ";".join(
                sorted({sample_group[s] for s in rep_samples})
            ),
            "representative_assembly_project_ids": ";".join(
                sorted({sample_project[s] for s in rep_samples})
            ),
            "representative_assembly_sample_project_pairs": ";".join(
                f"{s}:{sample_project[s]}" for s in sorted(rep_samples)
            ),
            "member_assembly_sample_count": str(len(member_samples)),
            "member_assembly_sample_ids": ";".join(sorted(member_samples)),
            "member_assembly_sample_groups": ";".join(
                sorted({sample_group[s] for s in member_samples})
            ),
            "member_assembly_project_ids": ";".join(
                sorted({sample_project[s] for s in member_samples})
            ),
            "member_assembly_sample_project_pairs": ";".join(
                f"{s}:{sample_project[s]}" for s in sorted(member_samples)
            ),
        })

    args.out.mkdir(parents=True, exist_ok=True)
    mag_rows.sort(key=lambda row: (row["sgb_id"], row["mag_id"]))
    write_tsv(args.out / "mag_source_metadata.tsv", list(mag_rows[0]), mag_rows)
    write_tsv(args.out / "sgb_summary.tsv", list(sgb_rows[0]), sgb_rows)
    print(f"SGBs={len(sgb_rows)} MAGs={len(mag_rows)}")
    print(args.out / "sgb_summary.tsv")
    print(args.out / "mag_source_metadata.tsv")


if __name__ == "__main__":
    main()
