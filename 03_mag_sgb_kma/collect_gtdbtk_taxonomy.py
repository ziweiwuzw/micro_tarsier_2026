#!/usr/bin/env python3
"""Combine GTDB-Tk bacterial/archaeal summaries into one SGB taxonomy table."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path


RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")


def genome_id(value: str) -> str:
    name = Path(value.strip()).name
    for suffix in (".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa"):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return name


def split_taxonomy(value: str) -> dict[str, str]:
    parsed = {rank: "" for rank in RANKS}
    prefixes = {
        "d__": "domain", "p__": "phylum", "c__": "class",
        "o__": "order", "f__": "family", "g__": "genus", "s__": "species",
    }
    for item in value.split(";"):
        item = item.strip()
        for prefix, rank in prefixes.items():
            if item.startswith(prefix):
                parsed[rank] = item
                break
    return parsed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--representatives", required=True, type=Path)
    parser.add_argument("--summary", required=True, action="append", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--method", default="GTDB-Tk-v2.4.0-R220")
    args = parser.parse_args()

    with args.representatives.open(newline="", encoding="utf-8") as handle:
        representatives = list(csv.DictReader(handle, delimiter="\t"))
    if not representatives or "sgb_id" not in representatives[0]:
        raise SystemExit(f"Unexpected or empty representative table: {args.representatives}")
    expected = [row["sgb_id"] for row in representatives]

    by_sgb: dict[str, dict[str, str]] = {}
    for summary in args.summary:
        if not summary.is_file() or summary.stat().st_size == 0:
            continue
        with summary.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        for row in rows:
            if "user_genome" not in row or "classification" not in row:
                raise SystemExit(f"Unexpected GTDB-Tk summary columns: {summary}")
            sgb_id = genome_id(row["user_genome"])
            if sgb_id in by_sgb:
                raise SystemExit(f"Duplicate GTDB-Tk classification for {sgb_id}")
            by_sgb[sgb_id] = row

    missing = [sgb_id for sgb_id in expected if sgb_id not in by_sgb]
    unexpected = sorted(set(by_sgb) - set(expected))
    if missing or unexpected:
        raise SystemExit(
            "GTDB-Tk/representative mismatch; "
            f"missing={','.join(missing) or 'NONE'}; "
            f"unexpected={','.join(unexpected) or 'NONE'}"
        )

    output_rows = []
    for sgb_id in expected:
        source = by_sgb[sgb_id]
        taxonomy = source["classification"].strip()
        row = {
            "sgb_id": sgb_id,
            "taxonomy_method": args.method,
            "taxonomy": taxonomy,
            **split_taxonomy(taxonomy),
            "classification_method": source.get("classification_method", ""),
            "fastani_reference": source.get("fastani_reference", ""),
            "fastani_ani": source.get("fastani_ani", ""),
            "fastani_af": source.get("fastani_af", ""),
            "closest_genome_reference": source.get("closest_genome_reference", ""),
            "closest_genome_ani": source.get("closest_genome_ani", ""),
            "closest_genome_af": source.get("closest_genome_af", ""),
            "note": source.get("note", ""),
            "warnings": source.get("warnings", ""),
        }
        output_rows.append(row)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(args.output.name + f".part.{os.getpid()}")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(output_rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(output_rows)
    temporary.replace(args.output)
    print(f"[SUCCESS] GTDB-Tk taxonomy rows={len(output_rows)} output={args.output}")


if __name__ == "__main__":
    main()
