#!/usr/bin/env python3
"""Convert the master read manifest to one KMA row per sample (lanes stay paired)."""

from __future__ import annotations

import argparse
import csv
import os
import re
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-samples", type=int, default=547)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"Reads manifest already exists; inspect before reuse: {args.output}")
    with args.input.open(newline="", encoding="utf-8") as handle:
        source = list(csv.DictReader(handle, delimiter="\t"))
    required = {"sample_id", "species", "r1", "r2", "mean_read_length"}
    if not source or not required <= set(source[0]):
        raise SystemExit(f"Missing required fields in {args.input}: {sorted(required)}")
    samples: dict[str, dict] = {}
    for line, row in enumerate(source, start=2):
        sid, species = row["sample_id"].strip(), row["species"].strip()
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", sid) or not re.fullmatch(r"[A-Za-z0-9_.-]+", species):
            raise SystemExit(f"Unsafe sample or species ID at line {line}")
        r1, r2 = Path(row["r1"].strip()), Path(row["r2"].strip())
        if not r1.is_file() or not r2.is_file() or not r1.stat().st_size or not r2.stat().st_size:
            raise SystemExit(f"Missing or empty FASTQ at line {line}: {r1}, {r2}")
        length = float(row["mean_read_length"])
        if not 0 < length <= 1000:
            raise SystemExit(f"Invalid mean read length at line {line}: {length}")
        item = samples.setdefault(sid, {
            "sample_id": sid, "group": species, "project": row.get("project", ""),
            "target": row.get("target", ""), "read1": [], "read2": [], "lengths": [],
        })
        if item["group"] != species or str(r1) in item["read1"] or str(r2) in item["read2"]:
            raise SystemExit(f"Inconsistent or duplicate lane for {sid}")
        item["read1"].append(str(r1))
        item["read2"].append(str(r2))
        item["lengths"].append(length)
    if len(samples) != args.expected_samples:
        raise SystemExit(f"Expected {args.expected_samples} samples, found {len(samples)}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(args.output.name + f".part.{os.getpid()}")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "sample_index", "sample_id", "group", "project", "target", "read1", "read2",
            "mean_read_length", "coverage_bases",
        ], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for index, sid in enumerate(sorted(samples), start=1):
            item = samples[sid]
            mean_length = sum(item["lengths"]) / len(item["lengths"])
            writer.writerow({
                "sample_index": index, "sample_id": sid, "group": item["group"],
                "project": item["project"], "target": item["target"],
                "read1": ";".join(item["read1"]), "read2": ";".join(item["read2"]),
                "mean_read_length": f"{mean_length:.8g}",
                "coverage_bases": f"{2 * mean_length:.8g}",
            })
    temporary.replace(args.output)
    print(f"[SUCCESS] samples={len(samples)} reads={args.output}")


if __name__ == "__main__":
    main()
