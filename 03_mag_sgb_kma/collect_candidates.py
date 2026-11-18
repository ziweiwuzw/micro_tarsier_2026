#!/usr/bin/env python3
"""Collect completed bins; optionally use Pan_troglodytes MetaWRAP output."""

from __future__ import annotations

import argparse
import csv
import os
import re
from pathlib import Path


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip()).strip("_") or "unnamed"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--groups", required=True, type=Path)
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--expected-groups", type=int, default=32)
    parser.add_argument("--completeness", type=int, default=50)
    parser.add_argument("--contamination", type=int, default=10)
    parser.add_argument(
        "--pan-source", choices=("semibin2", "metawrap"),
        default=os.environ.get("PAN_SOURCE", "semibin2"),
    )
    args = parser.parse_args()

    genomes = args.output_root / "genomes"
    manifest = args.output_root / "candidate_manifest.tsv"
    if manifest.is_file() and manifest.stat().st_size > 0:
        print(f"[REUSE] {manifest}")
        return
    if genomes.exists() and any(genomes.iterdir()):
        raise SystemExit(f"Candidate directory exists without completed manifest: {genomes}")
    genomes.mkdir(parents=True, exist_ok=True)

    with args.groups.open(newline="", encoding="utf-8") as handle:
        group_rows = list(csv.DictReader(handle, delimiter="\t"))
    slugs = [safe_name(row["group_slug"]) for row in group_rows]
    if len(slugs) != args.expected_groups or len(set(slugs)) != len(slugs):
        raise SystemExit(f"Expected {args.expected_groups} unique groups, got {len(slugs)} rows")
    if slugs.count("Pan_troglodytes") != 1:
        raise SystemExit("Expected exactly one Pan_troglodytes group")
    output_rows = []
    for group_row in group_rows:
        group_slug = safe_name(group_row["group_slug"])
        group_output = args.source_root / "groups" / group_slug
        if not group_output.is_dir():
            raise SystemExit(f"Missing source group: {group_output}")
        if group_slug == "Pan_troglodytes" and args.pan_source == "semibin2":
            bins = group_output / "inputs" / "semibin2"
            source_binner = "semibin2"
        else:
            refinement = group_output / "refinement"
            done = refinement / ".metawrap_refinement.done"
            bins = refinement / f"metawrap_{args.completeness}_{args.contamination}_bins"
            if not done.is_file() or done.stat().st_size == 0 or not bins.is_dir():
                raise SystemExit(f"Incomplete refinement for {group_slug}: {refinement}")
            source_binner = "metawrap"
        bin_files = sorted(bins.glob("*.fa"))
        if not bin_files:
            raise SystemExit(f"No refined FASTA bins: {bins}")
        for index, source in enumerate(bin_files, start=1):
            if not source.is_file() or source.stat().st_size == 0:
                raise SystemExit(f"Empty or broken bin FASTA: {source}")
            mag_id = f"{group_slug}__MAG_{index:06d}"
            filename = mag_id + ".fa"
            target = genomes / filename
            if target.exists() or target.is_symlink():
                raise SystemExit(f"Duplicate candidate target: {target}")
            target.symlink_to(source.resolve())
            output_rows.append({
                "mag_id": mag_id,
                "filename": filename,
                "group": group_row["group"],
                "group_slug": group_slug,
                "source_binner": source_binner,
                "source_bin": str(source),
                "candidate_path": str(target),
            })

    tmp = manifest.with_name(manifest.name + f".part.{os.getpid()}")
    with tmp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    tmp.replace(manifest)
    print(f"[SUCCESS] candidates={len(output_rows)} manifest={manifest}")


if __name__ == "__main__":
    main()
