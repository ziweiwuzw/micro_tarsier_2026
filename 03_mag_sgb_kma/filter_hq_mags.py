#!/usr/bin/env python3
"""Filter candidate MAGs using CheckM2 completeness and contamination."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
from typing import Dict, Iterable


def norm(value: str) -> str:
    name = Path(value.strip()).name
    lower = name.lower()
    for suffix in (".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa"):
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    return name


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return [{k: (v or "").strip() for k, v in row.items()} for row in csv.DictReader(handle, delimiter="\t")]


def find_column(row: dict[str, str], names: Iterable[str]) -> str:
    lookup = {k.casefold(): k for k in row}
    for name in names:
        if name.casefold() in lookup:
            return lookup[name.casefold()]
    raise KeyError(f"None of the columns found: {', '.join(names)}; have {', '.join(row)}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--checkm2", required=True, type=Path)
    parser.add_argument("--min-completeness", type=float, default=90.0)
    parser.add_argument("--max-contamination", type=float, default=5.0)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args()

    manifest = args.output_root / "hq_manifest.tsv"
    genomes_dir = args.output_root / "genomes"
    # dRep's --genomeInfo parser requires comma-separated CSV, not TSV.
    drep_info = args.output_root / "drep_genome_info.csv"
    done = args.output_root / ".filter.done"
    expected_config = (
        f"min_completeness={args.min_completeness}\n"
        f"max_contamination={args.max_contamination}\n"
        "quality_tool=CheckM2\n"
        "gunc=NOT_RUN\n"
    )
    if manifest.is_file() and manifest.stat().st_size > 0:
        if done.is_file() and done.read_text(encoding="utf-8") == expected_config:
            print(f"[REUSE] {manifest}")
            return
        raise SystemExit(
            f"Existing HQ results were produced with a different or unknown filter configuration: {args.output_root}. "
            "Use a new --output-root/PIPELINE_ROOT."
        )
    if genomes_dir.exists() and any(genomes_dir.iterdir()):
        raise SystemExit(f"HQ directory exists without completed manifest: {genomes_dir}")
    genomes_dir.mkdir(parents=True, exist_ok=True)

    candidates = read_tsv(args.candidates)
    checkm_rows = read_tsv(args.checkm2)
    if not checkm_rows:
        raise SystemExit("CheckM2 report has no data rows")
    name_col = find_column(checkm_rows[0], ["Name", "genome", "Bin Id"])
    comp_col = find_column(checkm_rows[0], ["Completeness"])
    cont_col = find_column(checkm_rows[0], ["Contamination"])
    quality: Dict[str, tuple[float, float]] = {}
    for row in checkm_rows:
        key = norm(row[name_col])
        if key in quality:
            raise SystemExit(f"Duplicate CheckM2 genome: {key}")
        quality[key] = (float(row[comp_col]), float(row[cont_col]))

    output_rows = []
    info_rows = []
    for candidate in candidates:
        mag_id = candidate["mag_id"]
        key = norm(candidate["filename"])
        if key not in quality:
            raise SystemExit(f"Candidate missing from CheckM2: {candidate['filename']}")
        completeness, contamination = quality[key]
        quality_ok = completeness >= args.min_completeness and contamination <= args.max_contamination
        selected = quality_ok
        reason = "PASS" if selected else ";".join(
            part for part, failed in (
                ("LOW_COMPLETENESS", completeness < args.min_completeness),
                ("HIGH_CONTAMINATION", contamination > args.max_contamination),
            ) if failed
        )
        row = dict(candidate)
        row.update({
            "completeness": f"{completeness:.6g}",
            "contamination": f"{contamination:.6g}",
            "gunc_pass": "NOT_RUN",
            "selected_hq": str(selected),
            "filter_reason": reason,
        })
        output_rows.append(row)
        if selected:
            source = Path(candidate["candidate_path"])
            target = genomes_dir / candidate["filename"]
            target.symlink_to(source.resolve())
            info_rows.append({
                "genome": candidate["filename"],
                "completeness": f"{completeness:.6g}",
                "contamination": f"{contamination:.6g}",
            })

    if not info_rows:
        raise SystemExit("No MAG passed the requested quality filters")
    tmp = manifest.with_name(manifest.name + f".part.{os.getpid()}")
    with tmp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    tmp.replace(manifest)
    with drep_info.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["genome", "completeness", "contamination"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(info_rows)
    done.write_text(expected_config, encoding="utf-8")
    print(f"[SUCCESS] selected={len(info_rows)}/{len(output_rows)} hq_dir={genomes_dir}")


if __name__ == "__main__":
    main()
