#!/usr/bin/env python3
"""Build the shared KMA SGB contig FASTA and one-to-one contig mapping."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
from pathlib import Path


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--representatives", required=True, type=Path)
    parser.add_argument("--reference-root", required=True, type=Path)
    args = parser.parse_args()
    with args.representatives.open(newline="", encoding="utf-8") as handle:
        reps = list(csv.DictReader(handle, delimiter="\t"))
    if not reps or not {"sgb_id", "representative_path"} <= set(reps[0]):
        raise SystemExit("Missing or invalid SGB representatives table")
    ids = [row["sgb_id"] for row in reps]
    if len(set(ids)) != len(ids):
        raise SystemExit("Duplicate SGB ID")
    out = args.reference_root
    if out.exists():
        raise SystemExit(f"Reference directory exists; inspect before rebuilding: {out}")
    out.mkdir(parents=True)
    catalogue = out / "catalogue.fa"
    mapping = out / "contig_to_sgb.tsv"
    count = 0
    with catalogue.open("w", encoding="utf-8") as fasta, mapping.open("w", newline="", encoding="utf-8") as tab:
        writer = csv.writer(tab, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig_id", "sgb_id", "original_id", "length"])
        for row in reps:
            sgb = row["sgb_id"]
            source = Path(row["representative_path"])
            if not source.is_file() or source.stat().st_size == 0:
                raise SystemExit(f"Missing SGB FASTA: {source}")
            number = 0
            header = None
            sequence: list[str] = []

            def emit() -> None:
                nonlocal number, count
                if header is None:
                    return
                bases = "".join(sequence).upper()
                if not bases or any(letter not in "ACGTURYSWKMBDHVN" for letter in bases):
                    raise ValueError(f"Invalid sequence in {source}: {header}")
                number += 1
                count += 1
                contig = f"{sgb}_contig_{number}"
                fasta.write(f">{contig}\n{bases}\n")
                writer.writerow([contig, sgb, header, len(bases)])

            with source.open(encoding="utf-8") as handle:
                for line in handle:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        emit()
                        header = line[1:].strip()
                        if not header:
                            raise ValueError(f"Empty FASTA header in {source}")
                        sequence = []
                    elif header is None:
                        raise ValueError(f"Sequence before header in {source}")
                    else:
                        sequence.append(line)
            emit()
            if number == 0:
                raise ValueError(f"No contigs in {source}")
    (out / ".reference.done").write_text(
        f"representatives_sha256={digest(args.representatives)}\n"
        f"catalogue_sha256={digest(catalogue)}\n"
        f"map_sha256={digest(mapping)}\n"
        f"sgbs={len(ids)}\ncontigs={count}\n",
        encoding="utf-8",
    )
    print(f"[SUCCESS] sgbs={len(ids)} contigs={count} reference={out}")


if __name__ == "__main__":
    main()
