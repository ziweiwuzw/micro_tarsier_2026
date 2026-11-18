#!/usr/bin/env python3
"""Validate VFDB/Toxinome references and retain source-specific protein hits."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path


HIT_FIELDS = ("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qstart", "qend",
              "sstart", "send", "evalue", "bitscore")
OUT_FIELDS = ("mag_id", "orf_id", "source", "reference_id", "class", "description", "reference_source",
              "pident", "query_coverage", "subject_coverage", "evalue", "bitscore")
SOURCES = ("vfdb_A", "vfdb_B", "toxinome")


def fasta_headers(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                # Keep valid UTF-8 and preserve isolated legacy bytes in VFDB descriptions.
                decoded = line[1:].decode("utf-8", errors="surrogateescape")
                header = "".join(chr(ord(char) - 0xDC00) if "\udc80" <= char <= "\udcff" else char
                                 for char in decoded).strip()
                if not header:
                    raise ValueError(f"Empty FASTA header: {path}")
                protein_id = header.split()[0]
                if not protein_id.isascii():
                    raise ValueError(f"Non-ASCII FASTA protein ID in {path}: {protein_id!r}")
                yield protein_id, header


def toxicity(value: str) -> str:
    normalized = value.strip().casefold().replace("_", "").replace("-", "").replace(" ", "")
    if normalized in ("tox", "toxin", "toxic"):
        return "toxin"
    if normalized in ("antitox", "antitoxin", "antitoxic", "immunity", "immune"):
        return "immunity"
    return "unclassified"


def prepare(args: argparse.Namespace) -> None:
    for path in (args.vfdb_a, args.vfdb_b, args.toxin_fasta, args.toxin_csv):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing reference file: {path}")
    args.out.mkdir(parents=True, exist_ok=True)
    db_path = args.out / "catalog.sqlite"
    if db_path.exists():
        raise ValueError(f"Reference catalog exists; use a new ANN_OUT: {db_path}")
    connection = sqlite3.connect(db_path)
    try:
        connection.execute("CREATE TABLE refs (source TEXT NOT NULL, protein_id TEXT NOT NULL, class TEXT NOT NULL, "
                           "description TEXT NOT NULL, reference_source TEXT NOT NULL, PRIMARY KEY(source, protein_id))")
        counts = Counter()
        for source, path in (("vfdb_A", args.vfdb_a), ("vfdb_B", args.vfdb_b)):
            with connection:
                for protein_id, header in fasta_headers(path):
                    connection.execute("INSERT INTO refs VALUES (?, ?, ?, ?, ?)",
                                       (source, protein_id, "virulence_candidate", header, "VFDB"))
                    counts[source] += 1
        toxin_ids = set()
        with connection:
            for protein_id, header in fasta_headers(args.toxin_fasta):
                if protein_id in toxin_ids:
                    raise ValueError(f"Duplicate Toxinome FASTA ID: {protein_id}")
                toxin_ids.add(protein_id)
        metadata_ids = set()
        with args.toxin_csv.open(newline="", encoding="utf-8-sig") as handle, connection:
            reader = csv.DictReader(handle)
            required = {"proteinID", "protein_toxicity", "product_name", "Data Base name"}
            if not required <= set(reader.fieldnames or []):
                raise ValueError(f"Unexpected Toxinome columns: {reader.fieldnames}")
            for row in reader:
                protein_id = row["proteinID"].strip()
                if protein_id not in toxin_ids:
                    continue
                if protein_id in metadata_ids:
                    raise ValueError(f"Duplicate Toxinome proteinID: {protein_id}")
                metadata_ids.add(protein_id)
                category = toxicity(row["protein_toxicity"])
                connection.execute("INSERT INTO refs VALUES (?, ?, ?, ?, ?)",
                                   ("toxinome", protein_id, category, row["product_name"].strip(),
                                    row["Data Base name"].strip()))
                counts["toxinome_" + category] += 1
        if not toxin_ids or len(metadata_ids) != len(toxin_ids):
            raise ValueError(f"Toxinome FASTA/CSV ID mismatch: fasta={len(toxin_ids)} matched={len(metadata_ids)}")
        if not counts["vfdb_A"] or not counts["vfdb_B"]:
            raise ValueError("Empty VFDB FASTA")
        connection.execute("CREATE INDEX refs_class_idx ON refs(source, class)")
        connection.commit()
        with (args.out / "catalog_counts.tsv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(("category", "proteins"))
            writer.writerows(sorted(counts.items()))
    except Exception:
        connection.close()
        db_path.unlink(missing_ok=True)
        raise
    connection.close()
    print(f"Reference catalog: {db_path}; counts={dict(counts)}")


def parse(args: argparse.Namespace) -> None:
    connection = sqlite3.connect(args.catalog)
    all_hits = []
    best: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for source in SOURCES:
        path = args.raw_dir / f"{source}.tsv"
        if not path.is_file():
            raise ValueError(f"Missing DIAMOND output: {path}")
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.reader(handle, delimiter="\t")
            for values in reader:
                if len(values) != len(HIT_FIELDS):
                    raise ValueError(f"Unexpected DIAMOND row in {path}: {values}")
                hit = dict(zip(HIT_FIELDS, values))
                orf = hit["qseqid"]
                if not orf.startswith(args.mag + "_contig_"):
                    raise ValueError(f"Foreign ORF in {path}: {orf}")
                protein_id = hit["sseqid"]
                meta = connection.execute("SELECT class, description, reference_source FROM refs "
                                          "WHERE source = ? AND protein_id = ?", (source, protein_id)).fetchone()
                if meta is None:
                    raise ValueError(f"Reference ID absent from catalog: {source}:{protein_id}")
                qlen, slen, length = float(hit["qlen"]), float(hit["slen"]), float(hit["length"])
                pident, evalue, bitscore = float(hit["pident"]), float(hit["evalue"]), float(hit["bitscore"])
                if min(qlen, slen, length) <= 0 or not all(map(math.isfinite, (pident, evalue, bitscore))):
                    raise ValueError(f"Invalid DIAMOND score in {path}: {values}")
                qcov = 100 * (abs(int(hit["qend"]) - int(hit["qstart"])) + 1) / qlen
                scov = 100 * (abs(int(hit["send"]) - int(hit["sstart"])) + 1) / slen
                row = {"mag_id": args.mag, "orf_id": orf, "source": source, "reference_id": protein_id,
                       "class": meta[0], "description": meta[1], "reference_source": meta[2],
                       "pident": f"{pident:.4g}", "query_coverage": f"{qcov:.4g}",
                       "subject_coverage": f"{scov:.4g}", "evalue": f"{evalue:.4g}",
                       "bitscore": f"{bitscore:.4g}"}
                all_hits.append(row)
                best[orf, source].append(row)
    args.out.mkdir(parents=True, exist_ok=True)
    def write(name: str, data: list[dict]) -> None:
        with (args.out / name).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=OUT_FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(data)
    write("all_hits.tsv", all_hits)
    for source in SOURCES:
        write(f"{source}_hits.tsv", [row for row in all_hits if row["source"] == source])
    selected = []
    for group in best.values():
        group.sort(key=lambda r: (-float(r["bitscore"]), float(r["evalue"]), r["reference_id"]))
        top = dict(group[0])
        if top["source"] == "toxinome":
            rival = next((r for r in group[1:] if r["class"] != top["class"]), None)
            if rival and float(rival["bitscore"]) >= 0.9 * float(top["bitscore"]):
                top["class"] = "ambiguous"
        selected.append(top)
    selected.sort(key=lambda r: (r["mag_id"], r["orf_id"], r["source"]))
    write("best_hits.tsv", selected)
    for source in SOURCES:
        write(f"{source}_best_hits.tsv", [row for row in selected if row["source"] == source])
    print(f"MAG={args.mag} all_hits={len(all_hits)} best_hits={len(selected)}")


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--vfdb-a", type=Path, required=True)
    p.add_argument("--vfdb-b", type=Path, required=True)
    p.add_argument("--toxin-fasta", type=Path, required=True)
    p.add_argument("--toxin-csv", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p = sub.add_parser("parse")
    p.add_argument("--mag", required=True)
    p.add_argument("--catalog", type=Path, required=True)
    p.add_argument("--raw-dir", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    {"prepare": prepare, "parse": parse}[args.command](args)


if __name__ == "__main__":
    main()
