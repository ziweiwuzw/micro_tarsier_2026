#!/usr/bin/env python3
"""Freeze an existing SGB set, call shared ORFs, and summarize dbCAN families."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


FIELDS = ("mag_id", "sgb_id", "is_representative", "genome_path", "genome_sha256")
FAMILY = re.compile(r"(?<![A-Za-z0-9])(?:GH|GT|PL|CE|AA|CBM)\d+(?:_\d+)?(?![A-Za-z0-9])")


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        result = list(reader)
        if not reader.fieldnames or not result:
            raise ValueError(f"Empty table: {path}")
        return result


def write_table(path: Path, fields: tuple[str, ...] | list[str], data: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", dir=path.parent, prefix=path.name + ".", delete=False,
                                     newline="", encoding="utf-8") as handle:
        temporary = Path(handle.name)
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(data)
    temporary.replace(path)


def check_abundance_columns(columns: list[str] | None, sgb_order: list[str], out: Path) -> None:
    if not columns or columns[0] != "sample_id" or len(columns) != len(set(columns)):
        raise ValueError("Salmon SGB columns do not match frozen manifest")
    available = set(columns[1:])
    expected = set(sgb_order)
    if not expected <= available or (not (out / "selected_sgbs.tsv").exists() and expected != available):
        raise ValueError("Salmon SGB columns do not match frozen manifest")


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def prepare(args: argparse.Namespace) -> None:
    root = args.pipeline_root
    memberships = rows(root / "sgb/sgb_membership.tsv")
    representatives = rows(root / "sgb/sgb_representatives.tsv")
    hq = root / "hq/genomes"
    required_members = {"mag_id", "sgb_id", "is_representative"}
    if not required_members <= memberships[0].keys():
        raise ValueError("Unexpected SGB membership columns")
    if not {"sgb_id", "representative_mag"} <= representatives[0].keys():
        raise ValueError("Unexpected SGB representative columns")
    selection = None
    if getattr(args, "sgb_list", None):
        selected = rows(args.sgb_list)
        if "sgb_id" not in selected[0]:
            raise ValueError("SGB selection requires an sgb_id column")
        selection = {r["sgb_id"] for r in selected}
        if len(selection) != len(selected) or not selection:
            raise ValueError("Duplicate or empty SGB selection")
    by_sgb = {}
    for row in representatives:
        sgb = row["sgb_id"]
        if sgb in by_sgb:
            raise ValueError(f"Duplicate representative: {sgb}")
        if selection is None or sgb in selection:
            rep_file = root / "sgb/representatives" / f"{sgb}.fa"
            if not rep_file.is_file() or rep_file.stat().st_size == 0:
                raise ValueError(f"Missing SGB representative: {rep_file}")
        by_sgb[sgb] = row["representative_mag"]
    if selection is not None and not selection <= set(by_sgb):
        raise ValueError("Unknown SGB in selection")
    expected = selection if selection is not None else set(by_sgb)
    manifest = []
    seen = set()
    reps_seen = set()
    for row in memberships:
        mag, sgb = row["mag_id"], row["sgb_id"]
        if sgb not in expected:
            continue
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", mag) or not re.fullmatch(r"SGB_\d+", sgb):
            raise ValueError(f"Unsafe MAG/SGB ID: {mag}, {sgb}")
        if mag in seen or sgb not in expected:
            raise ValueError(f"Duplicate MAG or unknown SGB: {mag}, {sgb}")
        seen.add(mag)
        rep = mag == by_sgb[sgb]
        if row["is_representative"].lower() != str(rep).lower():
            raise ValueError(f"Representative flag mismatch: {mag}")
        genome = hq / f"{mag}.fa"
        if not genome.is_file() or genome.stat().st_size == 0:
            raise ValueError(f"Missing HQ member MAG: {genome}")
        genome_sha256 = digest(genome)
        if rep:
            reps_seen.add(sgb)
            sgb_rep = root / "sgb/representatives" / f"{sgb}.fa"
            if genome_sha256 != digest(sgb_rep):
                raise ValueError(f"Representative differs from SGB reference: {sgb}")
        manifest.append(dict(zip(FIELDS, (mag, sgb, str(rep).lower(), str(genome.resolve()), genome_sha256))))
    if reps_seen != expected:
        raise ValueError(f"Missing representative members: {sorted(expected - reps_seen)}")
    manifest.sort(key=lambda r: (r["sgb_id"], r["mag_id"]))
    selection_path = args.out / "selected_sgbs.tsv"
    if selection is not None:
        selected_rows = [{"sgb_id": sgb} for sgb in sorted(selection)]
        if selection_path.exists():
            if rows(selection_path) != selected_rows:
                raise ValueError(f"SGB selection changed; choose a new ANN_OUT: {selection_path}")
        elif (args.out / "input_manifest.tsv").exists():
            raise ValueError(f"Existing ANN_OUT has no SGB selection: {args.out}")
        else:
            write_table(selection_path, ["sgb_id"], selected_rows)
    elif selection_path.exists():
        raise ValueError(f"Existing ANN_OUT is restricted to selected SGBs: {args.out}")
    output = args.out / "input_manifest.tsv"
    if output.exists():
        if rows(output) != manifest:
            raise ValueError(f"Input set changed; choose a new ANN_OUT: {output}")
    else:
        write_table(output, FIELDS, manifest)
    print(f"SGBs={len(expected)} MAGs={len(manifest)} manifest={output}")


def normalize_fasta(source: Path, target: Path, mag: str, map_path: Path) -> None:
    contigs = []
    original = None
    sequence = []

    def flush(handle) -> None:
        if original is None:
            return
        seq = "".join(sequence).upper()
        if not seq or re.search(r"[^ACGTURYSWKMBDHVN]", seq):
            raise ValueError(f"Invalid sequence in {source}: {original}")
        contig = f"{mag}_contig_{len(contigs) + 1:06d}"
        handle.write(f">{contig}\n{seq}\n")
        contigs.append({"contig_id": contig, "original_header": original, "length": len(seq)})

    with source.open(encoding="utf-8") as src, target.open("w", encoding="utf-8") as dest:
        for line in src:
            line = line.strip()
            if line.startswith(">"):
                flush(dest)
                original = line[1:]
                if not original:
                    raise ValueError(f"Empty FASTA header in {source}")
                sequence = []
            elif line:
                if original is None:
                    raise ValueError(f"Sequence before FASTA header in {source}")
                sequence.append(line)
        flush(dest)
    if not contigs:
        raise ValueError(f"No contigs in {source}")
    write_table(map_path, ["contig_id", "original_header", "length"], contigs)


def canonicalize_gff(gff: Path, faa: Path) -> set[str]:
    prodigal_to_orf = {}
    with faa.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            orf = line[1:].split()[0]
            match = re.search(r"(?:^|[;\s])ID=([^;\s]+)", line)
            if not match or match.group(1) in prodigal_to_orf:
                raise ValueError(f"Missing or repeated Prodigal ID in {faa}: {line.strip()}")
            prodigal_to_orf[match.group(1)] = orf
    seen = set()
    temp = gff.with_name(gff.name + ".canonical")
    with gff.open(encoding="utf-8") as source, temp.open("w", encoding="utf-8") as dest:
        for line in source:
            if line.startswith("#") or not line.strip():
                dest.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                raise ValueError(f"Unexpected Prodigal GFF record: {line.strip()}")
            match = re.search(r"(?:^|;)ID=([^;]+)", fields[8])
            if not match or match.group(1) not in prodigal_to_orf:
                raise ValueError(f"GFF/FAA ID mismatch: {line.strip()}")
            orf = prodigal_to_orf[match.group(1)]
            if orf in seen:
                raise ValueError(f"Repeated ORF in GFF: {orf}")
            seen.add(orf)
            fields[8] = re.sub(r"(^|;)ID=[^;]+", lambda m: m.group(1) + f"ID={orf}", fields[8], count=1)
            dest.write("\t".join(fields) + "\n")
    if seen != set(prodigal_to_orf.values()):
        raise ValueError(f"GFF lacks {len(set(prodigal_to_orf.values()) - seen)} FAA genes")
    temp.replace(gff)
    return seen


def genes(args: argparse.Namespace) -> None:
    manifest = {r["mag_id"]: r for r in rows(args.out / "input_manifest.tsv")}
    if args.mag not in manifest:
        raise ValueError(f"MAG absent from frozen manifest: {args.mag}")
    row = manifest[args.mag]
    if digest(Path(row["genome_path"])) != row["genome_sha256"]:
        raise ValueError(f"Genome changed after manifest was frozen: {args.mag}")
    final = args.out / "01_orf" / args.mag
    if (final / ".done").exists():
        for suffix in ("fa", "faa", "fna", "gff", "contigs.tsv"):
            if not (final / f"{args.mag}.{suffix}").is_file():
                raise ValueError(f"Incomplete ORF output: {final}")
        print(f"REUSE {final}")
        return
    if final.exists():
        raise ValueError(f"Partial ORF output: {final}")
    final.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=args.mag + ".", dir=final.parent))
    try:
        path = lambda suffix: temporary / f"{args.mag}.{suffix}"
        normalize_fasta(Path(row["genome_path"]), path("fa"), args.mag,
                        temporary / f"{args.mag}.contigs.tsv")
        subprocess.run(["prodigal", "-i", str(path("fa")), "-a", str(path("faa")),
                        "-d", str(path("fna")), "-f", "gff", "-o", str(path("gff")),
                        "-p", "meta", "-g", "11", "-q"], check=True)
        for suffix in ("faa", "fna", "gff"):
            output_path = path(suffix)
            if not output_path.is_file() or output_path.stat().st_size == 0:
                raise ValueError(f"Missing Prodigal output: {output_path}")
        faa_ids = canonicalize_gff(path("gff"), path("faa"))
        fna_ids = {line.split()[0][1:] for line in path("fna").open() if line.startswith(">")}
        if not faa_ids or faa_ids != fna_ids or not all(x.startswith(args.mag + "_contig_") for x in faa_ids):
            raise ValueError(f"Inconsistent Prodigal IDs: {args.mag}")
        (temporary / ".done").write_text(f"genes={len(faa_ids)}\n", encoding="utf-8")
        temporary.rename(final)
    except Exception:
        shutil.rmtree(temporary)
        raise
    print(f"ORFs={len(faa_ids)} {final}")


def summarize(args: argparse.Namespace) -> None:
    manifest = rows(args.out / "input_manifest.tsv")
    sgb_order = list(dict.fromkeys(r["sgb_id"] for r in manifest))
    by_sgb: dict[str, list[str]] = defaultdict(list)
    mag_families: dict[str, dict[str, set[str]]] = {}
    all_families = set()
    for row in manifest:
        mag, sgb = row["mag_id"], row["sgb_id"]
        by_sgb[sgb].append(mag)
        overview = args.out / "04_dbcan" / mag / "overview.tsv"
        if not (args.out / "04_dbcan" / mag / ".done").exists():
            raise ValueError(f"Missing completed dbCAN result: {mag}")
        with overview.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            columns = reader.fieldnames or []
            gene_column = next((c for c in ("Gene ID", "Gene_ID", "gene_id") if c in columns), None)
            evidence_columns = [c for c in columns if c.lower().replace("-", "_") in
                                ("hmmer", "dbcan_hmm", "dbcan_sub", "dbcan_sub_hmm")]
            recommended_column = "Recommend Results" if "Recommend Results" in columns else None
            tool_count_column = "#ofTools" if "#ofTools" in columns else None
            if not gene_column or not (evidence_columns or recommended_column):
                raise ValueError(f"Unexpected dbCAN overview columns in {overview}: {columns}")
            family_genes: dict[str, set[str]] = defaultdict(set)
            for hit in reader:
                gene = hit[gene_column].strip()
                if not gene.startswith(mag + "_contig_"):
                    raise ValueError(f"Foreign ORF ID in {overview}: {gene}")
                if tool_count_column and int(hit[tool_count_column]) < 2:
                    continue
                annotation = hit[recommended_column] if recommended_column else ""
                if not FAMILY.search(annotation):
                    annotation = " ".join(hit[column] for column in evidence_columns)
                for family in FAMILY.findall(annotation):
                    family_genes[family.split("_", 1)[0]].add(gene)
        mag_families[mag] = family_genes
        all_families.update(family_genes)
    families = sorted(all_families)
    mag_rows = []
    sgb_rows = []
    for row in manifest:
        for family in families:
            count = len(mag_families[row["mag_id"]].get(family, ()))
            mag_rows.append({"mag_id": row["mag_id"], "sgb_id": row["sgb_id"],
                             "family": family, "orf_count": count, "present": int(count > 0)})
    prevalence: dict[tuple[str, str], float] = {}
    for sgb in sgb_order:
        members = by_sgb[sgb]
        rep = next(r["mag_id"] for r in manifest if r["sgb_id"] == sgb and r["is_representative"] == "true")
        for family in families:
            positive = sum(bool(mag_families[m].get(family)) for m in members)
            value = positive / len(members)
            prevalence[sgb, family] = value
            sgb_rows.append({"sgb_id": sgb, "family": family, "members_with_family": positive,
                             "members_total": len(members), "prevalence": f"{value:.8g}",
                             "representative_orf_count": len(mag_families[rep].get(family, ()))})
    out = args.out / "10_trait_matrix"
    write_table(out / "mag_x_cazy_family_long.tsv", ["mag_id", "sgb_id", "family", "orf_count", "present"], mag_rows)
    write_table(out / "sgb_x_cazy_family_long.tsv", ["sgb_id", "family", "members_with_family", "members_total",
                                                 "prevalence", "representative_orf_count"], sgb_rows)
    abundance = args.pipeline_root / "tables/sample_x_sgb_tpm.tsv"
    with abundance.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        check_abundance_columns(reader.fieldnames, sgb_order, args.out)
        sample_rows = []
        seen_samples = set()
        for row in reader:
            sample = row["sample_id"]
            if not sample or sample in seen_samples:
                raise ValueError(f"Empty or duplicate sample ID: {sample}")
            seen_samples.add(sample)
            values = {sgb: float(row[sgb]) for sgb in sgb_order}
            if any(not math.isfinite(v) or v < 0 for v in values.values()):
                raise ValueError(f"Invalid TPM in sample {sample}")
            sample_rows.append({"sample_id": sample, **{
                family: f"{sum(values[sgb] * prevalence[sgb, family] for sgb in sgb_order):.8g}"
                for family in families}})
    if not sample_rows:
        raise ValueError(f"Empty Salmon abundance table: {abundance}")
    write_table(args.out / "11_abundance_weighted/sample_x_cazy_family_weighted_tpm.tsv",
                ["sample_id", *families], sample_rows)
    print(f"Families={len(families)} MAGs={len(manifest)} SGBs={len(sgb_order)} samples={len(sample_rows)}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pipeline-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    sub = parser.add_subparsers(dest="command", required=True)
    prepare_parser = sub.add_parser("prepare")
    prepare_parser.add_argument("--sgb-list", type=Path)
    gene_parser = sub.add_parser("genes")
    gene_parser.add_argument("--mag", required=True)
    sub.add_parser("summarize")
    args = parser.parse_args()
    {"prepare": prepare, "genes": genes, "summarize": summarize}[args.command](args)


if __name__ == "__main__":
    main()
