#!/usr/bin/env python3
"""KMA + MicrobeCensus runner for an existing SGB contig catalogue."""

import argparse
import csv
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path


DEFAULT_PIPELINE = Path(
    "/path/to/scratch/08_mag_pipeline/all_species_3binners_v2_kma"
)
PIPELINE = Path(os.environ.get("PIPELINE_ROOT", DEFAULT_PIPELINE))
ROOT = Path(os.environ.get("KMA_ROOT", PIPELINE / "kma"))
READS = Path(os.environ.get("READS_MANIFEST", PIPELINE / "reads.tsv"))
REFERENCE_ROOT = Path(os.environ.get("REFERENCE_ROOT", PIPELINE / "sgb_reference"))
CATALOGUE = REFERENCE_ROOT / "catalogue.fa"
MAP = REFERENCE_ROOT / "contig_to_sgb.tsv"
INDEX = ROOT / "index/sgb_contigs"


def rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def require(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit(f"Missing or empty: {path}")


def reference_token():
    marker = REFERENCE_ROOT / ".reference.done"
    require(marker)
    fields = dict(line.split("=", 1) for line in marker.read_text().splitlines() if "=" in line)
    token = fields.get("map_sha256")
    if not token or len(token) != 64:
        raise ValueError(f"Invalid reference marker: {marker}")
    return token


def sample_token(record):
    payload = json.dumps(record, sort_keys=True, separators=(",", ":")).encode()
    return f"reference_map_sha256={reference_token()}\nreads_row_sha256={hashlib.sha256(payload).hexdigest()}\n"


def catalog_ids():
    require(MAP)
    records = rows(MAP)
    if not records or not {"contig_id", "sgb_id", "length"} <= records[0].keys():
        raise ValueError(f"Unexpected catalogue mapping header: {MAP}")
    names = [r["contig_id"] for r in records]
    if len(names) != len(set(names)):
        raise ValueError("Duplicate contig IDs in catalogue mapping")
    if any(int(r["length"]) <= 0 for r in records):
        raise ValueError("Non-positive reference length")
    return records


def index():
    require(CATALOGUE)
    catalog_ids()
    token = reference_token()
    if (ROOT / "index/.done").exists():
        if f"reference_map_sha256={token}\n" not in (ROOT / "index/.done").read_text():
            raise ValueError("Existing KMA index belongs to another reference")
        print(f"Existing completed index: {INDEX}")
        return
    if INDEX.parent.exists() and any(INDEX.parent.iterdir()):
        raise SystemExit(f"Incomplete index directory exists: {INDEX.parent}")
    INDEX.parent.mkdir(parents=True)
    cmd = ["kma", "index", "-i", str(CATALOGUE), "-o", str(INDEX)]
    print("Building KMA index; no second FASTA copy is made", flush=True)
    subprocess.run(cmd, check=True)
    require(INDEX.with_suffix(".name"))
    (ROOT / "index/.done").write_text(f"reference_map_sha256={token}\ncommand={' '.join(cmd)}\n")


def parse_mapstat(path, known):
    require(path)
    header = None
    result = {}
    input_fragments = None
    with path.open() as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if not line.strip():
                continue
            if line.startswith("##"):
                meta = line[2:].strip().split(None, 1)
                if len(meta) == 2 and meta[0] == "fragmentCount":
                    input_fragments = int(meta[1])
                continue
            if line.startswith("#"):
                candidate = line.lstrip("# ").split("\t")
                if "refSequence" in candidate and "fragmentCount" in candidate:
                    header = candidate
                continue
            fields = line.split("\t")
            if header is None:
                if "refSequence" in fields and "fragmentCount" in fields:
                    header = fields
                    continue
                raise ValueError(f"Could not find mapstat header in {path}")
            if len(fields) != len(header):
                raise ValueError(f"Malformed mapstat row: {line[:100]}")
            record = dict(zip(header, fields))
            name = record["refSequence"]
            if name not in known:
                raise ValueError(f"KMA reference missing from catalogue map: {name}")
            if name in result:
                raise ValueError(f"Duplicate KMA reference: {name}")
            n = int(record["fragmentCount"])
            if n < 0:
                raise ValueError(f"Negative fragment count: {name}")
            if n:
                result[name] = n
    if header is None:
        raise ValueError(f"No mapstat table in {path}")
    # These are per-template counts, not a deduplicated count of input fragments.
    # Their sum can exceed the input fragment total even with -1t1.
    return result, input_fragments


def parse_ge(path):
    require(path)
    match = re.search(
        r"(?m)^\s*genome_equivalents\s*[:\t ]+\s*([0-9.eE+-]+)\s*$",
        path.read_text(),
    )
    if not match:
        raise ValueError(f"Cannot read genome_equivalents from {path}")
    value = float(match.group(1))
    if not 0 < value < float("inf"):
        raise ValueError(f"Invalid genome_equivalents: {value}")
    return value


def sample(sample_index):
    require(ROOT / "index/.done")
    if f"reference_map_sha256={reference_token()}\n" not in (ROOT / "index/.done").read_text():
        raise ValueError("KMA index/reference mismatch")
    known = {r["contig_id"] for r in catalog_ids()}
    require(READS)
    manifest = rows(READS)
    matches = [r for r in manifest if r["sample_index"] == str(sample_index)]
    if len(matches) != 1:
        raise ValueError(f"Expected one reads.tsv row for sample_index={sample_index}")
    record = matches[0]
    sample_id = record["sample_id"]
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", sample_id):
        raise ValueError(f"Unsafe sample_id: {sample_id}")
    r1 = record["read1"].split(";")
    r2 = record["read2"].split(";")
    if len(r1) != len(r2) or not r1 or any(not x for x in r1 + r2):
        raise ValueError(f"Invalid paired lanes for {sample_id}")
    pairs = [Path(x) for lane in zip(r1, r2) for x in lane]
    for path in pairs:
        require(path)

    out = ROOT / "samples" / sample_id
    expected_done = sample_token(record)
    if (out / ".done").exists():
        if (out / ".done").read_text() != expected_done:
            raise ValueError(f"Completed sample has a different reference or reads row: {sample_id}")
        print(f"Already completed: {sample_id}")
        return
    out.mkdir(parents=True, exist_ok=True)

    kma_prefix = out / "kma"
    cmd = [
        "kma", "-ipe", *map(str, pairs), "-t_db", str(INDEX),
        "-o", str(kma_prefix), "-t", os.environ.get("SLURM_CPUS_PER_TASK", "16"),
        "-mem_mode", "-ef", "-1t1", "-nf", "-nc",
    ]
    temp_dir = os.environ.get("KMA_TMPDIR") or os.environ.get("SLURM_TMPDIR")
    if temp_dir:
        if not Path(temp_dir).is_dir():
            raise ValueError(f"KMA temporary directory does not exist: {temp_dir}")
        cmd.extend(["-tmp", temp_dir])
    mapstat = kma_prefix.with_suffix(".mapstat")
    result_file = kma_prefix.with_suffix(".res")
    if mapstat.is_file() and result_file.is_file():
        # A prior job may have completed KMA before MicrobeCensus failed.
        counts, input_fragments = parse_mapstat(mapstat, known)
        require(result_file)
        print(f"Reusing completed KMA mapping for {sample_id}", flush=True)
    else:
        if any(out.iterdir()):
            raise SystemExit(
                f"Incomplete KMA output in {out}; inspect it before rerunning"
            )
        print("Running KMA:", " ".join(cmd), flush=True)
        subprocess.run(cmd, check=True)
        require(result_file)
        counts, input_fragments = parse_mapstat(mapstat, known)

    ge_file = out / "microbecensus.txt"
    ge_cmd = [
        "run_microbe_census.py", "-v", "-t", os.environ.get("SLURM_CPUS_PER_TASK", "16"),
        ",".join(map(str, pairs)), str(ge_file),
    ]
    try:
        ge = parse_ge(ge_file)
        mc_status = "OK"
        print(f"Reusing completed MicrobeCensus report for {sample_id}", flush=True)
    except (SystemExit, ValueError):
        print(f"Running MicrobeCensus for {sample_id}", flush=True)
        ge_env = dict(os.environ)
        if temp_dir:
            ge_env["TMPDIR"] = temp_dir
        run = subprocess.run(
            ge_cmd, env=ge_env, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        )
        if run.stdout:
            print(run.stdout, end="" if run.stdout.endswith("\n") else "\n", flush=True)
        if run.returncode == 0:
            ge = parse_ge(ge_file)
            mc_status = "OK"
        elif "Cannot compute AGS using reads shorter than 50 bp." in run.stdout:
            # KMA counts and within-catalogue TPM remain usable without GE.
            ge = None
            mc_status = "READS_TOO_SHORT_FOR_GE"
            print(f"MicrobeCensus unavailable for {sample_id}: short reads; GE-derived values will be NA", flush=True)
        else:
            raise subprocess.CalledProcessError(run.returncode, ge_cmd, output=run.stdout)
    count_path = out / "contig_counts.tsv"
    with count_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["contig_id", "fragment_count"])
        writer.writerows((name, n) for name, n in sorted(counts.items()))

    subprocess.run(
        [
            "Rscript", str(Path(__file__).with_name("profile_rarefaction.R")),
            str(count_path), str(MAP), str(ge) if ge is not None else "NA", record["coverage_bases"],
            sample_id, str(out),
        ],
        check=True,
    )
    with (out / "mapping_qc.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "sample_id", "kma_input_fragments", "template_fragment_assignments",
            "genome_equivalents", "assignments_per_input_fragment", "microbecensus_status",
        ])
        assigned = sum(counts.values())
        ratio = assigned / input_fragments if input_fragments else ""
        writer.writerow([
            sample_id, input_fragments if input_fragments is not None else "",
            assigned, ge if ge is not None else "", ratio, mc_status,
        ])
    (out / ".done").write_text(expected_done)
    print(f"Completed {sample_id}: assigned={sum(counts.values())} GE={ge} MC={mc_status}", flush=True)


def merge():
    require(READS)
    require(ROOT / "index/.done")
    if f"reference_map_sha256={reference_token()}\n" not in (ROOT / "index/.done").read_text():
        raise ValueError("KMA index/reference mismatch")
    refs = catalog_ids()
    sgbs = list(dict.fromkeys(r["sgb_id"] for r in refs))
    manifest = rows(READS)
    ids = [r["sample_id"] for r in manifest]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate sample IDs in reads manifest")
    profiles = {}
    for sample_id in ids:
        out = ROOT / "samples" / sample_id
        require(out / ".done")
        record = next(row for row in manifest if row["sample_id"] == sample_id)
        if (out / ".done").read_text() != sample_token(record):
            raise ValueError(f"Completed sample has a different reference or reads row: {sample_id}")
        profile = rows(out / "sgb_full.tsv")
        if [r["sgb_id"] for r in profile] != sgbs:
            raise ValueError(f"SGB reference mismatch: {sample_id}")
        profiles[sample_id] = profile
    tables = ROOT / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    for field in ("fragment_count", "fpkb_per_ge", "coverage_per_ge", "tpm", "relative_abundance"):
        with (tables / f"sample_x_sgb_{field}.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample_id", *sgbs])
            for sample_id in ids:
                writer.writerow([sample_id, *(r[field] for r in profiles[sample_id])])
    for filename in ("rarefied_sgb.tsv",):
        with (tables / filename).open("w") as output:
            for i, sample_id in enumerate(ids):
                with (ROOT / "samples" / sample_id / filename).open() as source:
                    header = next(source)
                    if i == 0:
                        output.write(header)
                    for line in source:
                        output.write(line)
    # Older completed samples retain their original per-sample files. Normalize
    # the merged headers and status labels without recomputing KMA or MicrobeCensus.
    with (tables / "rarefaction_status.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "sample_id", "depth", "replicate", "status", "detected_sgbs",
        ], delimiter="\t")
        writer.writeheader()
        for sample_id in ids:
            for row in rows(ROOT / "samples" / sample_id / "rarefaction_status.tsv"):
                if row["sample_id"] != sample_id:
                    raise ValueError(f"Rarefaction sample ID mismatch: {sample_id}")
                if row["status"] == "INSUFFICIENT_MAPPED_FRAGMENTS":
                    row["status"] = "INSUFFICIENT_TEMPLATE_ASSIGNMENTS"
                writer.writerow(row)
    with (tables / "mapping_qc.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "sample_id", "kma_input_fragments", "template_fragment_assignments",
            "genome_equivalents", "assignments_per_input_fragment", "microbecensus_status",
        ], delimiter="\t")
        writer.writeheader()
        for sample_id in ids:
            qc_rows = rows(ROOT / "samples" / sample_id / "mapping_qc.tsv")
            if len(qc_rows) != 1 or qc_rows[0]["sample_id"] != sample_id:
                raise ValueError(f"Expected one mapping QC row for {sample_id}")
            qc = qc_rows[0]
            writer.writerow({
                "sample_id": sample_id,
                "kma_input_fragments": qc["kma_input_fragments"],
                "template_fragment_assignments": qc.get(
                    "template_fragment_assignments", qc.get("catalogue_assigned_fragments")
                ),
                "genome_equivalents": qc["genome_equivalents"],
                "assignments_per_input_fragment": qc.get(
                    "assignments_per_input_fragment", qc.get("catalogue_mapping_fraction")
                ),
                "microbecensus_status": qc.get("microbecensus_status") or (
                    "OK" if qc["genome_equivalents"] else "UNKNOWN"
                ),
            })
    print(f"Merged {len(ids)} samples and {len(sgbs)} SGBs into {tables}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("index", "sample", "merge"))
    parser.add_argument("sample_index", nargs="?", type=int)
    args = parser.parse_args()
    if args.command == "index":
        index()
    elif args.command == "sample":
        if args.sample_index is None:
            parser.error("sample requires sample_index")
        sample(args.sample_index)
    else:
        merge()
