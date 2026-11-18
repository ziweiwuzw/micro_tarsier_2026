#!/usr/bin/env python3
"""Prepare a joint SemiBin2/VAMB per-group BAM workflow.

This script deliberately reuses an existing SemiBin2-by-group preparation:
``sample_inputs.tsv``, ``prepare.done`` and each group's
``concatenated.fa.gz``. It adds canonical read profiles either from an optional
VAMB manifest or directly by sampling each R1, then writes job manifests for
Strobealign BAM generation.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
from pathlib import Path


def require_file(path: Path, label: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit(f"Missing/empty {label}: {path}")


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]", "_", value.strip()) or "unnamed"


def safe_sample(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]", "_", value.strip()) or "unnamed"


def fasta_header_check(path: Path) -> tuple[int, int]:
    """Return (header_count, headers_without_colon) without loading sequences."""
    headers = 0
    missing_colon = 0
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:  # type: ignore[arg-type]
        for raw in handle:
            if raw.startswith(">"):
                headers += 1
                text = raw[1:].strip()
                identifier = text.split()[0] if text else ""
                if ":" not in identifier:
                    missing_colon += 1
    return headers, missing_colon


def read_profiles(path: Path) -> dict[str, dict[str, str]]:
    require_file(path, "VAMB manifest")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_id", "r1", "r2", "canonical_profile"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"VAMB manifest missing columns: {', '.join(sorted(missing))}")
        result: dict[str, dict[str, str]] = {}
        for line_no, raw in enumerate(reader, start=2):
            row = {key: (value or "").strip() for key, value in raw.items()}
            sample = row["sample_id"]
            profile = row["canonical_profile"]
            if not sample or sample in result:
                raise SystemExit(f"VAMB manifest line {line_no}: empty/duplicate sample_id: {sample!r}")
            if not profile.isdigit() or int(profile) < 1:
                raise SystemExit(f"VAMB manifest line {line_no}: invalid canonical_profile: {profile!r}")
            for key in ("r1", "r2"):
                require_file(Path(row[key]).expanduser(), f"{sample} {key}")
            result[sample] = row
    if not result:
        raise SystemExit("VAMB manifest contains no samples")
    return result


def read_semibin_inputs(root: Path, profiles: dict[str, dict[str, str]] | None) -> list[dict[str, str]]:
    require_file(root / "prepare.done", "SemiBin2 prepare.done")
    sample_inputs = root / "sample_inputs.tsv"
    require_file(sample_inputs, "SemiBin2 sample_inputs.tsv")
    with sample_inputs.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_id", "canonical_species", "contigs", "r1", "r2"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"SemiBin2 sample_inputs.tsv missing columns: {', '.join(sorted(missing))}")
        rows: list[dict[str, str]] = []
        seen: set[str] = set()
        for line_no, raw in enumerate(reader, start=2):
            row = {key: (value or "").strip() for key, value in raw.items()}
            sample = row["sample_id"]
            group = row["canonical_species"]
            if not sample or sample in seen:
                raise SystemExit(f"SemiBin2 sample_inputs.tsv line {line_no}: empty/duplicate sample_id: {sample!r}")
            if not group:
                raise SystemExit(f"SemiBin2 sample_inputs.tsv line {line_no}: empty canonical_species")
            if profiles is not None and sample not in profiles:
                raise SystemExit(f"Sample in SemiBin2 inputs but not VAMB manifest: {sample}")
            for key in ("contigs", "r1", "r2"):
                require_file(Path(row[key]).expanduser(), f"{sample} {key}")
            if profiles is not None:
                for key in ("r1", "r2"):
                    left = Path(row[key]).expanduser().resolve()
                    right = Path(profiles[sample][key]).expanduser().resolve()
                    if left != right:
                        raise SystemExit(f"Read path mismatch for {sample} {key}: {left} != {right}")
            seen.add(sample)
            row["canonical_group"] = group
            rows.append(row)
    missing = sorted(set(profiles) - seen) if profiles is not None else []
    if missing:
        raise SystemExit(f"VAMB manifest has {len(missing)} samples absent from SemiBin2 inputs; first: {missing[:10]}")
    if not rows:
        raise SystemExit("SemiBin2 sample_inputs.tsv contains no samples")
    return rows


CANONICAL_LENGTHS = (50, 75, 100, 125, 150, 250, 400)


def derive_profile(r1: Path, limit: int) -> tuple[int, float]:
    lengths: list[int] = []
    with gzip.open(r1, "rt", encoding="ascii", errors="strict") as handle:
        for _ in range(limit):
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not sequence or not plus or not quality or not header.startswith("@") or not plus.startswith("+"):
                raise SystemExit(f"Malformed/truncated FASTQ while estimating profile: {r1}")
            sequence = sequence.rstrip("\r\n")
            quality = quality.rstrip("\r\n")
            if not sequence or len(sequence) != len(quality):
                raise SystemExit(f"Invalid FASTQ record while estimating profile: {r1}")
            lengths.append(len(sequence))
    if not lengths:
        raise SystemExit(f"FASTQ contains no reads: {r1}")
    mean_length = sum(lengths) / len(lengths)
    profile = min(CANONICAL_LENGTHS, key=lambda value: (abs(value - mean_length), value))
    return profile, mean_length


def add_derived_profiles(rows: list[dict[str, str]], limit: int) -> None:
    if limit < 1:
        raise SystemExit("--profile-reads must be positive")
    for row in rows:
        profile, mean_length = derive_profile(Path(row["r1"]), limit)
        row["canonical_profile"] = str(profile)
        row["mean_read_length"] = f"{mean_length:.3f}"


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_suffix(path.suffix + ".part")
    with part.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    part.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--semibin-root", required=True, type=Path)
    parser.add_argument("--vamb-manifest", type=Path, help="optional; if omitted, estimate profiles from each R1")
    parser.add_argument("--profile-reads", type=int, default=500, help="R1 reads sampled when --vamb-manifest is omitted")
    parser.add_argument("--bam-dir-name", default="joint_bams")
    parser.add_argument("--semibin-out-name", default="joint_semibin2")
    parser.add_argument("--vamb-out-name", default="joint_vamb")
    parser.add_argument("--force", action="store_true", help="overwrite generated metadata/job manifests")
    args = parser.parse_args()
    root = args.semibin_root.expanduser().resolve()
    require_file(root / "prepare.done", "SemiBin2 prepare.done")
    profiles = read_profiles(args.vamb_manifest.expanduser().resolve()) if args.vamb_manifest else None
    rows = read_semibin_inputs(root, profiles)
    if profiles is None:
        add_derived_profiles(rows, args.profile_reads)
    else:
        for row in rows:
            row["canonical_profile"] = profiles[row["sample_id"]]["canonical_profile"]

    groups: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        groups.setdefault(row["canonical_group"], []).append(row)

    groups_out: list[dict[str, str]] = []
    for index, group in enumerate(sorted(groups), start=1):
        slug = safe_name(group)
        group_dir = root / slug
        reference = group_dir / "concatenated.fa.gz"
        require_file(reference, f"SemiBin2 group reference for {group}")
        header_count, missing_colon = fasta_header_check(reference)
        if header_count == 0:
            raise SystemExit(f"Reference has no FASTA headers: {reference}")
        if missing_colon:
            raise SystemExit(f"Reference has {missing_colon}/{header_count} headers without ':'; cannot use VAMB -o ':' safely: {reference}")

        group_rows = groups[group]
        bam_dir = group_dir / args.bam_dir_name
        semibin_out = group_dir / args.semibin_out_name
        vamb_out = group_dir / args.vamb_out_name
        jobs_path = group_dir / "joint_mapping_jobs.tsv"
        sample_names: set[str] = set()
        job_rows: list[dict[str, str]] = []
        for task, row in enumerate(group_rows, start=1):
            sample = row["sample_id"]
            bam_stem = safe_sample(sample)
            if bam_stem in sample_names:
                raise SystemExit(f"Sample filename collision in group {group}: {sample}")
            sample_names.add(bam_stem)
            job_rows.append({
                "array_index": str(task),
                "sample_id": sample,
                "group": group,
                "group_slug": slug,
                "group_dir": str(group_dir),
                "reference": str(reference),
                "r1": row["r1"],
                "r2": row["r2"],
                "canonical_profile": row["canonical_profile"],
                "bam": str(bam_dir / f"{bam_stem}.sorted.bam"),
            })
        if jobs_path.exists() and not args.force:
            raise SystemExit(f"Generated job file exists; use --force: {jobs_path}")
        write_tsv(jobs_path, list(job_rows[0]), job_rows)
        groups_out.append({
            "array_index": str(index),
            "group": group,
            "group_slug": slug,
            "group_dir": str(group_dir),
            "reference": str(reference),
            "mapping_jobs": str(jobs_path),
            "bam_dir": str(bam_dir),
            "semibin_out": str(semibin_out),
            "vamb_out": str(vamb_out),
            "sample_count": str(len(group_rows)),
            "reference_contigs": str(header_count),
        })

    groups_path = root / "joint_groups.tsv"
    if groups_path.exists() and not args.force:
        raise SystemExit(f"Generated groups manifest exists; use --force: {groups_path}")
    fields = ["array_index", "group", "group_slug", "group_dir", "reference", "mapping_jobs", "bam_dir", "semibin_out", "vamb_out", "sample_count", "reference_contigs"]
    write_tsv(groups_path, fields, groups_out)
    print(f"[PASS] groups={len(groups_out)} samples={len(rows)}")
    print(f"[PASS] groups manifest: {groups_path}")
    source = "VAMB manifest" if profiles is not None else f"R1 first {args.profile_reads} reads"
    print(f"[PASS] canonical profiles from {source}")
    print("[PASS] references use ':'-separated FASTA identifiers for VAMB binsplitting")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
