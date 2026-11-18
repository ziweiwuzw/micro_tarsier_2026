#!/usr/bin/env python3
"""Stage canonicalized SemiBin2, VAMB, and COMEBin FASTAs for MetaWRAP."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import re
import shutil
import sys
from pathlib import Path
from typing import Iterator, Optional


FASTA_RE = re.compile(r"\.(?:fa|fna|fasta)(?:\.gz)?$", re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--joint-root", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--pilot-count", type=int, default=3)
    parser.add_argument("--groups", default="", help="Comma-separated group/group_slug list")
    parser.add_argument(
        "--all-groups",
        action="store_true",
        help="Stage every eligible row in joint_groups.tsv and fail if any row is ineligible",
    )
    parser.add_argument(
        "--exclude-comebin-groups",
        default="",
        help="Comma-separated group/group_slug list for which COMEBin must be ignored",
    )
    parser.add_argument("--allow-unmarked", action="store_true")
    return parser.parse_args()


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip()).strip("_") or "unnamed"


def choose_bins_dir(output: Path) -> Optional[Path]:
    for candidate in (output / "bins", output / "output_bins"):
        if candidate.is_dir() and any(
            p.is_file() and FASTA_RE.search(p.name) for p in candidate.iterdir()
        ):
            return candidate
    return None


def fasta_files(directory: Path) -> list[Path]:
    return sorted(
        p for p in directory.iterdir() if p.is_file() and FASTA_RE.search(p.name)
    )


def open_text(path: Path):
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="strict")
    return path.open("rt", encoding="utf-8", errors="strict")


def fasta_records(path: Path) -> Iterator[tuple[str, str]]:
    identifier: Optional[str] = None
    sequence: list[str] = []
    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    yield identifier, "".join(sequence).upper()
                identifier = line[1:].split()[0] if line[1:].strip() else ""
                if not identifier:
                    raise RuntimeError(f"Empty FASTA header: {path}:{line_number}")
                sequence = []
            else:
                if identifier is None:
                    raise RuntimeError(f"Sequence before first header: {path}:{line_number}")
                sequence.append(re.sub(r"\s+", "", line).upper())
    if identifier is not None:
        yield identifier, "".join(sequence).upper()


def digest_sequence(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def validate_against_reference(reference: Path, staged_dirs: list[Path]) -> None:
    """Ensure every staged canonical contig is present and sequence-identical."""
    expected: dict[str, str] = {}
    for directory in staged_dirs:
        for path in fasta_files(directory):
            for identifier, sequence in fasta_records(path):
                digest = digest_sequence(sequence)
                previous = expected.get(identifier)
                if previous is not None and previous != digest:
                    raise RuntimeError(
                        f"Canonical contig has different sequences across binners: {identifier}"
                    )
                expected[identifier] = digest

    found: set[str] = set()
    seen_reference: set[str] = set()
    for identifier, sequence in fasta_records(reference):
        if identifier in seen_reference:
            raise RuntimeError(f"Duplicate contig ID in reference assembly: {identifier}")
        seen_reference.add(identifier)
        if identifier not in expected:
            continue
        found.add(identifier)
        if digest_sequence(sequence) != expected[identifier]:
            raise RuntimeError(
                f"Canonical contig sequence differs from reference assembly: {identifier}"
            )

    missing = sorted(set(expected) - found)
    if missing:
        preview = ", ".join(missing[:5])
        raise RuntimeError(
            f"{len(missing)} staged canonical contigs are absent from reference {reference}; "
            f"first: {preview}"
        )


def fasta_stem(path: Path) -> str:
    name = path.name
    lower = name.lower()
    for suffix in (".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa"):
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    return name


def sample_prefix(path: Path, method: str) -> str:
    stem = fasta_stem(path)
    if method == "semibin2":
        marker = "_SemiBin_"
        if marker not in stem:
            raise RuntimeError(f"Cannot infer SemiBin2 sample prefix from filename: {path.name}")
        return stem.split(marker, 1)[0]
    if method == "comebin":
        if ":" not in stem:
            raise RuntimeError(f"Cannot infer COMEBin sample prefix from filename: {path.name}")
        return stem.rsplit(":", 1)[0]
    return ""


def canonical_contig_id(identifier: str, path: Path, method: str) -> str:
    if method == "vamb":
        return identifier
    prefix = sample_prefix(path, method)
    marker = prefix + ":"
    return identifier if identifier.startswith(marker) else marker + identifier


def normalized_fasta_lines(path: Path, method: str) -> Iterator[str]:
    seen_header = False
    sequence_bases = 0
    seen_ids: set[str] = set()
    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                identifier = line[1:].strip().split()[0] if line[1:].strip() else ""
                if not identifier:
                    raise RuntimeError(f"Empty FASTA header: {path}:{line_number}")
                identifier = canonical_contig_id(identifier, path, method)
                if identifier in seen_ids:
                    raise RuntimeError(f"Duplicate contig assignment in {path}: {identifier}")
                seen_ids.add(identifier)
                seen_header = True
                yield f">{identifier}\n"
            else:
                if not seen_header:
                    raise RuntimeError(f"Sequence before first header: {path}:{line_number}")
                sequence = re.sub(r"\s+", "", line).upper()
                if not sequence:
                    continue
                if re.search(r"[^ACGTURYSWKMBDHVN.-]", sequence):
                    raise RuntimeError(f"Invalid FASTA sequence characters: {path}:{line_number}")
                sequence_bases += len(sequence)
                yield sequence + "\n"
    if not seen_header or sequence_bases == 0:
        raise RuntimeError(f"Empty FASTA content: {path}")


def stage_set(source_dir: Path, target_dir: Path, prefix: str, method: str) -> list[tuple[str, str]]:
    files = fasta_files(source_dir)
    if not files:
        raise RuntimeError(f"No FASTA bins in {source_dir}")
    target_dir.mkdir(parents=True, exist_ok=False)
    provenance: list[tuple[str, str]] = []
    for index, source in enumerate(files, start=1):
        staged_name = f"{prefix}_{index:06d}.fa"
        target = target_dir / staged_name
        with target.open("wt", encoding="utf-8", newline="\n") as handle:
            handle.writelines(normalized_fasta_lines(source, method))
        provenance.append((staged_name, str(source)))
    return provenance


def load_groups(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit(f"Missing joint_groups.tsv: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "array_index", "group", "group_slug", "reference", "semibin_out", "vamb_out"
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"joint_groups.tsv missing columns: {', '.join(sorted(missing))}")
        return [{k: (v or "").strip() for k, v in row.items()} for row in reader]


def main() -> None:
    args = parse_args()
    if args.pilot_count < 1:
        raise SystemExit("--pilot-count must be positive")
    if args.all_groups and args.groups.strip():
        raise SystemExit("--all-groups and --groups are mutually exclusive")
    manifest = args.output_root / "pilot_groups.tsv"

    requested = [x.strip() for x in args.groups.split(",") if x.strip()]
    requested_set = set(requested)
    exclude_comebin = {
        x.strip() for x in args.exclude_comebin_groups.split(",") if x.strip()
    }
    eligible: list[tuple[dict[str, str], Optional[Path], Optional[Path], Optional[Path]]] = []
    rejected: list[str] = []
    for row in load_groups(args.joint_root / "joint_groups.tsv"):
        if requested_set and row["group"] not in requested_set and row["group_slug"] not in requested_set:
            continue
        semibin_out = Path(row["semibin_out"])
        vamb_out = Path(row["vamb_out"])
        comebin_out = (
            Path(row["comebin_out"])
            if row.get("comebin_out")
            else semibin_out.parent / "comebin"
        )
        semibin_bins = choose_bins_dir(semibin_out)
        vamb_bins = choose_bins_dir(vamb_out)
        comebin_is_excluded = (
            row["group"] in exclude_comebin or row["group_slug"] in exclude_comebin
        )
        comebin_bins = None if comebin_is_excluded else choose_bins_dir(comebin_out)
        if comebin_is_excluded:
            print(f"[EXCLUDE] {row['group_slug']}: COMEBin disabled by request")
        marked = (
            (semibin_bins is None or (semibin_out / ".joint_semibin2.done").is_file())
            and (vamb_bins is None or (vamb_out / ".joint_vamb.done").is_file())
        )
        reference = Path(row["reference"])
        available_methods = sum(x is not None for x in (semibin_bins, vamb_bins, comebin_bins))
        if not reference.is_file() or reference.stat().st_size == 0:
            rejected.append(f"{row['group_slug']}: missing reference")
        elif available_methods == 0:
            rejected.append(
                f"{row['group_slug']}: no non-empty SemiBin2/VAMB/COMEBin bins directory"
            )
        elif not marked and not args.allow_unmarked:
            rejected.append(f"{row['group_slug']}: missing workflow completion marker")
        else:
            row = dict(row)
            row["comebin_out"] = str(comebin_out)
            eligible.append((row, semibin_bins, vamb_bins, comebin_bins))

    if args.all_groups:
        if rejected:
            raise SystemExit(
                "ALL_GROUPS requested, but one or more groups are ineligible:\n"
                + "\n".join(rejected)
            )
        selected = eligible
    elif requested_set:
        found = {row["group"] for row, _, _, _ in eligible} | {row["group_slug"] for row, _, _, _ in eligible}
        missing_requested = [name for name in requested if name not in found]
        if missing_requested:
            raise SystemExit(
                "Requested groups are not eligible: " + ", ".join(missing_requested)
                + "\n" + "\n".join(rejected)
            )
        selected = eligible
    else:
        selected = eligible[: args.pilot_count]
    expected_count = len(selected) if args.all_groups else (
        len(requested) if requested else args.pilot_count
    )
    if not selected:
        raise SystemExit("No eligible groups were selected")
    if len(selected) < expected_count:
        raise SystemExit(
            f"Not enough eligible groups: selected={len(selected)} requested="
            f"{expected_count}\n" + "\n".join(rejected)
        )

    expected_slugs = [safe_name(row["group_slug"]) for row, _, _, _ in selected]
    if manifest.is_file() and manifest.stat().st_size > 0:
        with manifest.open(newline="", encoding="utf-8") as handle:
            existing_rows = list(csv.DictReader(handle, delimiter="\t"))
        existing_slugs = [row.get("group_slug", "") for row in existing_rows]
        if existing_slugs != expected_slugs:
            raise SystemExit(
                f"Existing manifest selection differs from this request: {manifest}. "
                "Use a new --output-root for the full run."
            )
        for row in existing_rows:
            if (
                row.get("group", "") in exclude_comebin
                or row.get("group_slug", "") in exclude_comebin
            ) and row.get("comebin_input", ""):
                raise SystemExit(
                    f"Existing staged group still contains excluded COMEBin input: "
                    f"{row.get('group_slug', '')}. Use a new --output-root."
                )
        print(f"[REUSE] {manifest}")
        return

    args.output_root.mkdir(parents=True, exist_ok=True)
    groups_root = args.output_root / "groups"
    groups_root.mkdir(exist_ok=True)
    manifest_tmp = manifest.with_name(manifest.name + f".part.{os.getpid()}")
    rows: list[dict[str, str | int]] = []

    for pilot_index, (row, semibin_bins, vamb_bins, comebin_bins) in enumerate(selected, start=1):
        slug = safe_name(row["group_slug"])
        final_group = groups_root / slug
        prepare_done = final_group / ".prepare.done"
        if final_group.exists():
            if prepare_done.is_file():
                with (final_group / "input_provenance.tsv").open(newline="", encoding="utf-8") as handle:
                    provenance_rows = list(csv.DictReader(handle, delimiter="\t"))
                semibin_count = sum(x["method"] == "semibin2" for x in provenance_rows)
                vamb_count = sum(x["method"] == "vamb" for x in provenance_rows)
                comebin_count = sum(x["method"] == "comebin" for x in provenance_rows)
                if (
                    row["group"] in exclude_comebin
                    or row["group_slug"] in exclude_comebin
                ) and comebin_count:
                    raise SystemExit(
                        f"Existing staged group contains COMEBin but it is now excluded: "
                        f"{final_group}. Use a new --output-root."
                    )
            else:
                raise SystemExit(f"Incomplete staged group already exists: {final_group}")
        else:
            partial = groups_root / f".{slug}.part.{os.getpid()}"
            if partial.exists():
                shutil.rmtree(partial)
            partial.mkdir()
            try:
                semibin_prov = (
                    stage_set(semibin_bins, partial / "inputs" / "semibin2", "SB2", "semibin2")
                    if semibin_bins is not None else []
                )
                vamb_prov = (
                    stage_set(vamb_bins, partial / "inputs" / "vamb", "VAMB", "vamb")
                    if vamb_bins is not None else []
                )
                comebin_prov = (
                    stage_set(comebin_bins, partial / "inputs" / "comebin", "COMEBIN", "comebin")
                    if comebin_bins is not None else []
                )
                validate_against_reference(
                    Path(row["reference"]),
                    [directory for directory in (
                        partial / "inputs" / "semibin2",
                        partial / "inputs" / "vamb",
                        partial / "inputs" / "comebin",
                    ) if directory.is_dir()],
                )
                with (partial / "input_provenance.tsv").open("w", newline="", encoding="utf-8") as handle:
                    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                    writer.writerow(["method", "staged_filename", "source_path"])
                    for name, source in semibin_prov:
                        writer.writerow(["semibin2", name, source])
                    for name, source in vamb_prov:
                        writer.writerow(["vamb", name, source])
                    for name, source in comebin_prov:
                        writer.writerow(["comebin", name, source])
                (partial / ".prepare.done").write_text("status=SUCCESS\n", encoding="utf-8")
                partial.rename(final_group)
                semibin_count = len(semibin_prov)
                vamb_count = len(vamb_prov)
                comebin_count = len(comebin_prov)
            except Exception:
                shutil.rmtree(partial, ignore_errors=True)
                raise

        rows.append({
            "pilot_index": pilot_index,
            "source_array_index": row["array_index"],
            "group": row["group"],
            "group_slug": slug,
            "reference": row["reference"],
            "semibin_input": str(final_group / "inputs" / "semibin2") if semibin_count else "",
            "vamb_input": str(final_group / "inputs" / "vamb") if vamb_count else "",
            "comebin_input": str(final_group / "inputs" / "comebin") if comebin_count else "",
            "group_output": str(final_group),
            "semibin_count": semibin_count,
            "vamb_count": vamb_count,
            "comebin_count": comebin_count,
        })

    fieldnames = list(rows[0])
    with manifest_tmp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    manifest_tmp.replace(manifest)
    print(f"[SUCCESS] staged {len(rows)} groups: {manifest}")
    for row in rows:
        print(
            f"  {row['pilot_index']}\t{row['group_slug']}\t"
            f"SB2={row['semibin_count']}\tVAMB={row['vamb_count']}\t"
            f"COMEBin={row['comebin_count']}"
        )


if __name__ == "__main__":
    main()
