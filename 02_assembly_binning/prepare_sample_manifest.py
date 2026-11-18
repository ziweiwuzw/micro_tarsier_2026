#!/usr/bin/env python3
"""Discover MEGAHIT assemblies and validate the COMEBin per-sample manifest."""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from collections import Counter
from pathlib import Path


DEFAULT_ASSEMBLY_ROOT = Path(
    "/path/to/project/03_megahit_sample_assemblies/assemblies"
)
DEFAULT_OUTPUT = Path(
    "/path/to/project/"
    "06_comebin/per_sample_manifest.tsv"
)
DEFAULT_READ_ROOTS = (
    Path(
        "/path/to/scratch/"
        "02_host_removed_reads/reads"
    ),
    Path(
        "/path/to/legacy_data/"
        "02_host_removed_reads/reads"
    ),
    Path(
        "/path/to/project/"
        "02_host_removed_reads/reads"
    ),
    Path(
        "/path/to/legacy_scratch/"
        "02_host_removed_reads/reads"
    ),
)

ASSEMBLY_SUFFIX = "_final.contigs.fa"
R1_SUFFIX = "_1.nonhost.fastq.gz"
R2_SUFFIX = "_2.nonhost.fastq.gz"
SAFE_COMPONENT_RE = re.compile(r"^[A-Za-z0-9._-]+$")


def fail(message: str) -> None:
    raise RuntimeError(message)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Scan paired nonhost reads across several filesystem roots and join "
            "them to discovered per-sample MEGAHIT assemblies."
        )
    )
    assembly_source = parser.add_mutually_exclusive_group()
    assembly_source.add_argument(
        "--assembly-root", type=Path, default=DEFAULT_ASSEMBLY_ROOT,
        help="recursively discover per-sample MEGAHIT FASTA files under this root",
    )
    assembly_source.add_argument(
        "--assembly-list", type=Path,
        help="use an existing file with one absolute MEGAHIT assembly path per line",
    )
    parser.add_argument(
        "--read-root",
        action="append",
        type=Path,
        dest="read_roots",
        help=(
            "reads root with species/project/sample FASTQs; repeat for multiple "
            "roots (defaults to the four transfer destinations)"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"output TSV (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--expected-samples",
        type=int,
        default=547,
        help="required final sample count; use 0 to disable (default: 547)",
    )
    return parser.parse_args()


def validate_component(label: str, value: str, source: Path) -> None:
    if not SAFE_COMPONENT_RE.fullmatch(value):
        fail(f"Unsafe {label} {value!r} parsed from: {source}")


def validate_nonempty(label: str, path: Path) -> int:
    if not path.is_file():
        fail(f"{label} does not exist: {path}")
    size = path.stat().st_size
    if size == 0:
        fail(f"{label} is empty: {path}")
    return size


def validate_gzip_header(label: str, path: Path) -> None:
    with path.open("rb") as handle:
        magic = handle.read(2)
    if magic != b"\x1f\x8b":
        fail(f"{label} does not begin with the gzip magic bytes: {path}")


def count_fasta_records(path: Path) -> int:
    count = 0
    with path.open("rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                count += 1
    return count


def assembly_paths(root: Path, path: Path | None) -> tuple[list[Path], str]:
    if path is not None:
        validate_nonempty("Assembly list", path)
        with path.open("rt", encoding="utf-8") as handle:
            paths = [Path(line.strip()) for line in handle if line.strip()]
        source = str(path)
    else:
        if not root.is_dir():
            fail(f"Assembly root does not exist or is not a directory: {root}")
        paths = sorted(root.rglob(f"*{ASSEMBLY_SUFFIX}"))
        source = str(root)
    if not paths:
        fail(f"No MEGAHIT assemblies found in: {source}")
    return paths, source


def load_assemblies(paths: list[Path], source: str) -> dict[str, dict[str, object]]:
    assemblies: dict[str, dict[str, object]] = {}
    for line_number, assembly in enumerate(paths, start=1):
        if not assembly.is_absolute():
            fail(f"Assembly path is not absolute at {source}:{line_number}: {assembly}")
        if not assembly.name.endswith(ASSEMBLY_SUFFIX):
            fail(f"Unexpected assembly filename at {source}:{line_number}: {assembly}")

        sample_id = assembly.name[: -len(ASSEMBLY_SUFFIX)]
        project = assembly.parent.name
        species = assembly.parent.parent.name
        validate_component("sample_id", sample_id, assembly)
        validate_component("project", project, assembly)
        validate_component("species", species, assembly)

        if sample_id in assemblies:
            fail(
                "Duplicate assembly sample ID "
                f"{sample_id}:\n  {assemblies[sample_id]['path']}\n  {assembly}"
            )

        assembly_bytes = validate_nonempty("Assembly", assembly)
        contig_count = count_fasta_records(assembly)
        if contig_count == 0:
            fail(f"Assembly contains no FASTA records: {assembly}")

        assemblies[sample_id] = {
            "sample_id": sample_id,
            "species": species,
            "project": project,
            "path": assembly,
            "bytes": assembly_bytes,
            "contigs": contig_count,
        }
    return assemblies


def load_reads(read_roots: list[Path]) -> dict[str, dict[str, object]]:
    reads: dict[str, dict[str, object]] = {}

    for root in read_roots:
        if not root.is_dir():
            fail(f"Read root does not exist or is not a directory: {root}")

        r1_paths = set(root.rglob(f"*{R1_SUFFIX}"))
        r2_paths = set(root.rglob(f"*{R2_SUFFIX}"))
        paired_r2_paths: set[Path] = set()

        for r1 in sorted(r1_paths):
            relative = r1.relative_to(root)
            if len(relative.parts) != 3:
                fail(
                    "Expected read path species/project/filename below root, got: "
                    f"{relative} (root: {root})"
                )

            species, project, filename = relative.parts
            sample_id = filename[: -len(R1_SUFFIX)]
            validate_component("species", species, r1)
            validate_component("project", project, r1)
            validate_component("sample_id", sample_id, r1)

            r2 = r1.with_name(f"{sample_id}{R2_SUFFIX}")
            if r2 not in r2_paths:
                fail(f"R2 is missing for {sample_id}: expected {r2}")
            paired_r2_paths.add(r2)

            if sample_id in reads:
                fail(
                    "Sample has reads in more than one location "
                    f"({sample_id}):\n  {reads[sample_id]['r1']}\n  {r1}"
                )

            r1_bytes = validate_nonempty("R1", r1)
            r2_bytes = validate_nonempty("R2", r2)
            validate_gzip_header("R1", r1)
            validate_gzip_header("R2", r2)

            reads[sample_id] = {
                "sample_id": sample_id,
                "species": species,
                "project": project,
                "r1": r1,
                "r2": r2,
                "root": root,
                "read_bytes": r1_bytes + r2_bytes,
            }

        orphan_r2 = sorted(r2_paths - paired_r2_paths)
        if orphan_r2:
            preview = "\n".join(f"  {item}" for item in orphan_r2[:20])
            fail(f"Found {len(orphan_r2)} R2 files without matching R1:\n{preview}")

    if not reads:
        fail("No paired nonhost reads were found in the configured roots")
    return reads


def format_examples(values: list[str], limit: int = 20) -> str:
    return ", ".join(values[:limit])


def write_manifest(
    output: Path,
    assemblies: dict[str, dict[str, object]],
    reads: dict[str, dict[str, object]],
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + f".part.{os.getpid()}")
    fieldnames = [
        "array_index",
        "sample_id",
        "species",
        "project",
        "contigs",
        "r1",
        "r2",
        "assembly_contig_count",
        "assembly_bytes",
    ]

    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=fieldnames,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()

            for array_index, sample_id in enumerate(sorted(reads), start=1):
                assembly = assemblies[sample_id]
                read = reads[sample_id]
                writer.writerow(
                    {
                        "array_index": array_index,
                        "sample_id": sample_id,
                        "species": read["species"],
                        "project": read["project"],
                        "contigs": assembly["path"],
                        "r1": read["r1"],
                        "r2": read["r2"],
                        "assembly_contig_count": assembly["contigs"],
                        "assembly_bytes": assembly["bytes"],
                    }
                )

        temporary.replace(output)
    finally:
        temporary.unlink(missing_ok=True)


def main() -> None:
    args = parse_args()
    if args.expected_samples < 0:
        fail("--expected-samples must be 0 or a positive integer")

    read_roots = args.read_roots or list(DEFAULT_READ_ROOTS)
    if len(set(read_roots)) != len(read_roots):
        fail("The same --read-root was supplied more than once")

    print("Loading and checking MEGAHIT assemblies...", flush=True)
    paths, source = assembly_paths(args.assembly_root, args.assembly_list)
    assemblies = load_assemblies(paths, source)
    print("Scanning and checking paired reads...", flush=True)
    reads = load_reads(read_roots)

    assembly_only = sorted(set(assemblies) - set(reads))
    reads_only = sorted(set(reads) - set(assemblies))
    if assembly_only or reads_only:
        messages = [
            "Read and assembly sample IDs do not match exactly.",
            f"assembly only: {len(assembly_only)}",
            f"reads only: {len(reads_only)}",
        ]
        if assembly_only:
            messages.append("assembly-only examples: " + format_examples(assembly_only))
        if reads_only:
            messages.append("reads-only examples: " + format_examples(reads_only))
        fail("\n".join(messages))

    mismatches = []
    for sample_id in sorted(reads):
        read = reads[sample_id]
        assembly = assemblies[sample_id]
        if (
            read["species"] != assembly["species"]
            or read["project"] != assembly["project"]
        ):
            mismatches.append(
                f"{sample_id}: reads={read['species']}/{read['project']}; "
                f"assembly={assembly['species']}/{assembly['project']}"
            )
    if mismatches:
        fail(
            f"Species/project mismatch for {len(mismatches)} samples:\n  "
            + "\n  ".join(mismatches[:20])
        )

    if args.expected_samples and len(reads) != args.expected_samples:
        fail(
            f"Expected {args.expected_samples} matched samples, found {len(reads)}. "
            "No manifest was published."
        )

    write_manifest(args.output, assemblies, reads)

    root_counts = Counter(str(record["root"]) for record in reads.values())
    root_bytes = Counter()
    for record in reads.values():
        root_bytes[str(record["root"])] += int(record["read_bytes"])

    print("=" * 80)
    print("COMEBin multi-filesystem manifest created")
    print(f"Assemblies: {len(assemblies)}")
    print(f"Paired-read samples: {len(reads)}")
    print(f"Manifest rows: {len(reads)}")
    print(f"Output: {args.output}")
    print("Read distribution:")
    for root in read_roots:
        key = str(root)
        print(
            f"  {key}: {root_counts[key]} samples, "
            f"{root_bytes[key] / (1024 ** 3):.2f} GiB"
        )
    print("[OK] sample IDs, species and projects match exactly")
    print("[OK] all assemblies and paired reads exist and are nonempty")
    print("[OK] every read file begins with a gzip header")
    print("=" * 80)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"[FATAL] {error}", file=sys.stderr)
        sys.exit(1)
