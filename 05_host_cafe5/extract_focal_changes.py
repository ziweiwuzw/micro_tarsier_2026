#!/usr/bin/env python3
"""Extract expansion/contraction results for one focal CAFE5 terminal branch.

The script joins CAFE5's change, reconstructed-count, branch-probability and
family-wide result files to the original CAFE count matrix.  Branch-specific
probabilities and family-wide P values are intentionally kept separate.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


MISSING = {"", "na", "n/a", "n/s", "nan", "none", "-"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create focal-species expansion/contraction TSVs from CAFE5 output."
    )
    parser.add_argument("--results", required=True, type=Path,
                        help="CAFE5 model result directory containing Base_*. files")
    parser.add_argument("--families", required=True, type=Path,
                        help="The exact family matrix used in this CAFE5 run")
    parser.add_argument("--target", default="Cephalopachus_bancanus")
    parser.add_argument("--model", choices=["Base", "Gamma"], default="Base")
    parser.add_argument("--branch-p", type=float, default=0.05,
                        help="Branch-specific probability threshold (default: 0.05)")
    parser.add_argument("--strict-p", type=float, default=0.01,
                        help="Conservative branch threshold reported separately (default: 0.01)")
    parser.add_argument("--outdir", required=True, type=Path)
    return parser.parse_args()


def canonical_node(value: str) -> str:
    value = value.strip().strip("'\"")
    value = re.sub(r"\*?_-?\d+(?:\.\d+)?$", "", value)
    value = value.replace("*", "")
    return value


def species_part(value: str) -> str:
    return re.sub(r"<[^>]+>$", "", canonical_node(value))


def read_tsv(path: Path) -> Tuple[List[str], List[List[str]]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Missing or empty file: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    rows = [row for row in rows if row and any(cell.strip() for cell in row)]
    if not rows:
        raise ValueError(f"No rows found in {path}")
    if len(rows[0]) == 1:
        rows = [re.split(r"\s+", row[0].strip()) for row in rows]
    return [cell.strip() for cell in rows[0]], rows[1:]


def locate_target_column(header: Sequence[str], target: str) -> int:
    exact = [i for i, name in enumerate(header) if canonical_node(name) == target]
    if len(exact) == 1:
        return exact[0]
    species_matches = [
        i for i, name in enumerate(header) if species_part(name) == target
    ]
    if len(species_matches) == 1:
        return species_matches[0]
    raise ValueError(
        f"Could not uniquely locate target {target!r} in header: {list(header)!r}"
    )


def as_float(value: str) -> Optional[float]:
    clean = value.strip().lower()
    if clean in MISSING:
        return None
    try:
        number = float(clean)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def read_target_vector(path: Path, target: str) -> Dict[str, Optional[float]]:
    header, rows = read_tsv(path)
    column = locate_target_column(header, target)
    output: Dict[str, Optional[float]] = {}
    for row in rows:
        if len(row) <= column:
            continue
        family_id = row[0].strip()
        if family_id and not family_id.startswith("#"):
            output[family_id] = as_float(row[column])
    return output


def read_family_matrix(path: Path, target: str) -> Tuple[List[str], Dict[str, Tuple[str, int]]]:
    header, rows = read_tsv(path)
    if len(header) < 3:
        raise ValueError("Family matrix must have Desc, Family ID and species columns")
    target_column = locate_target_column(header, target)
    result: Dict[str, Tuple[str, int]] = {}
    order: List[str] = []
    for line_number, row in enumerate(rows, start=2):
        if len(row) != len(header):
            raise ValueError(
                f"Family matrix line {line_number}: expected {len(header)} fields, got {len(row)}"
            )
        family_id = row[1].strip()
        if family_id in result:
            raise ValueError(f"Duplicate family ID in matrix: {family_id}")
        try:
            observed = int(row[target_column])
        except ValueError as exc:
            raise ValueError(
                f"Non-integer target count for {family_id}: {row[target_column]!r}"
            ) from exc
        result[family_id] = (row[0].strip(), observed)
        order.append(family_id)
    return order, result


def read_family_results(path: Path) -> Dict[str, Tuple[Optional[float], str]]:
    header, rows = read_tsv(path)
    output: Dict[str, Tuple[Optional[float], str]] = {}
    for row in rows:
        if len(row) < 2:
            continue
        family_id = row[0].strip()
        if not family_id or family_id.startswith("#"):
            continue
        pvalue = as_float(row[1])
        flag = row[2].strip() if len(row) > 2 else ""
        output[family_id] = (pvalue, flag)
    return output


def read_clade_expected(path: Path, target: str) -> Optional[Tuple[int, int]]:
    if not path.is_file() or path.stat().st_size == 0:
        return None
    _, rows = read_tsv(path)
    matches = [row for row in rows if row and species_part(row[0]) == target]
    if len(matches) != 1 or len(matches[0]) < 3:
        return None
    try:
        return int(float(matches[0][1])), int(float(matches[0][2]))
    except ValueError:
        return None


def format_value(value: Optional[float]) -> str:
    if value is None:
        return "NA"
    if abs(value - round(value)) < 1e-9:
        return str(int(round(value)))
    return f"{value:.12g}"


def write_rows(path: Path, header: Sequence[str], rows: Sequence[Sequence[str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    if not 0 < args.strict_p <= args.branch_p <= 1:
        raise ValueError("Require 0 < --strict-p <= --branch-p <= 1")
    args.outdir.mkdir(parents=True, exist_ok=True)

    prefix = args.model
    change_path = args.results / f"{prefix}_change.tab"
    count_path = args.results / f"{prefix}_count.tab"
    probability_path = args.results / f"{prefix}_branch_probabilities.tab"
    family_results_path = args.results / f"{prefix}_family_results.txt"
    clade_results_path = args.results / f"{prefix}_clade_results.txt"

    order, matrix = read_family_matrix(args.families, args.target)
    changes = read_target_vector(change_path, args.target)
    counts = read_target_vector(count_path, args.target)
    probabilities = read_target_vector(probability_path, args.target)
    family_results = read_family_results(family_results_path)
    expected = read_clade_expected(clade_results_path, args.target)

    header = [
        "Family_ID", "Desc", "target_observed_count",
        "target_reconstructed_count", "parent_reconstructed_count",
        "change", "direction", "family_wide_p", "family_wide_flag",
        "branch_probability", f"branch_significant_p_le_{args.branch_p:g}",
        f"branch_strict_p_le_{args.strict_p:g}",
    ]
    all_rows: List[List[str]] = []
    missing_rows: List[List[str]] = []

    for family_id in order:
        description, observed = matrix[family_id]
        change = changes.get(family_id)
        reconstructed = counts.get(family_id)
        branch_p = probabilities.get(family_id)
        family_p, family_flag = family_results.get(family_id, (None, ""))

        if change is None or reconstructed is None:
            missing_rows.append([family_id, description, str(observed)])
            continue
        parent_count = reconstructed - change
        direction = "expanded" if change > 0 else "contracted" if change < 0 else "unchanged"
        row = [
            family_id,
            description,
            str(observed),
            format_value(reconstructed),
            format_value(parent_count),
            format_value(change),
            direction,
            format_value(family_p),
            family_flag or "NA",
            format_value(branch_p),
            "yes" if branch_p is not None and branch_p <= args.branch_p else "no",
            "yes" if branch_p is not None and branch_p <= args.strict_p else "no",
        ]
        all_rows.append(row)

    changed = [row for row in all_rows if row[6] != "unchanged"]
    expanded = [row for row in changed if row[6] == "expanded"]
    contracted = [row for row in changed if row[6] == "contracted"]
    significant = [row for row in changed if row[10] == "yes"]
    strict = [row for row in changed if row[11] == "yes"]

    outputs = {
        "all": all_rows,
        "changed": changed,
        "expanded": expanded,
        "contracted": contracted,
        "significant": significant,
        "strict": strict,
    }
    for name, rows in outputs.items():
        write_rows(args.outdir / f"{args.target}.{name}.tsv", header, rows)
    write_rows(
        args.outdir / f"{args.target}.unreconstructed.tsv",
        ["Family_ID", "Desc", "target_observed_count"],
        missing_rows,
    )

    observed_counts = (len(expanded), len(contracted))
    qc_status = "not_available"
    if expected is not None:
        qc_status = "pass" if observed_counts == expected else "WARNING_mismatch"
    summary = [
        ("model", prefix),
        ("target", args.target),
        ("families_in_input", str(len(order))),
        ("families_reconstructed", str(len(all_rows))),
        ("families_unreconstructed", str(len(missing_rows))),
        ("expanded_all", str(len(expanded))),
        ("contracted_all", str(len(contracted))),
        (f"changed_branch_p_le_{args.branch_p:g}", str(len(significant))),
        (f"changed_branch_p_le_{args.strict_p:g}", str(len(strict))),
        ("clade_results_expected_increase", str(expected[0]) if expected else "NA"),
        ("clade_results_expected_decrease", str(expected[1]) if expected else "NA"),
        ("change_vs_clade_results_qc", qc_status),
        ("note", "family_wide_p and branch_probability are different tests"),
    ]
    write_rows(args.outdir / "focal_summary.tsv", ["metric", "value"], summary)

    print(f"Target: {args.target}")
    print(f"Expanded / contracted: {len(expanded)} / {len(contracted)}")
    print(f"Branch P <= {args.branch_p:g}: {len(significant)}")
    print(f"Branch P <= {args.strict_p:g}: {len(strict)}")
    print(f"Unreconstructed: {len(missing_rows)}")
    print(f"QC against clade_results: {qc_status}")
    print(f"Output: {args.outdir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)

