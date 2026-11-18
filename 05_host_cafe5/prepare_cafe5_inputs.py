#!/usr/bin/env python3
"""Prepare a primate-only TOGA2 gene-family matrix and matched CAFE5 time tree.

The input family table must be a CAFE-style tab-delimited matrix whose first
two columns are description and family ID. Species start in column three.

The script deliberately starts from the unfiltered family table, removes the
requested outgroups, matches the remaining species to a large dated tree,
prunes and renames that tree, and then reapplies family filters.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

try:
    from Bio import Phylo
except ImportError as exc:
    raise SystemExit(
        "Biopython is required. Install it with: conda install -c conda-forge biopython"
    ) from exc


@dataclass
class Family:
    description: str
    family_id: str
    counts: List[int]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare primate-only TOGA2 counts and a matched CAFE5 tree."
    )
    parser.add_argument(
        "--families",
        required=True,
        type=Path,
        help="Unfiltered CAFE matrix, normally cafe5_gene_families.all.txt",
    )
    parser.add_argument(
        "--tree",
        required=True,
        type=Path,
        help="Large rooted, binary, ultrametric time tree in Newick format",
    )
    parser.add_argument(
        "--outdir", required=True, type=Path, help="Output directory"
    )
    parser.add_argument(
        "--drop",
        nargs="*",
        default=[],
        metavar="SPECIES",
        help=(
            "Optional species columns to remove. By default no species are removed."
        ),
    )
    parser.add_argument(
        "--tip-map",
        type=Path,
        help=(
            "Optional two-column TSV with header matrix_name,tree_tip. "
            "Manual mappings override automatic matching."
        ),
    )
    parser.add_argument(
        "--max-family-size",
        type=int,
        default=100,
        help="Separate families exceeding this copy number (default: 100)",
    )
    parser.add_argument(
        "--min-present-taxa",
        type=int,
        default=2,
        help="Minimum number of non-zero taxa in the general filtered set (default: 2)",
    )
    parser.add_argument(
        "--ultrametric-tolerance",
        type=float,
        default=1e-5,
        help="Relative root-to-tip distance tolerance (default: 1e-5)",
    )
    parser.add_argument(
        "--repair-rounded-ultrametric",
        action="store_true",
        help=(
            "For an otherwise valid dated tree with small rounding differences, "
            "extend terminal branches to the longest root-to-tip distance"
        ),
    )
    parser.add_argument(
        "--max-ultrametric-repair-ratio",
        type=float,
        default=5e-4,
        help=(
            "Maximum root-to-tip spread divided by tree height allowed for "
            "terminal-branch repair (default: 5e-4)"
        ),
    )
    parser.add_argument(
        "--focal",
        default="Cephalopachus_bancanus",
        help="Focal species required in the final matrix/tree",
    )
    return parser.parse_args()


def normalize_name(name: str) -> str:
    name = name.strip().strip("'\"")
    name = re.sub(r"\.prot$", "", name)
    name = re.sub(r"[\s.-]+", "_", name)
    name = re.sub(r"_+", "_", name)
    return name.strip("_")


def binomial_key(name: str) -> Optional[Tuple[str, str]]:
    parts = normalize_name(name).split("_")
    if len(parts) < 2:
        return None
    return parts[0].lower(), parts[1].lower()


def read_families(path: Path) -> Tuple[List[str], List[Family]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"Empty family table: {path}") from exc

        header = [field.lstrip("\ufeff").strip() for field in header]
        if len(header) < 3:
            raise ValueError("Family table must contain at least three columns")
        if len(set(header[2:])) != len(header[2:]):
            raise ValueError("Duplicate species columns found in family table")

        species = header[2:]
        families: List[Family] = []
        seen_ids: Set[str] = set()

        for line_number, row in enumerate(reader, start=2):
            if not row or all(not value.strip() for value in row):
                continue
            if len(row) != len(header):
                raise ValueError(
                    f"Line {line_number}: expected {len(header)} columns, got {len(row)}"
                )
            description, family_id = row[0].strip(), row[1].strip()
            if not family_id:
                raise ValueError(f"Line {line_number}: empty family ID")
            if family_id in seen_ids:
                raise ValueError(f"Line {line_number}: duplicate family ID {family_id}")
            seen_ids.add(family_id)

            try:
                counts = [int(value) for value in row[2:]]
            except ValueError as exc:
                raise ValueError(
                    f"Line {line_number}: all family counts must be integers"
                ) from exc
            if any(value < 0 for value in counts):
                raise ValueError(f"Line {line_number}: negative family count")

            families.append(Family(description, family_id, counts))

    return species, families


def read_manual_map(path: Optional[Path]) -> Dict[str, str]:
    if path is None:
        return {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"matrix_name", "tree_tip"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(
                "Tip map requires TSV header: matrix_name<TAB>tree_tip"
            )
        mapping = {}
        for row in reader:
            matrix_name = row["matrix_name"].strip()
            tree_tip = row["tree_tip"].strip()
            if matrix_name and tree_tip:
                mapping[matrix_name] = tree_tip
    return mapping


def match_species_to_tips(
    matrix_species: Sequence[str], tree_tips: Sequence[str], manual: Dict[str, str]
) -> Tuple[Dict[str, Tuple[str, str]], List[str]]:
    tip_set = set(tree_tips)
    normalized_index: Dict[str, List[str]] = {}
    binomial_index: Dict[Tuple[str, str], List[str]] = {}
    for tip in tree_tips:
        normalized_index.setdefault(normalize_name(tip).lower(), []).append(tip)
        key = binomial_key(tip)
        if key is not None:
            binomial_index.setdefault(key, []).append(tip)

    mapping: Dict[str, Tuple[str, str]] = {}
    unmatched: List[str] = []
    used_tips: Dict[str, str] = {}

    for species in matrix_species:
        chosen: Optional[str] = None
        method = ""

        if species in manual:
            candidate = manual[species]
            if candidate not in tip_set:
                raise ValueError(
                    f"Manual tree tip for {species!r} does not exist: {candidate!r}"
                )
            chosen, method = candidate, "manual"
        elif species in tip_set:
            chosen, method = species, "exact"
        else:
            normalized_candidates = normalized_index.get(
                normalize_name(species).lower(), []
            )
            if len(normalized_candidates) == 1:
                chosen, method = normalized_candidates[0], "normalized_exact"
            else:
                key = binomial_key(species)
                candidates = binomial_index.get(key, []) if key else []
                if len(candidates) == 1:
                    chosen, method = candidates[0], "unique_genus_species"

        if chosen is None:
            unmatched.append(species)
            continue
        if chosen in used_tips:
            raise ValueError(
                f"Tree tip {chosen!r} would be assigned to both "
                f"{used_tips[chosen]!r} and {species!r}. Use --tip-map to resolve it."
            )
        used_tips[chosen] = species
        mapping[species] = (chosen, method)

    return mapping, unmatched


def collapse_unary_root(tree) -> None:
    while len(tree.root.clades) == 1:
        child = tree.root.clades[0]
        tree.root = child


def validate_tree(tree, tolerance: float) -> Dict[str, float]:
    terminals = tree.get_terminals()
    if len(terminals) < 2:
        raise ValueError("The pruned tree has fewer than two tips")

    nonbinary = [clade for clade in tree.get_nonterminals() if len(clade.clades) != 2]
    if nonbinary:
        raise ValueError(
            f"The pruned tree is not binary: {len(nonbinary)} internal nodes do not have two children"
        )

    missing_lengths = []
    negative_lengths = []
    for clade in tree.find_clades(order="level"):
        if clade is tree.root:
            continue
        if clade.branch_length is None:
            missing_lengths.append(clade.name or "<internal>")
        elif clade.branch_length < 0:
            negative_lengths.append(clade.name or "<internal>")
    if missing_lengths:
        raise ValueError(
            f"Tree contains {len(missing_lengths)} branches without branch lengths"
        )
    if negative_lengths:
        raise ValueError(
            f"Tree contains {len(negative_lengths)} negative branch lengths"
        )

    distances = [tree.distance(tree.root, tip) for tip in terminals]
    minimum = min(distances)
    maximum = max(distances)
    allowed = tolerance * max(1.0, abs(maximum))
    spread = maximum - minimum
    if spread > allowed:
        raise ValueError(
            "Tree is not ultrametric within tolerance: "
            f"min root-tip={minimum:.10g}, max root-tip={maximum:.10g}, "
            f"spread={spread:.10g}, allowed={allowed:.10g}"
        )

    return {
        "tips": float(len(terminals)),
        "root_to_tip_min": minimum,
        "root_to_tip_max": maximum,
        "root_to_tip_spread": spread,
    }


def repair_rounded_ultrametric(tree, maximum_ratio: float) -> List[Tuple[str, float, float, float]]:
    """Make a nearly ultrametric rounded time tree exact by extending tip branches.

    This does not alter the topology or any internal branch.  It is deliberately
    limited to a small relative discrepancy so that a genuinely non-ultrametric
    tree cannot be silently converted.
    """
    if maximum_ratio < 0:
        raise ValueError("--max-ultrametric-repair-ratio must be non-negative")

    terminals = tree.get_terminals()
    distances = {tip: tree.distance(tree.root, tip) for tip in terminals}
    maximum = max(distances.values())
    minimum = min(distances.values())
    spread = maximum - minimum
    ratio = spread / max(1.0, abs(maximum))
    if ratio > maximum_ratio:
        raise ValueError(
            "Refusing ultrametric repair because the discrepancy is too large: "
            f"min root-tip={minimum:.10g}, max root-tip={maximum:.10g}, "
            f"spread={spread:.10g}, relative spread={ratio:.10g}, "
            f"maximum allowed={maximum_ratio:.10g}"
        )

    repairs: List[Tuple[str, float, float, float]] = []
    for tip in terminals:
        before = distances[tip]
        addition = maximum - before
        if addition < 1e-12:
            addition = 0.0
        else:
            tip.branch_length = float(tip.branch_length) + addition
        repairs.append((tip.name, before, addition, before + addition))
    return repairs


def write_family_table(path: Path, species: Sequence[str], families: Iterable[Family]) -> int:
    count = 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Desc", "Family ID", *species])
        for family in families:
            writer.writerow(
                [family.description, family.family_id, *family.counts]
            )
            count += 1
    return count


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    species, input_families = read_families(args.families)
    drop_set = set(args.drop)
    unknown_drops = sorted(drop_set.difference(species))
    if unknown_drops:
        print(
            "WARNING: requested drop species not found in matrix: "
            + ", ".join(unknown_drops),
            file=sys.stderr,
        )

    keep_indices = [index for index, name in enumerate(species) if name not in drop_set]
    kept_species = [species[index] for index in keep_indices]
    if args.focal not in kept_species:
        raise ValueError(
            f"Focal species {args.focal!r} is absent after outgroup removal"
        )

    families = [
        Family(
            family.description,
            family.family_id,
            [family.counts[index] for index in keep_indices],
        )
        for family in input_families
    ]

    tree = Phylo.read(str(args.tree), "newick")
    original_tips = [tip.name for tip in tree.get_terminals()]
    if any(name is None for name in original_tips):
        raise ValueError("Tree contains unnamed tips")

    manual_map = read_manual_map(args.tip_map)
    mapping, unmatched = match_species_to_tips(
        kept_species, original_tips, manual_map
    )

    unmatched_path = args.outdir / "unmatched_matrix_species.txt"
    unmatched_path.write_text(
        "".join(f"{name}\n" for name in unmatched), encoding="utf-8"
    )
    if unmatched:
        raise ValueError(
            f"{len(unmatched)} matrix species could not be matched to the tree. "
            f"See {unmatched_path} and provide --tip-map."
        )

    keep_tree_tips = {tree_tip for tree_tip, _ in mapping.values()}
    for terminal in list(tree.get_terminals()):
        if terminal.name not in keep_tree_tips:
            tree.prune(terminal)
    collapse_unary_root(tree)

    reverse_mapping = {
        tree_tip: matrix_name
        for matrix_name, (tree_tip, method) in mapping.items()
    }
    for terminal in tree.get_terminals():
        terminal.name = reverse_mapping[terminal.name]
    tree.rooted = True

    # First validate topology and branch lengths without rejecting small clock
    # discrepancies.  Rounded published time trees commonly accumulate a few
    # hundredths of a time unit along paths containing different numbers of edges.
    validate_tree(tree, float("inf"))
    repairs: List[Tuple[str, float, float, float]] = []
    if args.repair_rounded_ultrametric:
        repairs = repair_rounded_ultrametric(
            tree, args.max_ultrametric_repair_ratio
        )
    metrics = validate_tree(tree, args.ultrametric_tolerance)
    final_tree_tips = {tip.name for tip in tree.get_terminals()}
    if final_tree_tips != set(kept_species):
        raise ValueError("Internal error: final tree tips do not equal matrix species")

    root_children = tree.root.clades
    left_taxa = {tip.name for tip in root_children[0].get_terminals()}
    right_taxa = {tip.name for tip in root_children[1].get_terminals()}
    species_index = {name: index for index, name in enumerate(kept_species)}

    sparse: List[Family] = []
    large: List[Family] = []
    filtered: List[Family] = []
    root_supported: List[Family] = []
    root_unsupported: List[Family] = []

    for family in families:
        present_taxa = sum(value > 0 for value in family.counts)
        maximum = max(family.counts)
        if maximum > args.max_family_size:
            large.append(family)
            continue
        if present_taxa < args.min_present_taxa:
            sparse.append(family)
            continue

        filtered.append(family)
        left_present = any(
            family.counts[species_index[name]] > 0 for name in left_taxa
        )
        right_present = any(
            family.counts[species_index[name]] > 0 for name in right_taxa
        )
        if left_present and right_present:
            root_supported.append(family)
        else:
            root_unsupported.append(family)

    outputs = {
        "all": args.outdir / "cafe5_primate.all.txt",
        "filtered": args.outdir / "cafe5_primate.filtered.txt",
        "root_supported": args.outdir / "cafe5_primate.root_supported.txt",
        "root_unsupported": args.outdir / "cafe5_primate.root_unsupported.txt",
        "large": args.outdir / "cafe5_primate.large.txt",
        "sparse": args.outdir / "cafe5_primate.sparse.txt",
    }
    write_family_table(outputs["all"], kept_species, families)
    write_family_table(outputs["filtered"], kept_species, filtered)
    write_family_table(outputs["root_supported"], kept_species, root_supported)
    write_family_table(outputs["root_unsupported"], kept_species, root_unsupported)
    write_family_table(outputs["large"], kept_species, large)
    write_family_table(outputs["sparse"], kept_species, sparse)

    tree_path = args.outdir / "primate_timetree.nwk"
    Phylo.write(
        tree,
        str(tree_path),
        "newick",
        format_branch_length="%1.10f",
    )
    written_tree = Phylo.read(str(tree_path), "newick")
    written_metrics = validate_tree(written_tree, args.ultrametric_tolerance)
    written_tips = {tip.name for tip in written_tree.get_terminals()}
    if written_tips != set(kept_species):
        raise ValueError("Written tree tips do not equal matrix species")

    mapping_path = args.outdir / "tree_tip_mapping.tsv"
    with mapping_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["matrix_name", "original_tree_tip", "method"])
        for name in kept_species:
            tree_tip, method = mapping[name]
            writer.writerow([name, tree_tip, method])

    repair_path = args.outdir / "ultrametric_terminal_repairs.tsv"
    with repair_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "species",
                "root_to_tip_before",
                "terminal_branch_added",
                "root_to_tip_after",
            ]
        )
        for name, before, addition, after in sorted(repairs):
            writer.writerow(
                [name, f"{before:.12g}", f"{addition:.12g}", f"{after:.12g}"]
            )

    root_clades_path = args.outdir / "root_daughter_clades.tsv"
    with root_clades_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["root_side", "species"])
        for name in sorted(left_taxa):
            writer.writerow(["A", name])
        for name in sorted(right_taxa):
            writer.writerow(["B", name])

    summary_lines = [
        f"source_family_table\t{args.families.resolve()}",
        f"source_tree\t{args.tree.resolve()}",
        f"dropped_species\t{','.join(args.drop)}",
        f"retained_taxa\t{len(kept_species)}",
        f"focal_species\t{args.focal}",
        f"input_families\t{len(input_families)}",
        f"all_after_taxon_removal\t{len(families)}",
        f"general_filtered\t{len(filtered)}",
        f"root_supported\t{len(root_supported)}",
        f"root_unsupported\t{len(root_unsupported)}",
        f"large_gt_{args.max_family_size}\t{len(large)}",
        f"sparse_lt_{args.min_present_taxa}_taxa\t{len(sparse)}",
        f"root_to_tip_min\t{written_metrics['root_to_tip_min']:.12g}",
        f"root_to_tip_max\t{written_metrics['root_to_tip_max']:.12g}",
        f"root_to_tip_spread\t{written_metrics['root_to_tip_spread']:.12g}",
        f"ultrametric_terminal_repair\t{args.repair_rounded_ultrametric}",
        f"ultrametric_repaired_tips\t{sum(addition > 0 for _, _, addition, _ in repairs)}",
        f"ultrametric_max_terminal_addition\t{max((addition for _, _, addition, _ in repairs), default=0):.12g}",
        f"ultrametric_repair_report\t{repair_path.resolve()}",
        f"recommended_primary_input\t{outputs['root_supported'].resolve()}",
        f"sensitivity_input\t{outputs['filtered'].resolve()}",
        f"tree_output\t{tree_path.resolve()}",
    ]
    summary_path = args.outdir / "prepare_summary.tsv"
    summary_path.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    print(f"Prepared {len(kept_species)} taxa and {len(families)} families")
    print(f"Removed taxa: {', '.join(args.drop) if args.drop else '(none)'}")
    print(f"General filtered families: {len(filtered)}")
    print(f"Root-supported families: {len(root_supported)}")
    print(f"Primary CAFE5 input: {outputs['root_supported']}")
    print(f"Sensitivity input: {outputs['filtered']}")
    print(f"Matched time tree: {tree_path}")
    print(f"Summary: {summary_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
