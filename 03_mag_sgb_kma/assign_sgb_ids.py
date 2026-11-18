#!/usr/bin/env python3
"""Assign deterministic SGB IDs to dRep secondary clusters and representatives."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path


def stem(name: str) -> str:
    for suffix in (".fasta", ".fna", ".fa"):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return name


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--drep", required=True, type=Path)
    parser.add_argument("--hq-manifest", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args()
    reps_out = args.output_root / "representatives"
    membership_out = args.output_root / "sgb_membership.tsv"
    reps_table = args.output_root / "sgb_representatives.tsv"
    if reps_table.is_file() and reps_table.stat().st_size > 0:
        print(f"[REUSE] {reps_table}")
        return
    if reps_out.exists() and any(reps_out.iterdir()):
        raise SystemExit(f"Representative directory exists without completed table: {reps_out}")
    reps_out.mkdir(parents=True, exist_ok=True)

    with (args.drep / "data_tables" / "Cdb.csv").open(newline="", encoding="utf-8") as handle:
        clusters = list(csv.DictReader(handle))
    if not clusters or "genome" not in clusters[0] or "secondary_cluster" not in clusters[0]:
        raise SystemExit("Unexpected or empty dRep Cdb.csv")
    with args.hq_manifest.open(newline="", encoding="utf-8") as handle:
        hq_rows = list(csv.DictReader(handle, delimiter="\t"))
    hq = {stem(row["filename"]): row for row in hq_rows if row["selected_hq"].casefold() == "true"}

    winners = sorted(p for p in (args.drep / "dereplicated_genomes").iterdir() if p.is_file())
    winner_by_stem = {stem(p.name): p for p in winners}
    cluster_members = {}
    for row in clusters:
        cluster_members.setdefault(row["secondary_cluster"], []).append(row["genome"])
    winner_cluster = {}
    for cluster, members in cluster_members.items():
        matching = [m for m in members if stem(Path(m).name) in winner_by_stem]
        if len(matching) != 1:
            raise SystemExit(f"Expected one representative for cluster {cluster}; found {matching}")
        winner_cluster[cluster] = matching[0]

    ordered_clusters = sorted(cluster_members, key=lambda x: (str(x)))
    membership_rows = []
    representative_rows = []
    for number, cluster in enumerate(ordered_clusters, start=1):
        sgb_id = f"SGB_{number:06d}"
        winner = winner_cluster[cluster]
        winner_key = stem(Path(winner).name)
        winner_path = winner_by_stem[winner_key]
        target = reps_out / f"{sgb_id}.fa"
        target.symlink_to(winner_path.resolve())
        source_meta = hq.get(winner_key, {})
        representative_rows.append({
            "sgb_id": sgb_id,
            "drep_secondary_cluster": cluster,
            "representative_mag": winner_key,
            "representative_path": str(target),
            "group": source_meta.get("group", "NA"),
            "group_slug": source_meta.get("group_slug", "NA"),
            "completeness": source_meta.get("completeness", "NA"),
            "contamination": source_meta.get("contamination", "NA"),
            "gunc_pass": source_meta.get("gunc_pass", "NA"),
        })
        for member in sorted(cluster_members[cluster]):
            member_key = stem(Path(member).name)
            member_meta = hq.get(member_key, {})
            membership_rows.append({
                "sgb_id": sgb_id,
                "drep_secondary_cluster": cluster,
                "mag_id": member_key,
                "is_representative": str(member_key == winner_key),
                "group": member_meta.get("group", "NA"),
                "group_slug": member_meta.get("group_slug", "NA"),
            })

    for path, rows in ((membership_out, membership_rows), (reps_table, representative_rows)):
        tmp = path.with_name(path.name + f".part.{os.getpid()}")
        with tmp.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        tmp.replace(path)
    print(f"[SUCCESS] sgb_clusters={len(representative_rows)} representatives={reps_out}")


if __name__ == "__main__":
    main()

