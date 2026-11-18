#!/usr/bin/env python3
import argparse
import csv
import gzip
import os
import shutil
import tempfile
from collections import Counter, defaultdict
from pathlib import Path


def number(value):
    try:
        return f"{float(value):.12g}"
    except (TypeError, ValueError):
        return value


def write_tsv(path, header, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    out = args.out.resolve()

    with (out / "input_manifest.tsv").open(newline="", encoding="utf-8") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    mag_to_sgb = {row["mag_id"]: row["sgb_id"] for row in manifest}
    if len(mag_to_sgb) != len(manifest):
        raise ValueError("Duplicate mag_id values in input_manifest.tsv")

    sgb_members = defaultdict(set)
    for mag, sgb in mag_to_sgb.items():
        sgb_members[sgb].add(mag)

    destination = out / "12_final_tables" / "fegenie"
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix="fegenie.tmp.", dir=destination.parent))

    category_members = defaultdict(set)
    category_orfs = Counter()
    hmm_members = defaultdict(set)
    hmm_orfs = Counter()
    any_members = defaultdict(set)
    any_orfs = Counter()
    mag_category = Counter()
    mag_any = Counter()
    gene_rows = 0

    gene_path = temporary / "gene_x_fegenie.tsv.gz"
    with gzip.open(gene_path, "wt", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "sgb_id", "mag_id", "fegenie_orf_id", "category",
            "related_hmm", "hmm_bitscore", "gene_cluster_id",
        ])

        for mag, sgb in mag_to_sgb.items():
            result = out / "07_fegenie" / mag
            if not (result / ".done").is_file():
                raise ValueError(f"Incomplete FeGenie result: {mag}")
            source = result / "FinalSummary-dereplicated-clustered.csv"
            if not source.is_file():
                raise ValueError(f"Missing FeGenie summary: {source}")

            seen_orfs = set()
            with source.open(newline="", encoding="utf-8", errors="replace") as input_handle:
                for row in csv.reader(input_handle):
                    if not row or row[0] == "#" or row[0].lower() == "category":
                        continue
                    if len(row) < 6:
                        raise ValueError(f"Malformed FeGenie row in {source}: {row}")
                    category, cell, orf_id, hmm, bitscore, cluster_id = row[:6]
                    if orf_id in seen_orfs:
                        raise ValueError(f"Duplicate dereplicated ORF in {source}: {orf_id}")
                    seen_orfs.add(orf_id)
                    writer.writerow([
                        sgb, mag, orf_id, category, hmm,
                        number(bitscore), cluster_id,
                    ])
                    gene_rows += 1
                    category_members[(sgb, category)].add(mag)
                    category_orfs[(sgb, category)] += 1
                    hmm_members[(sgb, hmm)].add(mag)
                    hmm_orfs[(sgb, hmm)] += 1
                    any_members[sgb].add(mag)
                    any_orfs[sgb] += 1
                    mag_category[(sgb, mag, category)] += 1
                    mag_any[(sgb, mag)] += 1

    def prevalence_rows(member_map, count_map):
        result = []
        for (sgb, feature), members in sorted(member_map.items()):
            total = len(sgb_members[sgb])
            count = count_map[(sgb, feature)]
            result.append([
                sgb, feature, len(members), total,
                number(len(members) / total), count, number(count / total),
            ])
        return result

    long_header = [
        "sgb_id", "feature", "members_with_feature", "members_total",
        "prevalence", "orf_count_across_members", "mean_orfs_per_member",
    ]
    write_tsv(
        temporary / "sgb_x_fegenie_category_long.tsv",
        [value.replace("feature", "category") for value in long_header],
        prevalence_rows(category_members, category_orfs),
    )
    write_tsv(
        temporary / "sgb_x_fegenie_hmm_long.tsv",
        [value.replace("feature", "hmm") for value in long_header],
        prevalence_rows(hmm_members, hmm_orfs),
    )

    any_rows = []
    for sgb in sorted(sgb_members):
        total = len(sgb_members[sgb])
        members = len(any_members[sgb])
        count = any_orfs[sgb]
        any_rows.append([
            sgb, members, total, number(members / total),
            count, number(count / total),
        ])
    write_tsv(
        temporary / "sgb_x_fegenie_any.tsv",
        [
            "sgb_id", "members_with_any_hit", "members_total", "prevalence",
            "orf_count_across_members", "mean_orfs_per_member",
        ],
        any_rows,
    )

    mag_category_rows = [
        [sgb, mag, category, count]
        for (sgb, mag, category), count in sorted(mag_category.items())
    ]
    write_tsv(
        temporary / "mag_x_fegenie_category_long.tsv",
        ["sgb_id", "mag_id", "category", "orf_count"],
        mag_category_rows,
    )
    write_tsv(
        temporary / "mag_x_fegenie_any.tsv",
        ["sgb_id", "mag_id", "orf_count"],
        [[row["sgb_id"], row["mag_id"], mag_any[(row["sgb_id"], row["mag_id"])]]
         for row in manifest],
    )

    with (temporary / ".done").open("w", encoding="utf-8") as handle:
        handle.write("status=SUCCESS\n")
        handle.write(f"MAGs={len(manifest)}\n")
        handle.write(f"SGBs={len(sgb_members)}\n")
        handle.write(f"gene_rows={gene_rows}\n")

    if destination.exists():
        shutil.rmtree(destination)
    os.replace(temporary, destination)
    print(
        f"[SUCCESS] MAGs={len(manifest)} SGBs={len(sgb_members)} "
        f"FeGenie_hits={gene_rows}"
    )


if __name__ == "__main__":
    main()
