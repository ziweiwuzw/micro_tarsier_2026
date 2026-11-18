#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import pandas as pd

BAD = {"", "nan", "none", "na", "n/a", "not", "not provided", "missing", "-"}

SPECIES_SYNONYMS = {
    "Xanthonycticebus_pygmaeus": "Nycticebus_pygmaeus",
    "Cephalophacus_bancanus": "Cephalopachus_bancanus",
}

MACHADO_SYNONYMS = {
    "xanthonycticebus pygmaeus": "nycticebus pygmaeus",
    "cephalophacus bancanus": "cephalopachus bancanus",
    "lagothrix lagotricha": "lagothrix lagothricha",
    "rhinopithecus roxellanae": "rhinopithecus roxellana",
}

MANUAL_DIET = {
    "Cephalopachus_bancanus": "Insectivore",
}

FASTQ_RE = re.compile(
    r"^(?P<sample>.+?)(?:_R(?P<r>[12])(?:_001)?|_(?P<u>[12]))"
    r"\.nonhost\.(?:fastq|fq)(?:\.gz)?$",
    re.IGNORECASE,
)


def clean(value) -> str:
    if pd.isna(value):
        return ""

    text = str(value).strip()

    if text.lower() in BAD:
        return ""

    return text


def canonical_species(value) -> str:
    text = clean(value).replace(" ", "_")
    text = re.sub(r"[^A-Za-z0-9_]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")

    text = SPECIES_SYNONYMS.get(text, text)

    parts = text.split("_")

    if len(parts) >= 2:
        return "_".join(parts[:2])

    return text


def machado_key(value) -> str:
    text = clean(value).replace("_", " ")
    text = re.sub(r"\s+", " ", text).strip().lower()

    text = MACHADO_SYNONYMS.get(text, text)

    parts = text.split()

    if len(parts) >= 2:
        return " ".join(parts[:2])

    return text


def join_unique(values) -> str:
    out = []

    for value in values:
        for token in clean(value).split(";"):
            token = token.strip()

            if token and token not in out:
                out.append(token)

    return "; ".join(out)


def pick_column(
    df: pd.DataFrame,
    candidates: list[str],
) -> str | None:

    lookup = {
        re.sub(r"[^a-z0-9]+", "", str(col).lower()): col
        for col in df.columns
    }

    for candidate in candidates:
        key = re.sub(
            r"[^a-z0-9]+",
            "",
            candidate.lower(),
        )

        if key in lookup:
            return lookup[key]

    return None


def read_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in {".tsv", ".txt"}:
        return pd.read_csv(
            path,
            sep="\t",
            dtype=str,
        )

    return pd.read_csv(
        path,
        dtype=str,
    )


def load_species_meta(
    path: Path | None,
) -> pd.DataFrame:

    columns = [
        "species_merged",
        "metadata_species",
        "diet_species",
        "diet_species_source",
    ]

    if path is None or not path.exists():
        return pd.DataFrame(columns=columns)

    df = read_table(path)
    df.columns = [
        str(c).strip()
        for c in df.columns
    ]

    species_col = pick_column(
        df,
        [
            "Binomial name",
            "Species",
            "Scientific name",
            "Primate species",
            "taxon_name",
        ],
    )

    diet_col = pick_column(
        df,
        [
            "Diet",
            "Trophic guild",
            "trophic_guild",
        ],
    )

    if species_col is None:
        raise SystemExit(
            f"ERROR: no species column in {path}"
        )

    out = pd.DataFrame(
        {
            "species_merged": (
                df[species_col]
                .map(canonical_species)
            ),
            "metadata_species": (
                df[species_col]
                .map(clean)
            ),
            "diet_species": (
                df[diet_col].map(clean)
                if diet_col
                else ""
            ),
            "diet_species_source": path.name,
        }
    )

    out = (
        out
        .groupby(
            "species_merged",
            dropna=False,
        )
        .agg(
            {
                "metadata_species": join_unique,
                "diet_species": join_unique,
                "diet_species_source": join_unique,
            }
        )
        .reset_index()
    )

    return out


def load_run_meta(
    path: Path | None,
) -> pd.DataFrame:

    columns = [
        "sample_id",
        "run_meta_species",
        "diet_run",
        "diet_run_source",
    ]

    if path is None or not path.exists():
        return pd.DataFrame(columns=columns)

    df = read_table(path)
    df.columns = [
        str(c).strip()
        for c in df.columns
    ]

    run_col = pick_column(
        df,
        [
            "run",
            "Run accession",
            "run_accession",
            "accession",
        ],
    )

    species_col = pick_column(
        df,
        [
            "taxon_name",
            "Taxon name",
            "host scientific name",
            "scientific_name",
        ],
    )

    diet_col = pick_column(
        df,
        [
            "diet",
            "Diet",
        ],
    )

    if run_col is None:
        raise SystemExit(
            f"ERROR: no run column in {path}"
        )

    out = pd.DataFrame(
        {
            "sample_id": (
                df[run_col]
                .map(clean)
            ),
            "run_meta_species": (
                df[species_col].map(clean)
                if species_col
                else ""
            ),
            "diet_run": (
                df[diet_col].map(clean)
                if diet_col
                else ""
            ),
            "diet_run_source": path.name,
        }
    )

    return out.drop_duplicates(
        "sample_id"
    )


def load_ddi(
    path: Path | None,
) -> pd.DataFrame:

    columns = [
        "ddi_key",
        "machado_species",
        "diet_diversity_DDI",
    ]

    if path is None or not path.exists():
        return pd.DataFrame(columns=columns)

    df = read_table(path)

    df.columns = [
        str(c).strip()
        for c in df.columns
    ]

    species_col = pick_column(
        df,
        [
            "Species",
            "Scientific name",
        ],
    )

    ddi_col = pick_column(
        df,
        [
            "Diet diversity (DDI)",
            "Diet diversity DDI",
            "DDI",
        ],
    )

    if species_col is None or ddi_col is None:
        raise SystemExit(
            "ERROR: Machado table needs Species "
            f"and Diet diversity (DDI): {path}"
        )

    out = pd.DataFrame(
        {
            "ddi_key": (
                df[species_col]
                .map(machado_key)
            ),
            "machado_species": (
                df[species_col]
                .map(clean)
            ),
            "diet_diversity_DDI": (
                df[ddi_col]
                .map(clean)
            ),
        }
    )

    out = out[
        out["ddi_key"].ne("")
    ].copy()

    return out.drop_duplicates(
        "ddi_key"
    )


def parse_done(
    path: Path,
) -> tuple[int | None, str]:

    if not path.exists():
        return None, "MISSING"

    try:
        values = {}

        for line in path.read_text(
            errors="replace"
        ).splitlines():

            if "=" not in line:
                continue

            key, value = line.split(
                "=",
                1,
            )

            values[key.strip()] = value.strip()

        if not values.get("nonhost_pairs"):
            return None, "UNPARSEABLE"

        return (
            int(values["nonhost_pairs"]),
            "OK",
        )

    except Exception as exc:
        return None, f"ERROR:{exc}"


def build_samples(
    stats_path: Path,
    reads_root: Path,
    require_done: bool,
) -> tuple[pd.DataFrame, list[str]]:

    stats = pd.read_csv(
        stats_path,
        sep="\t",
        dtype=str,
    )

    stats.columns = [
        str(c).strip()
        for c in stats.columns
    ]

    required = [
        "file",
        "num_seqs",
        "sum_len",
    ]

    missing = [
        c
        for c in required
        if c not in stats.columns
    ]

    if missing:
        raise SystemExit(
            "ERROR: seqkit table missing columns: "
            + ", ".join(missing)
        )

    rows = []
    skipped = []

    for _, row in stats.iterrows():

        path = Path(
            clean(row["file"])
        ).resolve()

        try:
            rel = path.relative_to(
                reads_root
            )

        except ValueError:
            skipped.append(str(path))
            continue

        if len(rel.parts) < 3:
            skipped.append(str(path))
            continue

        match = FASTQ_RE.match(
            path.name
        )

        if not match:
            skipped.append(str(path))
            continue

        num_seqs = clean(
            row["num_seqs"]
        ).replace(",", "")

        sum_len = clean(
            row["sum_len"]
        ).replace(",", "")

        rows.append(
            {
                "species_original": rel.parts[0],
                "species_merged": canonical_species(
                    rel.parts[0]
                ),
                "project": rel.parts[1],
                "sample_id": match.group("sample"),
                "mate": (
                    match.group("r")
                    or match.group("u")
                ),
                "file_path": str(path),
                "num_seqs": pd.to_numeric(
                    num_seqs,
                    errors="coerce",
                ),
                "sum_len": pd.to_numeric(
                    sum_len,
                    errors="coerce",
                ),
            }
        )

    file_df = pd.DataFrame(rows)

    if file_df.empty:
        raise SystemExit(
            "ERROR: no recognizable nonhost "
            "FASTQ rows in seqkit table"
        )

    sample_rows = []

    keys = [
        "species_original",
        "species_merged",
        "project",
        "sample_id",
    ]

    for (
        sp_original,
        sp_merged,
        project,
        sample,
    ), group in file_df.groupby(
        keys,
        dropna=False,
    ):

        r1 = group[
            group["mate"].eq("1")
        ]

        r2 = group[
            group["mate"].eq("2")
        ]

        duplicate = (
            len(r1) > 1
            or len(r2) > 1
        )

        r1_row = (
            r1.iloc[0]
            if len(r1)
            else None
        )

        r2_row = (
            r2.iloc[0]
            if len(r2)
            else None
        )

        r1_reads = (
            int(r1_row["num_seqs"])
            if (
                r1_row is not None
                and pd.notna(
                    r1_row["num_seqs"]
                )
            )
            else 0
        )

        r2_reads = (
            int(r2_row["num_seqs"])
            if (
                r2_row is not None
                and pd.notna(
                    r2_row["num_seqs"]
                )
            )
            else 0
        )

        r1_bases = (
            int(r1_row["sum_len"])
            if (
                r1_row is not None
                and pd.notna(
                    r1_row["sum_len"]
                )
            )
            else 0
        )

        r2_bases = (
            int(r2_row["sum_len"])
            if (
                r2_row is not None
                and pd.notna(
                    r2_row["sum_len"]
                )
            )
            else 0
        )

        sample_dir = (
            reads_root
            / str(sp_original)
            / str(project)
        )

        done_path = (
            sample_dir
            / f"{sample}.host_removed.done"
        )

        flagstat_path = (
            sample_dir
            / f"{sample}.host_mapping.flagstat.txt"
        )

        done_pairs, done_status = parse_done(
            done_path
        )

        if duplicate:
            status = "DUPLICATE_MATE_FILE"

        elif r1_row is None:
            status = "MISSING_R1"

        elif r2_row is None:
            status = "MISSING_R2"

        elif (
            pd.isna(r1_row["num_seqs"])
            or pd.isna(r2_row["num_seqs"])
        ):
            status = "SEQKIT_PARSE_ERROR"

        elif r1_reads != r2_reads:
            status = "R1_R2_READ_COUNT_MISMATCH"

        elif require_done and not done_path.exists():
            status = "DONE_MARKER_MISSING"

        elif (
            done_path.exists()
            and done_pairs is None
        ):
            status = "DONE_MARKER_UNPARSEABLE"

        elif (
            done_path.exists()
            and done_pairs != r1_reads
        ):
            status = "DONE_PAIR_COUNT_MISMATCH"

        else:
            status = "COMPLETE_EQUAL"

        reads_total = (
            r1_reads
            + r2_reads
        )

        bases_total = (
            r1_bases
            + r2_bases
        )

        sample_rows.append(
            {
                "species_original": sp_original,
                "species_merged": sp_merged,
                "project": project,
                "sample_id": sample,
                "layout": "PE",

                "r1_reads": r1_reads,
                "r2_reads": r2_reads,

                "nonhost_read_pairs": (
                    r1_reads
                    if r1_reads == r2_reads
                    else ""
                ),

                "nonhost_reads_total": reads_total,

                "r1_bases": r1_bases,
                "r2_bases": r2_bases,

                "nonhost_bases_total_bp": bases_total,

                "nonhost_size_Gb": round(
                    bases_total / 1e9,
                    6,
                ),

                "mean_read_length_bp": (
                    round(
                        bases_total / reads_total,
                        3,
                    )
                    if reads_total
                    else ""
                ),

                "r1_r2_read_difference": (
                    r1_reads
                    - r2_reads
                ),

                "host_removed_done": (
                    "YES"
                    if done_path.exists()
                    else "NO"
                ),

                "done_marker_status": done_status,

                "done_nonhost_pairs": (
                    done_pairs
                    if done_pairs is not None
                    else ""
                ),

                "seqkit_minus_done_pairs": (
                    r1_reads - done_pairs
                    if done_pairs is not None
                    else ""
                ),

                "flagstat_exists": (
                    "YES"
                    if flagstat_path.exists()
                    else "NO"
                ),

                "pair_status": status,

                "valid_for_summary": (
                    "YES"
                    if status == "COMPLETE_EQUAL"
                    else "NO"
                ),

                "r1_file": (
                    clean(r1_row["file_path"])
                    if r1_row is not None
                    else ""
                ),

                "r2_file": (
                    clean(r2_row["file_path"])
                    if r2_row is not None
                    else ""
                ),
            }
        )

    return (
        pd.DataFrame(sample_rows),
        skipped,
    )


def add_metadata(
    samples: pd.DataFrame,
    species_meta: pd.DataFrame,
    run_meta: pd.DataFrame,
    ddi: pd.DataFrame,
) -> pd.DataFrame:

    out = samples.copy()

    if not species_meta.empty:
        out = out.merge(
            species_meta,
            how="left",
            on="species_merged",
        )

    if not run_meta.empty:
        out = out.merge(
            run_meta,
            how="left",
            on="sample_id",
        )

    metadata_columns = [
        "metadata_species",
        "diet_species",
        "diet_species_source",
        "run_meta_species",
        "diet_run",
        "diet_run_source",
    ]

    for col in metadata_columns:

        if col not in out.columns:
            out[col] = ""

        out[col] = (
            out[col]
            .fillna("")
            .map(clean)
        )

    out["diet"] = out["diet_species"]

    use_run = (
        out["diet"].eq("")
        & out["diet_run"].ne("")
    )

    out.loc[
        use_run,
        "diet",
    ] = out.loc[
        use_run,
        "diet_run",
    ]

    out["diet_source"] = ""

    has_species_diet = (
        out["diet_species"].ne("")
    )

    out.loc[
        has_species_diet,
        "diet_source",
    ] = out.loc[
        has_species_diet,
        "diet_species_source",
    ]

    out.loc[
        use_run,
        "diet_source",
    ] = out.loc[
        use_run,
        "diet_run_source",
    ]

    for species, diet in MANUAL_DIET.items():

        mask = (
            out["species_merged"].eq(species)
            & out["diet"].eq("")
        )

        out.loc[
            mask,
            "diet",
        ] = diet

        out.loc[
            mask,
            "diet_source",
        ] = "manual_this_study"

    out["ddi_key"] = (
        out["species_merged"]
        .map(machado_key)
    )

    if not ddi.empty:
        out = out.merge(
            ddi,
            how="left",
            on="ddi_key",
        )

    for col in [
        "machado_species",
        "diet_diversity_DDI",
    ]:

        if col not in out.columns:
            out[col] = ""

        out[col] = (
            out[col]
            .fillna("")
            .map(clean)
        )

    return out


def build_species_summary(
    samples: pd.DataFrame,
) -> pd.DataFrame:

    base = (
        samples
        .groupby(
            "species_merged",
            dropna=False,
        )
        .agg(
            source_species_names=(
                "species_original",
                join_unique,
            ),
            projects=(
                "project",
                join_unique,
            ),
            project_count=(
                "project",
                "nunique",
            ),
            samples_detected=(
                "sample_id",
                "count",
            ),
            diet=(
                "diet",
                join_unique,
            ),
            diet_source=(
                "diet_source",
                join_unique,
            ),
            diet_diversity_DDI=(
                "diet_diversity_DDI",
                join_unique,
            ),
        )
        .reset_index()
    )

    valid = samples[
        samples["valid_for_summary"].eq("YES")
    ].copy()

    numeric_columns = [
        "nonhost_read_pairs",
        "nonhost_reads_total",
        "nonhost_bases_total_bp",
        "nonhost_size_Gb",
    ]

    for col in numeric_columns:
        valid[col] = pd.to_numeric(
            valid[col],
            errors="coerce",
        )

    if valid.empty:
        agg = pd.DataFrame(
            columns=["species_merged"]
        )

    else:
        agg = (
            valid
            .groupby(
                "species_merged",
                dropna=False,
            )
            .agg(
                valid_complete_samples=(
                    "sample_id",
                    "count",
                ),
                total_nonhost_read_pairs=(
                    "nonhost_read_pairs",
                    "sum",
                ),
                total_nonhost_reads=(
                    "nonhost_reads_total",
                    "sum",
                ),
                total_nonhost_bases_bp=(
                    "nonhost_bases_total_bp",
                    "sum",
                ),
                total_nonhost_size_Gb=(
                    "nonhost_size_Gb",
                    "sum",
                ),
                mean_nonhost_read_pairs_per_sample=(
                    "nonhost_read_pairs",
                    "mean",
                ),
                mean_nonhost_bases_per_sample_bp=(
                    "nonhost_bases_total_bp",
                    "mean",
                ),
                mean_nonhost_size_per_sample_Gb=(
                    "nonhost_size_Gb",
                    "mean",
                ),
                median_nonhost_size_per_sample_Gb=(
                    "nonhost_size_Gb",
                    "median",
                ),
                min_nonhost_size_per_sample_Gb=(
                    "nonhost_size_Gb",
                    "min",
                ),
                max_nonhost_size_per_sample_Gb=(
                    "nonhost_size_Gb",
                    "max",
                ),
            )
            .reset_index()
        )

    out = base.merge(
        agg,
        how="left",
        on="species_merged",
    )

    if "valid_complete_samples" not in out.columns:
        out["valid_complete_samples"] = 0

    out["valid_complete_samples"] = (
        pd.to_numeric(
            out["valid_complete_samples"],
            errors="coerce",
        )
        .fillna(0)
        .astype(int)
    )

    out["excluded_problem_samples"] = (
        out["samples_detected"]
        - out["valid_complete_samples"]
    )

    return out.sort_values(
        [
            "total_nonhost_bases_bp",
            "species_merged",
        ],
        ascending=[
            False,
            True,
        ],
        na_position="last",
    )


def choose_balanced_top(
    group: pd.DataFrame,
    top_n: int,
) -> pd.DataFrame:

    group = group.copy()

    group["project_key"] = (
        group["project"]
        .map(
            lambda x: (
                clean(x)
                or "Unknown_project"
            )
        )
    )

    queues = {}

    for project, sub in group.groupby(
        "project_key",
        dropna=False,
    ):

        sub = sub.sort_values(
            [
                "nonhost_bases_total_bp",
                "sample_id",
            ],
            ascending=[
                False,
                True,
            ],
        )

        queues[str(project)] = list(
            sub.index
        )

    selected = []

    while len(selected) < top_n:

        candidates = []

        for project, queue in queues.items():

            if not queue:
                continue

            idx = queue[0]

            candidates.append(
                (
                    project,
                    idx,
                    float(
                        group.at[
                            idx,
                            "nonhost_bases_total_bp",
                        ]
                    ),
                )
            )

        if not candidates:
            break

        candidates.sort(
            key=lambda x: (
                -x[2],
                x[0],
            )
        )

        for project, idx, _ in candidates:

            if len(selected) >= top_n:
                break

            selected.append(idx)
            queues[project].pop(0)

    if not selected:
        return group.iloc[0:0].copy()

    return group.loc[selected].copy()


def build_top5_summary(
    samples: pd.DataFrame,
    top_n: int,
) -> pd.DataFrame:

    valid = samples[
        samples["valid_for_summary"].eq("YES")
    ].copy()

    numeric_columns = [
        "nonhost_read_pairs",
        "nonhost_reads_total",
        "nonhost_bases_total_bp",
        "nonhost_size_Gb",
    ]

    for col in numeric_columns:
        valid[col] = pd.to_numeric(
            valid[col],
            errors="coerce",
        )

    rows = []

    for species, group in valid.groupby(
        "species_merged",
        dropna=False,
    ):

        selected = choose_balanced_top(
            group,
            top_n,
        )

        if selected.empty:
            continue

        project_counts = (
            selected["project"]
            .map(
                lambda x: (
                    clean(x)
                    or "Unknown_project"
                )
            )
            .value_counts()
            .sort_index()
        )

        rows.append(
            {
                "species_merged": species,

                "source_species_names": join_unique(
                    group["species_original"]
                ),

                "diet": join_unique(
                    group["diet"]
                ),

                "diet_source": join_unique(
                    group["diet_source"]
                ),

                "diet_diversity_DDI": join_unique(
                    group["diet_diversity_DDI"]
                ),

                "available_valid_samples": len(group),

                "available_projects": (
                    group["project"]
                    .nunique()
                ),

                "top_n_requested": top_n,

                "top_n_used": len(selected),

                "balanced_project_distribution": (
                    "; ".join(
                        f"{project}={number}"
                        for project, number
                        in project_counts.items()
                    )
                ),

                "top_sample_ids": (
                    "; ".join(
                        f"{row.project}:{row.sample_id}"
                        for row
                        in selected.itertuples()
                    )
                ),

                "selection_method": (
                    "project_balanced_round_robin_"
                    "then_nonhost_bases_desc"
                ),

                "rank_metric": (
                    "nonhost_bases_total_bp"
                ),

                "top_n_avg_nonhost_read_pairs": round(
                    selected[
                        "nonhost_read_pairs"
                    ].mean()
                ),

                "top_n_avg_nonhost_reads_total": round(
                    selected[
                        "nonhost_reads_total"
                    ].mean()
                ),

                "top_n_avg_nonhost_bases_bp": round(
                    selected[
                        "nonhost_bases_total_bp"
                    ].mean()
                ),

                "top_n_avg_nonhost_size_Gb": round(
                    selected[
                        "nonhost_size_Gb"
                    ].mean(),
                    6,
                ),

                "top_n_median_nonhost_size_Gb": round(
                    selected[
                        "nonhost_size_Gb"
                    ].median(),
                    6,
                ),

                "top_n_min_nonhost_size_Gb": round(
                    selected[
                        "nonhost_size_Gb"
                    ].min(),
                    6,
                ),

                "top_n_max_nonhost_size_Gb": round(
                    selected[
                        "nonhost_size_Gb"
                    ].max(),
                    6,
                ),
            }
        )

    if not rows:
        return pd.DataFrame()

    return (
        pd.DataFrame(rows)
        .sort_values(
            [
                "top_n_avg_nonhost_bases_bp",
                "species_merged",
            ],
            ascending=[
                False,
                True,
            ],
        )
    )


def optional_path(
    value: str,
) -> Path | None:

    if not clean(value):
        return None

    return Path(
        value
    ).expanduser().resolve()


def main() -> None:

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--reads-root",
        required=True,
    )

    parser.add_argument(
        "--seqkit-stats",
        required=True,
    )

    parser.add_argument(
        "--species-meta",
        default="",
    )

    parser.add_argument(
        "--run-meta",
        default="",
    )

    parser.add_argument(
        "--machado",
        default="",
    )

    parser.add_argument(
        "--outdir",
        required=True,
    )

    parser.add_argument(
        "--top-n",
        type=int,
        default=5,
    )

    parser.add_argument(
        "--include-without-done",
        action="store_true",
    )

    args = parser.parse_args()

    reads_root = Path(
        args.reads_root
    ).expanduser().resolve()

    stats_path = Path(
        args.seqkit_stats
    ).expanduser().resolve()

    outdir = Path(
        args.outdir
    ).expanduser().resolve()

    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not reads_root.is_dir():
        raise SystemExit(
            f"ERROR: reads root not found: {reads_root}"
        )

    if not stats_path.is_file():
        raise SystemExit(
            f"ERROR: seqkit stats not found: {stats_path}"
        )

    samples, skipped = build_samples(
        stats_path,
        reads_root,
        not args.include_without_done,
    )

    samples = add_metadata(
        samples,
        load_species_meta(
            optional_path(args.species_meta)
        ),
        load_run_meta(
            optional_path(args.run_meta)
        ),
        load_ddi(
            optional_path(args.machado)
        ),
    )

    samples = samples.sort_values(
        [
            "species_merged",
            "species_original",
            "project",
            "sample_id",
        ]
    )

    samples.to_csv(
        outdir / "01_nonhost_sample_summary.tsv",
        sep="\t",
        index=False,
    )

    build_species_summary(
        samples
    ).to_csv(
        outdir / "02_nonhost_species_summary.tsv",
        sep="\t",
        index=False,
    )

    build_top5_summary(
        samples,
        args.top_n,
    ).to_csv(
        outdir / "03_nonhost_species_top5_average.tsv",
        sep="\t",
        index=False,
    )

    samples[
        samples["valid_for_summary"].ne("YES")
    ].to_csv(
        outdir / "04_nonhost_problem_samples.tsv",
        sep="\t",
        index=False,
    )

    samples[
        samples["metadata_species"]
        .fillna("")
        .eq("")
    ][
        [
            "species_original",
            "species_merged",
        ]
    ].drop_duplicates().to_csv(
        outdir / "05_unmatched_species_metadata.tsv",
        sep="\t",
        index=False,
    )

    samples[
        samples["diet_diversity_DDI"]
        .fillna("")
        .eq("")
    ][
        [
            "species_original",
            "species_merged",
            "diet",
        ]
    ].drop_duplicates().to_csv(
        outdir / "06_unmatched_species_ddi.tsv",
        sep="\t",
        index=False,
    )

    pd.DataFrame(
        {
            "skipped_file": skipped
        }
    ).to_csv(
        outdir / "07_skipped_nonhost_fastq_names.tsv",
        sep="\t",
        index=False,
    )

    print(
        f"Samples detected: {len(samples)}"
    )

    print(
        "Valid samples: "
        f"{(samples['valid_for_summary'] == 'YES').sum()}"
    )

    print(
        "Problem samples: "
        f"{(samples['valid_for_summary'] != 'YES').sum()}"
    )

    print(
        f"Results: {outdir}"
    )


if __name__ == "__main__":
    main()
