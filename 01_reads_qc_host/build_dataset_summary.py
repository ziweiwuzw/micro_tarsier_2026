#!/usr/bin/env python3
import re
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


BAD_VALUES = {
    "",
    "nan",
    "none",
    "na",
    "n/a",
    "not",
    "not provided",
    "missing",
    "missing: not provided",
    "missing: not collected",
    "not collected",
    "-",
}


def clean_value(x):
    if pd.isna(x):
        return ""
    x = str(x).strip()
    if x.lower() in BAD_VALUES:
        return ""
    return x


def norm_species_name(x):
    x = clean_value(x)
    x = x.replace(" ", "_")
    x = re.sub(r"[^A-Za-z0-9_]+", "_", x)
    x = re.sub(r"_+", "_", x).strip("_")
    return x


def canonical_binomial_species(x):
    x = norm_species_name(x)

    synonym_map = {
        "Xanthonycticebus_pygmaeus": "Nycticebus_pygmaeus",
        "Cephalophacus_bancanus": "Cephalopachus_bancanus",
    }

    x = synonym_map.get(x, x)

    parts = x.split("_")
    if len(parts) >= 2:
        return "_".join(parts[:2])

    return x


def normalize_species_to_machado_key(x):
    x = clean_value(x)
    x = x.replace("_", " ")
    x = re.sub(r"\s+", " ", x).strip().lower()

    synonym_map = {
        "xanthonycticebus pygmaeus": "nycticebus pygmaeus",
        "cephalophacus bancanus": "cephalopachus bancanus",
        "cephalopachus bancanus": "cephalopachus bancanus",
        "lagothrix lagotricha": "lagothrix lagothricha",
        "rhinopithecus roxellanae": "rhinopithecus roxellana",
    }

    if x in synonym_map:
        return synonym_map[x]

    parts = x.split()
    if len(parts) >= 2:
        return " ".join(parts[:2])

    return x


def bytes_to_gb(x):
    return round(float(x) / 1024 / 1024 / 1024, 3)


def split_and_clean_semicolon_values(x):
    x = clean_value(x)
    if not x:
        return []

    parts = [p.strip() for p in str(x).split(";")]
    out = []

    for p in parts:
        p = clean_value(p)
        if p and p not in out:
            out.append(p)

    return out


def join_unique(values):
    out = []
    for v in values:
        for token in split_and_clean_semicolon_values(v):
            if token and token not in out:
                out.append(token)
    return "; ".join(out)


def first_nonempty(values):
    for v in values:
        v = clean_value(v)
        if v:
            return v
    return ""


def to_num(s):
    return pd.to_numeric(s, errors="coerce").fillna(0)


def format_int_or_blank(x):
    if pd.isna(x):
        return ""
    try:
        x = float(x)
    except Exception:
        return ""
    if not np.isfinite(x) or x <= 0:
        return ""
    return int(round(x, 0))


def format_gb_or_blank(x, ndigits=3):
    if pd.isna(x):
        return ""
    try:
        x = float(x)
    except Exception:
        return ""
    if not np.isfinite(x) or x <= 0:
        return ""
    return round(x, ndigits)


def parse_run_and_mate(filename):
    name = filename

    for suffix in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break

    name = re.sub(r"\.clean$", "", name)
    name = re.sub(r"_clean$", "", name)

    run = name
    mate = "SE"

    patterns = [
        (r"(.+)_R1_001$", "1"),
        (r"(.+)_R2_001$", "2"),
        (r"(.+)_1$", "1"),
        (r"(.+)_2$", "2"),
        (r"(.+)_R1$", "1"),
        (r"(.+)_R2$", "2"),
        (r"(.+)\.1$", "1"),
        (r"(.+)\.2$", "2"),
    ]

    for pat, mlabel in patterns:
        m = re.match(pat, name)
        if m:
            run = m.group(1)
            mate = mlabel
            break

    return run, mate


def scan_fastq_files(base_dir):
    base_dir = Path(base_dir).resolve()

    roots = [
        base_dir / "kingfisher_download_by_host",
        base_dir / "special_downloads",
        base_dir / "01_clean_data" / "01_clean_reads",
    ]

    rows = []

    for root in roots:
        if not root.exists():
            continue

        for p in root.rglob("*"):
            if not p.is_file():
                continue

            if not re.search(r"\.(fastq|fq)(\.gz)?$", p.name, re.IGNORECASE):
                continue

            try:
                rel = p.relative_to(base_dir)
            except ValueError:
                continue

            parts = rel.parts

            source_class = ""
            source_detail = ""
            source_species = ""
            project = ""

            if len(parts) >= 4 and parts[0] == "kingfisher_download_by_host":
                if parts[1].startswith("_kingfisher_"):
                    continue
                source_class = "raw"
                source_detail = "raw_kingfisher"
                source_species = parts[1]
                project = parts[2]

            elif len(parts) >= 4 and parts[0] == "special_downloads":
                source_class = "raw"
                source_detail = "raw_special"
                project = parts[1]
                source_species = parts[2]

            elif len(parts) >= 5 and parts[0] == "01_clean_data" and parts[1] == "01_clean_reads":
                source_class = "clean"
                source_detail = "clean_fastp"
                source_species = parts[2]
                project = parts[3]

            else:
                continue

            run, mate = parse_run_and_mate(p.name)

            rows.append({
                "source_species": source_species,
                "source_species_key": norm_species_name(source_species),
                "binomial_name": canonical_binomial_species(source_species),
                "project": project,
                "run": run,
                "mate": mate,
                "source_class": source_class,
                "source_detail": source_detail,
                "file_path": str(p),
                "file_name": p.name,
                "size_bytes": p.stat().st_size,
                "size_GB": bytes_to_gb(p.stat().st_size),
            })

    return pd.DataFrame(rows)


def normalize_columns(df):
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]
    return df


def pick_column(df, candidates):
    norm_map = {}
    for c in df.columns:
        key = re.sub(r"[^a-z0-9]+", "", c.lower())
        norm_map[key] = c

    for cand in candidates:
        key = re.sub(r"[^a-z0-9]+", "", cand.lower())
        if key in norm_map:
            return norm_map[key]

    return None


def load_run_metadata(path):
    path = Path(path)
    if not path.exists():
        print(f"WARNING: run metadata not found: {path}")
        return pd.DataFrame()

    df = pd.read_csv(path, dtype=str)
    df = normalize_columns(df)

    run_col = pick_column(df, ["run", "Run", "run_accession", "Run accession", "accession"])
    if run_col is None:
        print(f"WARNING: run metadata has no run column: {path}")
        return pd.DataFrame()

    out = pd.DataFrame()
    out["run"] = df[run_col].map(clean_value)

    col_map = {
        "meta_bioproject": ["bioproject", "BioProject", "study_accession", "Study accession", "project"],
        "meta_species": ["taxon_name", "Taxon name", "host scientific name", "Host scientific name", "scientific_name"],
        "meta_model": ["model", "Model", "instrument_model", "instrument model", "sequencing platform", "platform"],
        "meta_country": ["country", "Country"],
        "meta_geo_loc": ["geo loc name", "geo_loc_name", "geographic location", "geographic_location"],
        "meta_state": ["state", "State", "province", "region", "location"],
        "meta_habitat": ["habitat2", "habitat", "source_name"],
        "meta_diet": ["diet", "Diet"],
        "meta_biosample": ["biosample", "BioSample", "sample_accession"],
        "meta_sample_name": ["sample_name", "Sample name", "sample alias"],
        "meta_library_layout": ["library_layout", "Library layout"],
        "meta_bases": ["bases", "Bases"],
        "meta_spots": ["spots", "Spots"],
    }

    for new_col, candidates in col_map.items():
        old_col = pick_column(df, candidates)
        if old_col is not None:
            out[new_col] = df[old_col].map(clean_value)
        else:
            out[new_col] = ""

    out = out.drop_duplicates("run")
    return out


def load_species_metadata(path):
    path = Path(path)
    if not path.exists():
        print(f"WARNING: species metadata not found: {path}")
        return pd.DataFrame()

    df = pd.read_csv(path, dtype=str)
    df = normalize_columns(df)

    species_col = pick_column(df, ["Binomial name", "Species", "species", "Scientific name", "Primate species"])
    if species_col is None:
        print(f"WARNING: species metadata has no species column: {path}")
        return pd.DataFrame()

    out = pd.DataFrame()
    out["binomial_name"] = df[species_col].map(canonical_binomial_species)

    col_map = {
        "species_meta_country": ["Country", "country"],
        "species_meta_continent": ["Continent", "continent"],
        "species_meta_diet": ["Diet", "diet"],
        "species_meta_living_environment": ["Captive or Wild", "Living environment", "environment"],
        "species_meta_article": ["Article", "Publication", "Document"],
        "species_meta_project": ["Bioproject/DDBJ", "BioProject", "Project"],
        "species_meta_data_type": ["Data type", "data type"],
    }

    for new_col, candidates in col_map.items():
        old_col = pick_column(df, candidates)
        if old_col is not None:
            out[new_col] = df[old_col].map(clean_value)
        else:
            out[new_col] = ""

    out = (
        out.groupby("binomial_name", dropna=False)
        .agg({
            "species_meta_country": join_unique,
            "species_meta_continent": join_unique,
            "species_meta_diet": join_unique,
            "species_meta_living_environment": join_unique,
            "species_meta_article": join_unique,
            "species_meta_project": join_unique,
            "species_meta_data_type": join_unique,
        })
        .reset_index()
    )

    return out


def format_project_with_source_species(df):
    items = []

    for project, g in df.groupby("project", dropna=False):
        project = clean_value(project)
        source_species = sorted(set(g["source_species"].map(clean_value)))
        source_species = [x for x in source_species if x]

        if not project:
            continue

        if source_species:
            items.append(f"{project}({';'.join(source_species)})")
        else:
            items.append(project)

    return "; ".join(items)



def load_allowed_sample_keys_from_sample_fastp(path, allowed_layout=""):
    """
    Build allowed PE/SE sample keys from sample_fastp_summary.tsv.

    Key format:
    canonical_binomial_species + project + sample_id

    This is used to filter the whole final table, including:
    Sample number, disk size, reported size, and Top5 sample selection.
    """
    if not allowed_layout:
        return None

    path = Path(path)

    if not path.exists():
        raise SystemExit(f"ERROR: strict layout filtering was requested, but sample table was not found: {path}")

    df = pd.read_csv(path, sep="\t", dtype=str)

    if len(df.columns) <= 3:
        df = pd.read_csv(path, sep=r"\s+", dtype=str, engine="python")

    df.columns = [c.strip() for c in df.columns]

    required = ["project", "sample_id", "layout"]
    for c in required:
        if c not in df.columns:
            raise SystemExit(f"ERROR: strict layout filtering requires column '{c}' in {path}")

    if "species_merged" in df.columns:
        species_col = "species_merged"
    elif "species_original" in df.columns:
        species_col = "species_original"
    else:
        raise SystemExit(f"ERROR: strict layout filtering requires species_merged or species_original in {path}")

    allowed_set = {x.strip().upper() for x in allowed_layout.split(",") if x.strip()}

    before = len(df)
    df["layout"] = df["layout"].astype(str).str.upper()
    df = df[df["layout"].isin(allowed_set)].copy()
    after = len(df)

    df["binomial_key"] = df[species_col].map(canonical_binomial_species)
    df["project"] = df["project"].map(clean_value)
    df["sample_id"] = df["sample_id"].map(clean_value)

    df = df[
        df["binomial_key"].map(clean_value).ne("") &
        df["project"].map(clean_value).ne("") &
        df["sample_id"].map(clean_value).ne("")
    ].copy()

    keys = set(
        df["binomial_key"].astype(str) + "\t" +
        df["project"].astype(str) + "\t" +
        df["sample_id"].astype(str)
    )

    print(f"Strict layout filter from {path}: {before} rows -> {after} rows; allowed sample keys = {len(keys)}")

    return keys


def filter_inventory_by_allowed_sample_keys(inv, allowed_keys, label="inventory"):
    if allowed_keys is None:
        return inv

    before_files = len(inv)
    before_samples = inv[["binomial_name", "project", "run"]].drop_duplicates().shape[0]

    inv = inv.copy()
    inv["_allowed_key"] = (
        inv["binomial_name"].map(canonical_binomial_species).astype(str) + "\t" +
        inv["project"].map(clean_value).astype(str) + "\t" +
        inv["run"].map(clean_value).astype(str)
    )

    inv = inv[inv["_allowed_key"].isin(allowed_keys)].copy()
    inv = inv.drop(columns=["_allowed_key"], errors="ignore")

    after_files = len(inv)
    after_samples = inv[["binomial_name", "project", "run"]].drop_duplicates().shape[0] if not inv.empty else 0

    print(f"Strict layout filter applied to {label}: files {before_files} -> {after_files}; samples {before_samples} -> {after_samples}")

    return inv

def build_base_species_summary(base_dir, meta_run, meta_species, outdir, sample_fastp_summary_path=None, allowed_layout="", write_debug=False):
    base_dir = Path(base_dir).resolve()
    outdir = Path(outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    inv = scan_fastq_files(base_dir)

    if inv.empty:
        raise SystemExit("ERROR: 没有扫描到 fastq/fq 文件，请检查路径。")

    allowed_keys = load_allowed_sample_keys_from_sample_fastp(
        sample_fastp_summary_path,
        allowed_layout=allowed_layout
    )

    inv = filter_inventory_by_allowed_sample_keys(
        inv,
        allowed_keys,
        label="main FASTQ inventory"
    )

    if inv.empty:
        raise SystemExit("ERROR: strict layout filtering removed all FASTQ files. Please check sample_fastp_summary.tsv sample_id/project names.")

    source_species_with_raw = set(
        inv.loc[inv["source_class"] == "raw", "source_species_key"].unique()
    )

    inv["use_for_report"] = False
    inv["size_source_decision"] = ""

    inv.loc[
        inv["source_species_key"].isin(source_species_with_raw) &
        (inv["source_class"] == "raw"),
        "use_for_report"
    ] = True

    inv.loc[
        inv["source_species_key"].isin(source_species_with_raw),
        "size_source_decision"
    ] = "raw_exists_for_source_species_use_raw_only"

    inv.loc[
        (~inv["source_species_key"].isin(source_species_with_raw)) &
        (inv["source_class"] == "clean"),
        "use_for_report"
    ] = True

    inv.loc[
        ~inv["source_species_key"].isin(source_species_with_raw),
        "size_source_decision"
    ] = "raw_absent_for_source_species_use_clean"

    selected = inv[inv["use_for_report"]].copy()

    if write_debug:
        inv.to_csv(outdir / "debug_00_all_fastq_file_inventory_raw_and_clean.csv", index=False)
        selected.to_csv(outdir / "debug_01_selected_fastq_files_for_report.csv", index=False)

    presence = (
        inv.groupby(["binomial_name", "source_species", "project", "run", "source_class"], dropna=False)
        .agg(
            file_number=("file_path", "count"),
            size_GB=("size_GB", "sum"),
            mates=("mate", lambda x: ",".join(sorted(set(map(str, x)))))
        )
        .reset_index()
    )

    raw_status = (
        presence[presence["source_class"] == "raw"]
        .groupby("binomial_name", dropna=False)
        .agg(
            Raw_sample_number_on_disk=("run", "nunique"),
            Raw_fastq_size_GB_on_disk=("size_GB", "sum"),
        )
        .reset_index()
    )

    clean_status = (
        presence[presence["source_class"] == "clean"]
        .groupby("binomial_name", dropna=False)
        .agg(
            Clean_sample_number_on_disk=("run", "nunique"),
            Clean_fastq_size_GB_on_disk=("size_GB", "sum"),
        )
        .reset_index()
    )

    sample_table = (
        selected.groupby(["binomial_name", "source_species", "project", "run"], dropna=False)
        .agg(
            Fastq_file_number=("file_path", "count"),
            Current_representative_fastq_size_GB=("size_GB", "sum"),
            Source_used=("source_class", join_unique),
            Source_detail=("source_detail", join_unique),
            Size_source=("size_source_decision", join_unique),
            Mates_detected=("mate", lambda x: ",".join(sorted(set(map(str, x)))))
        )
        .reset_index()
    )

    sample_table["Current_representative_fastq_size_GB"] = (
        sample_table["Current_representative_fastq_size_GB"].round(3)
    )

    run_meta = load_run_metadata(meta_run)
    species_meta = load_species_metadata(meta_species)

    if not run_meta.empty:
        sample_table = sample_table.merge(run_meta, how="left", on="run")

    if not species_meta.empty:
        sample_table = sample_table.merge(species_meta, how="left", on="binomial_name")

    sample_table["Country"] = sample_table.apply(
        lambda r: first_nonempty([
            r.get("meta_country", ""),
            r.get("species_meta_country", "")
        ]),
        axis=1
    )

    sample_table["State/Location"] = sample_table.apply(
        lambda r: join_unique([
            r.get("meta_state", ""),
            r.get("meta_geo_loc", ""),
            r.get("meta_habitat", "")
        ]),
        axis=1
    )

    sample_table["Continent"] = sample_table.apply(
        lambda r: first_nonempty([
            r.get("species_meta_continent", "")
        ]),
        axis=1
    )

    sample_table["Diet"] = sample_table.apply(
        lambda r: first_nonempty([
            r.get("species_meta_diet", ""),
            r.get("meta_diet", "")
        ]),
        axis=1
    )

    sample_table["Diet source"] = sample_table.apply(
        lambda r: "ENA_gut_mirco_source_Primates.csv"
        if clean_value(r.get("species_meta_diet", ""))
        else ("run_metadata_after_removing_Not" if clean_value(r.get("meta_diet", "")) else ""),
        axis=1
    )

    sample_table["Living environment"] = sample_table.apply(
        lambda r: first_nonempty([
            r.get("species_meta_living_environment", "")
        ]),
        axis=1
    )

    sample_table["Article"] = sample_table.apply(
        lambda r: first_nonempty([
            r.get("species_meta_article", "")
        ]),
        axis=1
    )

    sample_table["Sequencing platform/model"] = sample_table.apply(
        lambda r: first_nonempty([
            r.get("meta_model", "")
        ]),
        axis=1
    )

    if write_debug:
        sample_table.to_csv(outdir / "debug_02_sample_level_after_metadata_merge.csv", index=False)

    project_annot_rows = []

    for binomial_name, g in selected.groupby("binomial_name", dropna=False):
        project_annot_rows.append({
            "binomial_name": binomial_name,
            "Bioproject/DDBJ": format_project_with_source_species(g)
        })

    project_annot = pd.DataFrame(project_annot_rows)

    species_summary = (
        sample_table.groupby("binomial_name", dropna=False)
        .agg(
            Project_number=("project", "nunique"),
            Source_species_names=("source_species", join_unique),
            Sample_number=("run", "nunique"),
            Reported_overall_data_size_GB=("Current_representative_fastq_size_GB", "sum"),
            Size_source_for_reported_overall_data=("Size_source", join_unique),
            Source_used=("Source_used", join_unique),
            Country=("Country", join_unique),
            State_Location=("State/Location", join_unique),
            Continent=("Continent", join_unique),
            Diet=("Diet", join_unique),
            Diet_source=("Diet source", join_unique),
            Living_environment=("Living environment", join_unique),
            Article=("Article", join_unique),
            Sequencing_platform_model=("Sequencing platform/model", join_unique),
        )
        .reset_index()
    )

    species_summary = species_summary.merge(project_annot, how="left", on="binomial_name")
    species_summary = species_summary.merge(raw_status, how="left", on="binomial_name")
    species_summary = species_summary.merge(clean_status, how="left", on="binomial_name")

    for c in [
        "Raw_sample_number_on_disk",
        "Clean_sample_number_on_disk",
        "Raw_fastq_size_GB_on_disk",
        "Clean_fastq_size_GB_on_disk",
    ]:
        if c in species_summary.columns:
            species_summary[c] = species_summary[c].fillna(0)

    species_summary = species_summary.rename(columns={
        "binomial_name": "Binomial name",
        "Project_number": "Project number",
        "Source_species_names": "Source species names",
        "Sample_number": "Sample number",
        "Reported_overall_data_size_GB": "Reported overall data size (GB)",
        "Size_source_for_reported_overall_data": "Size source for reported overall data",
        "State_Location": "State/Location",
        "Living_environment": "Living environment",
        "Diet_source": "Diet source",
        "Sequencing_platform_model": "Sequencing platform/model",
        "Raw_sample_number_on_disk": "Raw sample number on disk",
        "Clean_sample_number_on_disk": "Clean sample number on disk",
        "Raw_fastq_size_GB_on_disk": "Raw fastq size on disk (GB)",
        "Clean_fastq_size_GB_on_disk": "Clean fastq size on disk (GB)",
    })

    species_summary["Reported overall data size (GB)"] = (
        species_summary["Reported overall data size (GB)"].round(3)
    )
    species_summary["Raw fastq size on disk (GB)"] = (
        species_summary["Raw fastq size on disk (GB)"].round(3)
    )
    species_summary["Clean fastq size on disk (GB)"] = (
        species_summary["Clean fastq size on disk (GB)"].round(3)
    )

    mask_ceph = species_summary["Binomial name"].eq("Cephalopachus_bancanus")
    if mask_ceph.any():
        species_summary.loc[mask_ceph, "Diet"] = "Insectivore"
        species_summary.loc[mask_ceph, "Diet source"] = "manual_this_study"
        species_summary.loc[mask_ceph, "Country"] = "Indonesia"
        species_summary.loc[mask_ceph, "State/Location"] = "Sriwijaya University, Indonesia"
        species_summary.loc[mask_ceph, "Continent"] = "Asia"
        species_summary.loc[mask_ceph, "Living environment"] = "Wild"
        species_summary.loc[mask_ceph, "Article"] = "This study"
        species_summary.loc[mask_ceph, "Sequencing platform/model"] = "Illumina NovaSeqX"

    return species_summary


def read_fastp_species_summary(path):
    path = Path(path)

    if not path.exists():
        print(f"WARNING: fastp species summary not found: {path}")
        return pd.DataFrame(columns=["binomial_key"])

    df = pd.read_csv(path, sep="\t", dtype=str)

    if len(df.columns) <= 3:
        df = pd.read_csv(path, sep=r"\s+", dtype=str, engine="python")

    df.columns = [c.strip() for c in df.columns]

    if "species" not in df.columns:
        print("WARNING: species_summary_merged_level.tsv 中没有 species 列，跳过 fastp reads/bases 添加。")
        return pd.DataFrame(columns=["binomial_key"])

    df["binomial_key"] = df["species"].map(canonical_binomial_species)

    needed_cols = [
        "samples",
        "missing_fastp_json_samples",
        "raw_reads_total",
        "raw_bases_total",
        "clean_reads_total",
        "clean_bases_total",
        "raw_read_pairs_or_se_reads",
        "clean_read_pairs_or_se_reads",
    ]

    for c in needed_cols:
        if c not in df.columns:
            df[c] = 0
        df[c] = to_num(df[c])

    ss = (
        df.groupby("binomial_key", dropna=False)
        .agg(
            fastp_summary_samples=("samples", "sum"),
            missing_fastp_json_samples=("missing_fastp_json_samples", "sum"),
            raw_reads_total=("raw_reads_total", "sum"),
            raw_bases_total=("raw_bases_total", "sum"),
            clean_reads_total=("clean_reads_total", "sum"),
            clean_bases_total=("clean_bases_total", "sum"),
            raw_read_pairs_or_se_reads=("raw_read_pairs_or_se_reads", "sum"),
            clean_read_pairs_or_se_reads=("clean_read_pairs_or_se_reads", "sum"),
        )
        .reset_index()
    )

    ss["fastp_json_available_samples"] = (
        ss["fastp_summary_samples"] - ss["missing_fastp_json_samples"]
    ).clip(lower=0)

    return ss


def safe_median_ratio(numerator, denominator):
    numerator = pd.to_numeric(numerator, errors="coerce")
    denominator = pd.to_numeric(denominator, errors="coerce")

    ratio = numerator / denominator
    ratio = ratio.replace([np.inf, -np.inf], np.nan).dropna()
    ratio = ratio[ratio > 0]

    if len(ratio) == 0:
        return np.nan

    return ratio.median()


def estimate_one_value(row, raw_disk_col, clean_disk_col, reported_col,
                       ratio_rawdisk, ratio_cleandisk, ratio_reported,
                       fallback_label_prefix):
    raw_disk = row.get(raw_disk_col, 0)
    clean_disk = row.get(clean_disk_col, 0)
    reported = row.get(reported_col, 0)

    if raw_disk > 0 and pd.notna(ratio_rawdisk) and ratio_rawdisk > 0:
        return raw_disk * ratio_rawdisk, f"{fallback_label_prefix}_from_raw_gz_disk_median_ratio"

    if clean_disk > 0 and pd.notna(ratio_cleandisk) and ratio_cleandisk > 0:
        return clean_disk * ratio_cleandisk, f"{fallback_label_prefix}_from_clean_gz_disk_median_ratio"

    if reported > 0 and pd.notna(ratio_reported) and ratio_reported > 0:
        return reported * ratio_reported, f"{fallback_label_prefix}_from_reported_gz_size_median_ratio"

    return np.nan, "fallback_failed_no_usable_disk_size_or_ratio"


def add_fastp_reads_bases(species_summary, fastp_summary_path):
    df = species_summary.copy()
    df["binomial_key"] = df["Binomial name"].map(canonical_binomial_species)

    ss = read_fastp_species_summary(fastp_summary_path)

    if ss.empty or "binomial_key" not in ss.columns:
        for c in [
            "Raw reads total",
            "Raw read pairs or SE reads",
            "Raw bases total (bp)",
            "Raw sequencing size (Gb bases)",
            "Clean reads total",
            "Clean read pairs or SE reads",
            "Clean bases total (bp)",
            "Clean sequencing size (Gb bases)",
            "Fastp / size estimate source",
        ]:
            df[c] = ""
        return df.drop(columns=["binomial_key"], errors="ignore")

    merged = df.merge(ss, how="left", on="binomial_key")

    for c in [
        "Raw fastq size on disk (GB)",
        "Clean fastq size on disk (GB)",
        "Reported overall data size (GB)",
        "fastp_summary_samples",
        "missing_fastp_json_samples",
        "fastp_json_available_samples",
        "raw_reads_total",
        "raw_bases_total",
        "clean_reads_total",
        "clean_bases_total",
        "raw_read_pairs_or_se_reads",
        "clean_read_pairs_or_se_reads",
    ]:
        if c not in merged.columns:
            merged[c] = 0
        merged[c] = to_num(merged[c])

    total_samples = merged["fastp_summary_samples"]
    available_samples = merged["fastp_json_available_samples"]
    missing_samples = merged["missing_fastp_json_samples"]

    has_fastp_data = (
        (merged["raw_bases_total"] > 0) |
        (merged["clean_bases_total"] > 0) |
        (merged["raw_reads_total"] > 0) |
        (merged["clean_reads_total"] > 0)
    )

    extrapolate_factor = np.where(
        (has_fastp_data) &
        (available_samples > 0) &
        (total_samples > available_samples),
        total_samples / available_samples,
        1.0
    )

    merged["Raw reads final"] = np.where(
        merged["raw_reads_total"] > 0,
        merged["raw_reads_total"] * extrapolate_factor,
        np.nan
    )

    merged["Raw read pairs or SE reads final"] = np.where(
        merged["raw_read_pairs_or_se_reads"] > 0,
        merged["raw_read_pairs_or_se_reads"] * extrapolate_factor,
        np.nan
    )

    merged["Raw bases total final (bp)"] = np.where(
        merged["raw_bases_total"] > 0,
        merged["raw_bases_total"] * extrapolate_factor,
        np.nan
    )

    merged["Clean reads final"] = np.where(
        merged["clean_reads_total"] > 0,
        merged["clean_reads_total"] * extrapolate_factor,
        np.nan
    )

    merged["Clean read pairs or SE reads final"] = np.where(
        merged["clean_read_pairs_or_se_reads"] > 0,
        merged["clean_read_pairs_or_se_reads"] * extrapolate_factor,
        np.nan
    )

    merged["Clean bases total final (bp)"] = np.where(
        merged["clean_bases_total"] > 0,
        merged["clean_bases_total"] * extrapolate_factor,
        np.nan
    )

    merged["Fastp estimate source"] = np.where(
        has_fastp_data & (missing_samples == 0),
        "fastp_summary_direct",
        np.where(
            has_fastp_data & (missing_samples > 0) & (available_samples > 0),
            "fastp_summary_extrapolated_from_completed_samples",
            "needs_fallback_from_other_species"
        )
    )

    training = merged[has_fastp_data].copy()

    raw_disk_col = "Raw fastq size on disk (GB)"
    clean_disk_col = "Clean fastq size on disk (GB)"
    reported_col = "Reported overall data size (GB)"

    ratio_map = {}

    targets = [
        ("Raw reads final", "raw_reads"),
        ("Raw read pairs or SE reads final", "raw_pairs_or_se"),
        ("Raw bases total final (bp)", "raw_bases"),
        ("Clean reads final", "clean_reads"),
        ("Clean read pairs or SE reads final", "clean_pairs_or_se"),
        ("Clean bases total final (bp)", "clean_bases"),
    ]

    for target_col, short_name in targets:
        ratio_map[(short_name, "rawdisk")] = safe_median_ratio(
            training[target_col],
            training[raw_disk_col]
        )
        ratio_map[(short_name, "cleandisk")] = safe_median_ratio(
            training[target_col],
            training[clean_disk_col]
        )
        ratio_map[(short_name, "reported")] = safe_median_ratio(
            training[target_col],
            training[reported_col]
        )

    fallback_rows = merged["Fastp estimate source"].eq("needs_fallback_from_other_species")

    for idx in merged.index[fallback_rows]:
        row = merged.loc[idx]

        raw_reads, src1 = estimate_one_value(
            row, raw_disk_col, clean_disk_col, reported_col,
            ratio_map[("raw_reads", "rawdisk")],
            ratio_map[("raw_reads", "cleandisk")],
            ratio_map[("raw_reads", "reported")],
            "raw_reads_estimated"
        )

        raw_pairs, src2 = estimate_one_value(
            row, raw_disk_col, clean_disk_col, reported_col,
            ratio_map[("raw_pairs_or_se", "rawdisk")],
            ratio_map[("raw_pairs_or_se", "cleandisk")],
            ratio_map[("raw_pairs_or_se", "reported")],
            "raw_pairs_or_se_estimated"
        )

        raw_bases, src3 = estimate_one_value(
            row, raw_disk_col, clean_disk_col, reported_col,
            ratio_map[("raw_bases", "rawdisk")],
            ratio_map[("raw_bases", "cleandisk")],
            ratio_map[("raw_bases", "reported")],
            "raw_bases_estimated"
        )

        clean_reads, src4 = estimate_one_value(
            row, raw_disk_col, clean_disk_col, reported_col,
            ratio_map[("clean_reads", "rawdisk")],
            ratio_map[("clean_reads", "cleandisk")],
            ratio_map[("clean_reads", "reported")],
            "clean_reads_estimated"
        )

        clean_pairs, src5 = estimate_one_value(
            row, raw_disk_col, clean_disk_col, reported_col,
            ratio_map[("clean_pairs_or_se", "rawdisk")],
            ratio_map[("clean_pairs_or_se", "cleandisk")],
            ratio_map[("clean_pairs_or_se", "reported")],
            "clean_pairs_or_se_estimated"
        )

        clean_bases, src6 = estimate_one_value(
            row, raw_disk_col, clean_disk_col, reported_col,
            ratio_map[("clean_bases", "rawdisk")],
            ratio_map[("clean_bases", "cleandisk")],
            ratio_map[("clean_bases", "reported")],
            "clean_bases_estimated"
        )

        merged.at[idx, "Raw reads final"] = raw_reads
        merged.at[idx, "Raw read pairs or SE reads final"] = raw_pairs
        merged.at[idx, "Raw bases total final (bp)"] = raw_bases

        merged.at[idx, "Clean reads final"] = clean_reads
        merged.at[idx, "Clean read pairs or SE reads final"] = clean_pairs
        merged.at[idx, "Clean bases total final (bp)"] = clean_bases

        merged.at[idx, "Fastp estimate source"] = "; ".join(sorted(set([
            src1, src2, src3, src4, src5, src6
        ])))

    merged["Raw sequencing size (Gb bases)"] = (
        merged["Raw bases total final (bp)"] / 1e9
    )

    merged["Clean sequencing size (Gb bases)"] = (
        merged["Clean bases total final (bp)"] / 1e9
    )

    merged["Raw reads total"] = merged["Raw reads final"].map(format_int_or_blank)
    merged["Raw read pairs or SE reads"] = merged["Raw read pairs or SE reads final"].map(format_int_or_blank)
    merged["Raw bases total (bp)"] = merged["Raw bases total final (bp)"].map(format_int_or_blank)
    merged["Raw sequencing size (Gb bases)"] = merged["Raw sequencing size (Gb bases)"].map(format_gb_or_blank)

    merged["Clean reads total"] = merged["Clean reads final"].map(format_int_or_blank)
    merged["Clean read pairs or SE reads"] = merged["Clean read pairs or SE reads final"].map(format_int_or_blank)
    merged["Clean bases total (bp)"] = merged["Clean bases total final (bp)"].map(format_int_or_blank)
    merged["Clean sequencing size (Gb bases)"] = merged["Clean sequencing size (Gb bases)"].map(format_gb_or_blank)

    merged["Fastp / size estimate source"] = merged["Fastp estimate source"]

    drop_cols = [
        "binomial_key",
        "fastp_summary_samples",
        "missing_fastp_json_samples",
        "fastp_json_available_samples",
        "raw_reads_total",
        "raw_bases_total",
        "clean_reads_total",
        "clean_bases_total",
        "raw_read_pairs_or_se_reads",
        "clean_read_pairs_or_se_reads",
        "Raw reads final",
        "Raw read pairs or SE reads final",
        "Raw bases total final (bp)",
        "Clean reads final",
        "Clean read pairs or SE reads final",
        "Clean bases total final (bp)",
        "Fastp estimate source",
    ]

    merged = merged.drop(columns=[c for c in drop_cols if c in merged.columns])

    return merged


def read_sample_fastp_summary(path, allowed_layout=""):
    path = Path(path)

    if not path.exists():
        print(f"WARNING: sample_fastp_summary.tsv not found: {path}")
        return pd.DataFrame()

    df = pd.read_csv(path, sep="\t", dtype=str)

    if len(df.columns) <= 3:
        df = pd.read_csv(path, sep=r"\s+", dtype=str, engine="python")

    df.columns = [c.strip() for c in df.columns]

    if allowed_layout:
        if "layout" not in df.columns:
            raise SystemExit("ERROR: --strict-layout was set, but sample_fastp_summary.tsv has no layout column")

        allowed_set = {x.strip().upper() for x in allowed_layout.split(",") if x.strip()}
        before = len(df)
        df = df[df["layout"].astype(str).str.upper().isin(allowed_set)].copy()
        after = len(df)
        print(f"Filtered sample_fastp_summary by layout={allowed_layout}: {before} -> {after}")

    required = ["sample_id", "project"]
    for c in required:
        if c not in df.columns:
            raise SystemExit(f"ERROR: sample_fastp_summary.tsv does not contain column: {c}")

    if "species_merged" in df.columns:
        species_col = "species_merged"
    elif "species_original" in df.columns:
        species_col = "species_original"
    else:
        raise SystemExit("ERROR: sample_fastp_summary.tsv does not contain species_merged or species_original")

    df["binomial_key"] = df[species_col].map(canonical_binomial_species)
    df["sample_id"] = df["sample_id"].map(clean_value)
    df["project"] = df["project"].map(clean_value)

    numeric_cols = [
        "raw_reads_total",
        "raw_bases_total",
        "clean_reads_total",
        "clean_bases_total",
        "raw_read_pairs_or_se_reads",
        "clean_read_pairs_or_se_reads",
    ]

    for c in numeric_cols:
        if c not in df.columns:
            df[c] = 0
        df[c] = to_num(df[c])

    out = (
        df.groupby(["binomial_key", "project", "sample_id"], dropna=False)
        .agg(
            raw_reads_total=("raw_reads_total", "sum"),
            raw_bases_total=("raw_bases_total", "sum"),
            clean_reads_total=("clean_reads_total", "sum"),
            clean_bases_total=("clean_bases_total", "sum"),
            raw_read_pairs_or_se_reads=("raw_read_pairs_or_se_reads", "sum"),
            clean_read_pairs_or_se_reads=("clean_read_pairs_or_se_reads", "sum"),
        )
        .reset_index()
    )

    return out


def build_selected_gz_sample_table(base_dir, sample_fastp_summary_path=None, allowed_layout=""):
    inv = scan_fastq_files(base_dir)

    if inv.empty:
        return pd.DataFrame()

    allowed_keys = load_allowed_sample_keys_from_sample_fastp(
        sample_fastp_summary_path,
        allowed_layout=allowed_layout
    )

    inv = filter_inventory_by_allowed_sample_keys(
        inv,
        allowed_keys,
        label="Top5 gz fallback inventory"
    )

    if inv.empty:
        return pd.DataFrame()

    source_species_with_raw = set(
        inv.loc[inv["source_class"] == "raw", "source_species_key"].unique()
    )

    inv["use_for_top_sample"] = False

    inv.loc[
        inv["source_species_key"].isin(source_species_with_raw) &
        (inv["source_class"] == "raw"),
        "use_for_top_sample"
    ] = True

    inv.loc[
        (~inv["source_species_key"].isin(source_species_with_raw)) &
        (inv["source_class"] == "clean"),
        "use_for_top_sample"
    ] = True

    selected = inv[inv["use_for_top_sample"]].copy()

    if selected.empty:
        return pd.DataFrame()

    out = (
        selected.groupby(["binomial_name", "project", "run"], dropna=False)
        .agg(
            selected_gz_fastq_size_GB=("size_GB", "sum"),
            selected_gz_source=("source_class", join_unique),
        )
        .reset_index()
    )

    out = out.rename(columns={
        "binomial_name": "binomial_key",
        "run": "sample_id",
    })

    out["binomial_key"] = out["binomial_key"].map(canonical_binomial_species)
    out["project"] = out["project"].map(clean_value)
    out["sample_id"] = out["sample_id"].map(clean_value)

    return out


def choose_balanced_top_samples(g, top_n=5, rank_col="rank_value"):
    g = g.copy()
    g["project_for_select"] = g["project"].map(lambda x: clean_value(x) or "Unknown_project")

    project_lists = {}

    for project, sub in g.groupby("project_for_select", dropna=False):
        sub = sub.sort_values(rank_col, ascending=False)
        project_lists[project] = list(sub.index)

    selected_indices = []

    while len(selected_indices) < top_n:
        candidates = []

        for project, idx_list in project_lists.items():
            if len(idx_list) == 0:
                continue
            idx = idx_list[0]
            candidates.append((project, idx, g.loc[idx, rank_col]))

        if len(candidates) == 0:
            break

        candidates = sorted(candidates, key=lambda x: x[2], reverse=True)

        for project, idx, value in candidates:
            if len(selected_indices) >= top_n:
                break

            selected_indices.append(idx)
            project_lists[project].pop(0)

    selected = g.loc[selected_indices].copy()
    return selected


def add_balanced_top_sample_averages(
    species_summary,
    base_dir,
    sample_fastp_summary_path,
    top_n=5,
    rank_by="raw_bases",
    allowed_layout="",
    outdir=None,
    write_debug=False,
):
    df = species_summary.copy()
    df["binomial_key"] = df["Binomial name"].map(canonical_binomial_species)

    fastp_sample = read_sample_fastp_summary(
        sample_fastp_summary_path,
        allowed_layout=allowed_layout
    )

    gz_sample = build_selected_gz_sample_table(
        base_dir,
        sample_fastp_summary_path=sample_fastp_summary_path,
        allowed_layout=allowed_layout
    )

    if fastp_sample.empty and gz_sample.empty:
        print("WARNING: no sample-level fastp or gz table available; top sample averages will use species-level fallback only.")
        sample_table = pd.DataFrame(columns=[
            "binomial_key", "project", "sample_id",
            "selected_gz_fastq_size_GB",
            "raw_reads_total", "raw_bases_total",
            "clean_reads_total", "clean_bases_total",
            "raw_read_pairs_or_se_reads", "clean_read_pairs_or_se_reads"
        ])
    elif fastp_sample.empty:
        sample_table = gz_sample.copy()
    elif gz_sample.empty:
        sample_table = fastp_sample.copy()
    else:
        sample_table = gz_sample.merge(
            fastp_sample,
            how="outer",
            on=["binomial_key", "project", "sample_id"]
        )

    for c in [
        "selected_gz_fastq_size_GB",
        "raw_reads_total",
        "raw_bases_total",
        "clean_reads_total",
        "clean_bases_total",
        "raw_read_pairs_or_se_reads",
        "clean_read_pairs_or_se_reads",
    ]:
        if c not in sample_table.columns:
            sample_table[c] = 0
        sample_table[c] = to_num(sample_table[c])

    if "project" not in sample_table.columns:
        sample_table["project"] = ""
    if "sample_id" not in sample_table.columns:
        sample_table["sample_id"] = ""

    sample_table["project"] = sample_table["project"].map(clean_value)
    sample_table["sample_id"] = sample_table["sample_id"].map(clean_value)

    has_real_fastp = (
        (sample_table["raw_bases_total"] > 0) |
        (sample_table["clean_bases_total"] > 0) |
        (sample_table["raw_reads_total"] > 0) |
        (sample_table["clean_reads_total"] > 0)
    )

    has_gz = sample_table["selected_gz_fastq_size_GB"] > 0

    training = sample_table[has_real_fastp & has_gz].copy()

    metric_cols = [
        "raw_reads_total",
        "raw_read_pairs_or_se_reads",
        "raw_bases_total",
        "clean_reads_total",
        "clean_read_pairs_or_se_reads",
        "clean_bases_total",
    ]

    ratio = {}

    for c in metric_cols:
        ratio[c] = safe_median_ratio(
            training[c],
            training["selected_gz_fastq_size_GB"]
        )

    sample_table["sample_metric_source"] = np.where(
        has_real_fastp,
        "sample_fastp_summary_direct",
        "needs_sample_gz_fallback"
    )

    for c in metric_cols:
        final_col = c + "_final"

        sample_table[final_col] = np.where(
            sample_table[c] > 0,
            sample_table[c],
            np.nan
        )

        missing_metric = sample_table[final_col].isna() & has_gz

        if pd.notna(ratio.get(c, np.nan)) and ratio[c] > 0:
            sample_table.loc[missing_metric, final_col] = (
                sample_table.loc[missing_metric, "selected_gz_fastq_size_GB"] * ratio[c]
            )

    any_estimated = (
        sample_table["sample_metric_source"].eq("needs_sample_gz_fallback") &
        sample_table[[c + "_final" for c in metric_cols]].notna().any(axis=1)
    )

    sample_table.loc[any_estimated, "sample_metric_source"] = (
        "sample_gz_size_estimated_using_median_fastp_ratio"
    )

    sample_table.loc[
        sample_table[[c + "_final" for c in metric_cols]].isna().all(axis=1),
        "sample_metric_source"
    ] = "no_sample_level_metric_available"

    rank_map = {
        "raw_bases": "raw_bases_total_final",
        "clean_bases": "clean_bases_total_final",
        "raw_reads": "raw_reads_total_final",
        "clean_reads": "clean_reads_total_final",
    }

    rank_col = rank_map.get(rank_by, "raw_bases_total_final")

    sample_table["rank_value"] = to_num(sample_table[rank_col])

    sample_table.loc[
        sample_table["rank_value"].le(0) & sample_table["clean_bases_total_final"].notna(),
        "rank_value"
    ] = sample_table.loc[
        sample_table["rank_value"].le(0) & sample_table["clean_bases_total_final"].notna(),
        "clean_bases_total_final"
    ]

    sample_table.loc[
        sample_table["rank_value"].le(0) & sample_table["selected_gz_fastq_size_GB"].gt(0),
        "rank_value"
    ] = sample_table.loc[
        sample_table["rank_value"].le(0) & sample_table["selected_gz_fastq_size_GB"].gt(0),
        "selected_gz_fastq_size_GB"
    ]

    sample_table_valid = sample_table[
        sample_table["rank_value"].notna() &
        sample_table["rank_value"].gt(0)
    ].copy()

    top_rows = []

    for binomial_key, g in sample_table_valid.groupby("binomial_key", dropna=False):
        selected = choose_balanced_top_samples(
            g=g,
            top_n=top_n,
            rank_col="rank_value"
        )

        if selected.empty:
            continue

        project_counts = (
            selected["project"]
            .map(lambda x: clean_value(x) or "Unknown_project")
            .value_counts()
            .sort_index()
        )

        project_distribution = "; ".join([
            f"{project}={count}" for project, count in project_counts.items()
        ])

        sample_ids = []
        for _, r in selected.iterrows():
            project = clean_value(r.get("project", "")) or "Unknown_project"
            sample_id = clean_value(r.get("sample_id", ""))
            if sample_id:
                sample_ids.append(f"{project}:{sample_id}")

        source_note = join_unique(selected["sample_metric_source"])

        top_rows.append({
            "binomial_key": binomial_key,
            "Top5 balanced samples used": len(selected),
            "Top5 balanced project distribution": project_distribution,
            "Top5 balanced sample IDs": "; ".join(sample_ids),

            "Top5 avg raw reads": selected["raw_reads_total_final"].mean(),
            "Top5 avg raw read pairs or SE reads": selected["raw_read_pairs_or_se_reads_final"].mean(),
            "Top5 avg raw bases total (bp)": selected["raw_bases_total_final"].mean(),
            "Top5 avg raw sequencing size (Gb bases)": selected["raw_bases_total_final"].mean() / 1e9,

            "Top5 avg clean reads": selected["clean_reads_total_final"].mean(),
            "Top5 avg clean read pairs or SE reads": selected["clean_read_pairs_or_se_reads_final"].mean(),
            "Top5 avg clean bases total (bp)": selected["clean_bases_total_final"].mean(),
            "Top5 avg clean sequencing size (Gb bases)": selected["clean_bases_total_final"].mean() / 1e9,

            "Top5 avg source": source_note,
        })

    if len(top_rows) > 0:
        top_df = pd.DataFrame(top_rows)
    else:
        top_df = pd.DataFrame(columns=[
            "binomial_key",
            "Top5 balanced samples used",
            "Top5 balanced project distribution",
            "Top5 balanced sample IDs",
            "Top5 avg raw reads",
            "Top5 avg raw read pairs or SE reads",
            "Top5 avg raw bases total (bp)",
            "Top5 avg raw sequencing size (Gb bases)",
            "Top5 avg clean reads",
            "Top5 avg clean read pairs or SE reads",
            "Top5 avg clean bases total (bp)",
            "Top5 avg clean sequencing size (Gb bases)",
            "Top5 avg source",
        ])

    merged = df.merge(top_df, how="left", on="binomial_key")

    fallback_mask = merged["Top5 balanced samples used"].isna()

    for c in [
        "Sample number",
        "Raw reads total",
        "Raw read pairs or SE reads",
        "Raw bases total (bp)",
        "Clean reads total",
        "Clean read pairs or SE reads",
        "Clean bases total (bp)",
    ]:
        if c not in merged.columns:
            merged[c] = 0
        merged[c] = to_num(merged[c])

    sample_n = merged["Sample number"].replace(0, np.nan)

    merged.loc[fallback_mask, "Top5 balanced samples used"] = merged.loc[fallback_mask, "Sample number"]
    merged.loc[fallback_mask, "Top5 balanced project distribution"] = "species_total_average_fallback"
    merged.loc[fallback_mask, "Top5 balanced sample IDs"] = ""

    merged.loc[fallback_mask, "Top5 avg raw reads"] = (
        merged.loc[fallback_mask, "Raw reads total"] / sample_n[fallback_mask]
    )

    merged.loc[fallback_mask, "Top5 avg raw read pairs or SE reads"] = (
        merged.loc[fallback_mask, "Raw read pairs or SE reads"] / sample_n[fallback_mask]
    )

    merged.loc[fallback_mask, "Top5 avg raw bases total (bp)"] = (
        merged.loc[fallback_mask, "Raw bases total (bp)"] / sample_n[fallback_mask]
    )

    merged.loc[fallback_mask, "Top5 avg raw sequencing size (Gb bases)"] = (
        merged.loc[fallback_mask, "Top5 avg raw bases total (bp)"] / 1e9
    )

    merged.loc[fallback_mask, "Top5 avg clean reads"] = (
        merged.loc[fallback_mask, "Clean reads total"] / sample_n[fallback_mask]
    )

    merged.loc[fallback_mask, "Top5 avg clean read pairs or SE reads"] = (
        merged.loc[fallback_mask, "Clean read pairs or SE reads"] / sample_n[fallback_mask]
    )

    merged.loc[fallback_mask, "Top5 avg clean bases total (bp)"] = (
        merged.loc[fallback_mask, "Clean bases total (bp)"] / sample_n[fallback_mask]
    )

    merged.loc[fallback_mask, "Top5 avg clean sequencing size (Gb bases)"] = (
        merged.loc[fallback_mask, "Top5 avg clean bases total (bp)"] / 1e9
    )

    merged.loc[fallback_mask, "Top5 avg source"] = (
        "species_total_average_fallback_no_sample_level_ranking"
    )

    int_cols = [
        "Top5 balanced samples used",
        "Top5 avg raw reads",
        "Top5 avg raw read pairs or SE reads",
        "Top5 avg raw bases total (bp)",
        "Top5 avg clean reads",
        "Top5 avg clean read pairs or SE reads",
        "Top5 avg clean bases total (bp)",
    ]

    for c in int_cols:
        if c in merged.columns:
            merged[c] = merged[c].map(format_int_or_blank)

    gb_cols = [
        "Top5 avg raw sequencing size (Gb bases)",
        "Top5 avg clean sequencing size (Gb bases)",
    ]

    for c in gb_cols:
        if c in merged.columns:
            merged[c] = merged[c].map(format_gb_or_blank)

    merged = merged.drop(columns=["binomial_key"], errors="ignore")

    if write_debug and outdir is not None:
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        sample_table.to_csv(outdir / "debug_top5_sample_level_metrics.csv", index=False)
        sample_table_valid.to_csv(outdir / "debug_top5_valid_ranked_samples.csv", index=False)
        top_df.to_csv(outdir / "debug_top5_selected_summary.csv", index=False)

    return merged


def add_machado_ddi_only(df, machado_path):
    machado_path = Path(machado_path)

    if not machado_path.exists():
        print(f"WARNING: Machado file not found: {machado_path}")
        df["Diet diversity (DDI)"] = ""
        unmatched = df.copy()
        unmatched["_match_key"] = unmatched["Binomial name"].map(normalize_species_to_machado_key)
        return df, unmatched

    machado_df = pd.read_csv(machado_path, dtype=str)
    machado_df.columns = [c.strip() for c in machado_df.columns]

    if "Species" not in machado_df.columns:
        raise SystemExit("ERROR: Machado table does not contain column: Species")

    ddi_col = None
    for c in machado_df.columns:
        norm = re.sub(r"[^a-z0-9]+", "", c.lower())
        if norm == "dietdiversityddi":
            ddi_col = c
            break

    if ddi_col is None:
        raise SystemExit("ERROR: Machado table does not contain column: Diet diversity (DDI)")

    out = df.copy()

    old_cols = [
        "Diet diversity (DDI)",
        "Machado matched species",
        "Machado diet breadth",
        "Machado diet disparity (PSV)",
        "Machado trophic guild",
    ]
    out = out.drop(columns=[c for c in old_cols if c in out.columns], errors="ignore")

    out["_match_key"] = out["Binomial name"].map(normalize_species_to_machado_key)
    machado_df["_match_key"] = machado_df["Species"].map(normalize_species_to_machado_key)

    machado_sub = machado_df[["_match_key", ddi_col]].copy()
    machado_sub = machado_sub.rename(columns={ddi_col: "Diet diversity (DDI)"})
    machado_sub["Diet diversity (DDI)"] = machado_sub["Diet diversity (DDI)"].map(clean_value)

    machado_sub = machado_sub[machado_sub["_match_key"].map(clean_value).ne("")]
    machado_sub = machado_sub.sort_values(["_match_key", "Diet diversity (DDI)"])
    machado_sub = machado_sub.drop_duplicates("_match_key", keep="first")

    merged = out.merge(machado_sub, how="left", on="_match_key")

    unmatched = merged[
        merged["Diet diversity (DDI)"].isna() |
        merged["Diet diversity (DDI)"].map(clean_value).eq("")
    ].copy()

    merged = merged.drop(columns=["_match_key"], errors="ignore")

    cols = list(merged.columns)
    if "Diet diversity (DDI)" in cols:
        cols.remove("Diet diversity (DDI)")
        if "Diet" in cols:
            pos = cols.index("Diet") + 1
            cols.insert(pos, "Diet diversity (DDI)")
        else:
            cols.append("Diet diversity (DDI)")
        merged = merged[cols]

    return merged, unmatched


def reorder_final_columns(df):
    preferred = [
        "Binomial name",
        "Source species names",
        "Project number",
        "Bioproject/DDBJ",
        "Sequencing platform/model",
        "Sample number",
        "Raw sample number on disk",
        "Clean sample number on disk",
        "Country",
        "State/Location",
        "Continent",
        "Diet",
        "Diet diversity (DDI)",
        "Living environment",
        "Reported overall data size (GB)",
        "Size source for reported overall data",
        "Source_used",
        "Raw fastq size on disk (GB)",
        "Clean fastq size on disk (GB)",
        "Article",
        "Diet source",

        "Raw reads total",
        "Raw read pairs or SE reads",
        "Raw bases total (bp)",
        "Raw sequencing size (Gb bases)",
        "Clean reads total",
        "Clean read pairs or SE reads",
        "Clean bases total (bp)",
        "Clean sequencing size (Gb bases)",
        "Fastp / size estimate source",

        "Top5 balanced samples used",
        "Top5 balanced project distribution",
        "Top5 balanced sample IDs",
        "Top5 avg raw reads",
        "Top5 avg raw read pairs or SE reads",
        "Top5 avg raw bases total (bp)",
        "Top5 avg raw sequencing size (Gb bases)",
        "Top5 avg clean reads",
        "Top5 avg clean read pairs or SE reads",
        "Top5 avg clean bases total (bp)",
        "Top5 avg clean sequencing size (Gb bases)",
        "Top5 avg source",
    ]

    existing = [c for c in preferred if c in df.columns]
    remaining = [c for c in df.columns if c not in existing]

    return df[existing + remaining]


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--base-dir", default=".")
    parser.add_argument("--meta-run", default="all_bioproject_runs_metadata_sort.csv")
    parser.add_argument("--meta-species", default="ENA_gut_mirco_source_Primates.csv")
    parser.add_argument("--fastp-summary", default="01_clean_data/05_summary_tables/species_summary_merged_level.tsv")
    parser.add_argument("--sample-fastp-summary", default="01_clean_data/05_summary_tables/sample_fastp_summary.tsv")
    parser.add_argument("--top-n-samples", type=int, default=5)
    parser.add_argument("--top-sample-rank-by", default="raw_bases", choices=["raw_bases", "clean_bases", "raw_reads", "clean_reads"])
    parser.add_argument("--strict-layout", default="", help="Strictly keep only samples with this layout across the whole final table, e.g. PE. Empty means no strict filtering.")
    parser.add_argument("--machado", default="Machado_Diet_Breadth_TableS1.csv")
    parser.add_argument("--outdir", default="00_dataset_summary_species_priority")
    parser.add_argument("--out", default="09_final_species_dataset_summary.csv")
    parser.add_argument("--unmatched-ddi-out", default="09_unmatched_species_for_Machado_DDI.csv")
    parser.add_argument("--write-debug", action="store_true")

    args = parser.parse_args()

    base_dir = Path(args.base_dir).resolve()

    meta_run = Path(args.meta_run)
    if not meta_run.is_absolute():
        meta_run = base_dir / meta_run

    meta_species = Path(args.meta_species)
    if not meta_species.is_absolute():
        meta_species = base_dir / meta_species

    fastp_summary = Path(args.fastp_summary)
    if not fastp_summary.is_absolute():
        fastp_summary = base_dir / fastp_summary

    sample_fastp_summary = Path(args.sample_fastp_summary)
    if not sample_fastp_summary.is_absolute():
        sample_fastp_summary = base_dir / sample_fastp_summary

    machado = Path(args.machado)
    if not machado.is_absolute():
        machado = base_dir / machado

    outdir = Path(args.outdir)
    if not outdir.is_absolute():
        outdir = base_dir / outdir

    outdir.mkdir(parents=True, exist_ok=True)

    out_path = outdir / args.out
    unmatched_path = outdir / args.unmatched_ddi_out

    print("==============================================")
    print("Building final species dataset summary")
    print("==============================================")
    print(f"Base directory: {base_dir}")
    print(f"Run metadata: {meta_run}")
    print(f"Species metadata: {meta_species}")
    print(f"Fastp summary: {fastp_summary}")
    print(f"Sample fastp summary: {sample_fastp_summary}")
    print(f"Machado DDI table: {machado}")
    print(f"Top-N samples: {args.top_n_samples}")
    print(f"Top sample rank by: {args.top_sample_rank_by}")
    print(f"Strict layout filter: {args.strict_layout if args.strict_layout else 'OFF'}")
    print("")

    species_summary = build_base_species_summary(
        base_dir=base_dir,
        meta_run=meta_run,
        meta_species=meta_species,
        outdir=outdir,
        sample_fastp_summary_path=sample_fastp_summary,
        allowed_layout=args.strict_layout,
        write_debug=args.write_debug,
    )

    species_summary = add_fastp_reads_bases(
        species_summary=species_summary,
        fastp_summary_path=fastp_summary,
    )

    species_summary = add_balanced_top_sample_averages(
        species_summary=species_summary,
        base_dir=base_dir,
        sample_fastp_summary_path=sample_fastp_summary,
        top_n=args.top_n_samples,
        rank_by=args.top_sample_rank_by,
        allowed_layout=args.strict_layout,
        outdir=outdir,
        write_debug=args.write_debug,
    )

    final_df, unmatched = add_machado_ddi_only(
        df=species_summary,
        machado_path=machado,
    )

    final_df = reorder_final_columns(final_df)

    if "Reported overall data size (GB)" in final_df.columns:
        final_df["_sort_size"] = to_num(final_df["Reported overall data size (GB)"])
        final_df = final_df.sort_values("_sort_size", ascending=False)
        final_df = final_df.drop(columns=["_sort_size"], errors="ignore")

    final_df.to_csv(out_path, index=False)

    unmatched_cols = [
        "Binomial name",
        "Source species names",
        "Bioproject/DDBJ",
        "Diet",
        "_match_key",
    ]
    unmatched_cols = [c for c in unmatched_cols if c in unmatched.columns]

    unmatched_report = unmatched[unmatched_cols].copy()
    unmatched_report.to_csv(unmatched_path, index=False)

    print("==============================================")
    print("Final summary finished")
    print("==============================================")
    print(f"Final table: {out_path}")
    print(f"Unmatched Machado DDI report: {unmatched_path}")
    print("")
    print(f"Total species: {len(final_df)}")
    print(f"Matched Machado DDI: {len(final_df) - len(unmatched_report)}")
    print(f"Unmatched Machado DDI: {len(unmatched_report)}")
    print("")

    if len(unmatched_report) > 0:
        print("WARNING: species not matched to Machado Diet diversity (DDI):")
        print("----------------------------------------------")
        for _, r in unmatched_report.iterrows():
            print(
                f"{r.get('Binomial name', '')}\t"
                f"match_key={r.get('_match_key', '')}\t"
                f"source_species={r.get('Source species names', '')}"
            )
    else:
        print("All species matched to Machado Diet diversity (DDI).")

    print("==============================================")


if __name__ == "__main__":
    main()
