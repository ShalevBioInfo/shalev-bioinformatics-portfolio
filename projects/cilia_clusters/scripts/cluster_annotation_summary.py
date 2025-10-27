#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# File: cluster_annotation_summary.py
# Purpose: Produce clear category-level and cluster-level summaries from a
#          master table (e.g., cluster -> category annotations). Saves
#          timestamped CSVs and an Excel workbook with summary sheets, plus
#          console + logfile output suitable for publication on GitHub.
# Created: 2025-10-27
# Author: Shalev Yaacov (refactored for GitHub publication)
#
# Short description:
#  - Reads a CSV/XLSX master table containing at least a cluster ID column and
#    a category/annotation column (column names can be provided or auto-guessed).
#  - Produces:
#      1) category_summary: counts, number of supporting clusters, and percent of clusters
#      2) cluster_category_matrix: presence/absence matrix (clusters x categories)
#      3) optional top-N categories CSV
#  - Writes outputs to `--outdir` with timestamped filenames and a run log.
#
# Usage example:
#  python cluster_category_summary.py -i data/master_table.csv \
#      --cluster-col Cluster_ID --category-col Inclusion_criterion \
#      --outdir output_data --top-n 20
#
# -----------------------------------------------------------------------------

# 1. Imports and environment --------------------------------------------------
import argparse
import logging
import sys
from pathlib import Path
from datetime import datetime
from typing import Optional, List

import pandas as pd
import numpy as np

# 2. Utilities: timestamp, logging, directories -------------------------------
def timestamp_now() -> str:
    """2.1 Return compact timestamp for filenames."""
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def ensure_dir(p: Path):
    """2.2 Ensure directory exists."""
    p.mkdir(parents=True, exist_ok=True)


def setup_logger(log_path: Path):
    """2.3 Setup logging to console and file (INFO)."""
    logger = logging.getLogger("cluster_category_summary")
    logger.setLevel(logging.INFO)
    if logger.hasHandlers():
        logger.handlers.clear()

    fmt = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s", "%Y-%m-%d %H:%M:%S")

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)
    logger.addHandler(ch)

    fh = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    fh.setLevel(logging.INFO)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    return logger


# 3. I/O helpers & guessing --------------------------------------------------
def read_table(path: Path, sheet: Optional[str], logger) -> pd.DataFrame:
    """3.1 Read CSV/TSV/XLSX into DataFrame; trim column names."""
    logger.info(f"Reading input file: {path}")
    if not path.exists():
        logger.error(f"Input file not found: {path}")
        raise FileNotFoundError(path)
    sfx = path.suffix.lower()
    if sfx in (".csv", ".tsv", ".txt"):
        sep = "," if sfx == ".csv" else "\t"
        df = pd.read_csv(path, sep=sep, dtype=str, encoding="utf-8", low_memory=False)
    elif sfx in (".xls", ".xlsx"):
        df = pd.read_excel(path, sheet_name=sheet, dtype=str)
    else:
        df = pd.read_csv(path, dtype=str, encoding="utf-8", low_memory=False)
    df.columns = [str(c).strip() for c in df.columns]
    logger.info(f"Loaded table with shape: {df.shape}")
    return df


def guess_column(df: pd.DataFrame, candidates: List[str], logger=None) -> Optional[str]:
    """3.2 Guess a column name from candidate names (case-insensitive, partial)."""
    cols_lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols_lower:
            return cols_lower[cand.lower()]
    # partial match
    for cand in candidates:
        for col in df.columns:
            if cand.lower() in col.lower():
                if logger:
                    logger.info(f"Guessed column '{col}' for candidate '{cand}'.")
                return col
    return None


def normalize_strings_df(df: pd.DataFrame, logger) -> pd.DataFrame:
    """3.3 Normalize all object columns (trim, fillna -> empty string)."""
    logger.info("Normalizing string columns (strip whitespace, fill NaN -> '').")
    df2 = df.copy()
    obj_cols = df2.select_dtypes(include=["object", "string"]).columns
    for col in obj_cols:
        df2[col] = df2[col].fillna("").astype(str).str.strip()
    return df2


# 4. Summary builders --------------------------------------------------------
def build_presence_matrix(df: pd.DataFrame, cluster_col: str, category_col: str, logger) -> pd.DataFrame:
    """
    4.1 Build boolean presence matrix (clusters x categories).
    - Removes empty categories/clusters.
    """
    logger.info("Building presence matrix (cluster x category).")
    temp = df[[cluster_col, category_col]].copy()
    temp[cluster_col] = temp[cluster_col].astype(str).str.strip()
    temp[category_col] = temp[category_col].astype(str).str.strip()
    temp = temp[(temp[cluster_col] != "") & (temp[category_col] != "")]
    if temp.empty:
        logger.error("No valid rows remain after cleaning cluster/category columns.")
        raise ValueError("No valid rows for building presence matrix.")
    pres = (temp.drop_duplicates([cluster_col, category_col])
            .assign(pres=1)
            .pivot_table(index=cluster_col, columns=category_col, values='pres', fill_value=0))
    logger.info(f"Presence matrix built with shape: {pres.shape}")
    return pres


def build_category_summary(presence_df: pd.DataFrame, logger) -> pd.DataFrame:
    """
    4.2 Build category-level summary:
      - total_occurrences (sum of presence across clusters)
      - clusters_supporting (same as total_occurrences)
      - percent_of_clusters
    """
    logger.info("Computing category-level summary.")
    total_clusters = presence_df.shape[0]
    counts = presence_df.sum(axis=0).astype(int).rename("clusters_supporting")
    summary = counts.reset_index().rename(columns={"index": "category"})
    summary["clusters_supporting"] = summary["clusters_supporting"].astype(int)
    summary["percent_of_clusters"] = (summary["clusters_supporting"] / total_clusters * 100).round(2)
    summary = summary.sort_values(by="clusters_supporting", ascending=False).reset_index(drop=True)
    logger.info(f"Category summary computed for {len(summary)} categories.")
    return summary


def build_cluster_counts(presence_df: pd.DataFrame, logger) -> pd.DataFrame:
    """
    4.3 Build cluster-level counts:
      - number of categories per cluster
      - list of categories per cluster (semicolon separated)
    """
    logger.info("Computing cluster-level counts and category lists.")
    counts = presence_df.sum(axis=1).astype(int).rename("num_categories")
    # create semicolon-separated lists of categories present
    categories_list = presence_df.apply(lambda row: ";".join(row.index[row.astype(bool)].tolist()), axis=1)
    cluster_df = pd.DataFrame({
        "cluster": presence_df.index,
        "num_categories": counts.values,
        "categories_list": categories_list.values
    })
    logger.info(f"Cluster counts computed for {cluster_df.shape[0]} clusters.")
    return cluster_df


# 5. Save helpers ------------------------------------------------------------
def save_csv(df: pd.DataFrame, path: Path, logger):
    """5.1 Save DataFrame to CSV (utf-8)."""
    logger.info(f"Saving CSV: {path.name}")
    df.to_csv(path, index=False, encoding="utf-8")


def save_workbook(sheets: dict, path: Path, logger):
    """5.2 Save multiple DataFrames to an Excel workbook (sheetname -> df)."""
    try:
        logger.info(f"Saving Excel workbook: {path.name}")
        with pd.ExcelWriter(path, engine="openpyxl") as writer:
            for sheetname, df in sheets.items():
                safe_name = sheetname[:31]  # Excel sheetname limit
                df.to_excel(writer, sheet_name=safe_name, index=False)
        logger.info("Workbook saved successfully.")
    except Exception as e:
        logger.warning(f"Failed to save workbook ({e}). CSVs will be available.")


# 6. CLI ---------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Summarize categories across clusters (counts and presence matrix).")
    p.add_argument("--input", "-i", required=True, help="Path to input CSV/XLSX master table.")
    p.add_argument("--sheet", default=None, help="Sheet name/index if input is Excel.")
    p.add_argument("--cluster-col", default=None, help="Column name identifying cluster (auto-guessed if omitted).")
    p.add_argument("--category-col", default=None, help="Column name for category/annotation (auto-guessed if omitted).")
    p.add_argument("--outdir", "-o", default="output_data", help="Directory for outputs and logs.")
    p.add_argument("--top-n", type=int, default=0, help="If >0, also save top-N categories CSV.")
    p.add_argument("--min-support", type=int, default=0, help="Minimum clusters supporting a category to keep it in summary.")
    return p.parse_args()


# 7. Main --------------------------------------------------------------------
def main():
    args = parse_args()
    input_path = Path(args.input)
    outdir = Path(args.outdir)
    ensure_dir(outdir)
    ts = timestamp_now()
    log_path = outdir / f"cluster_category_summary_{ts}.log"
    logger = setup_logger(log_path)

    logger.info("=== START: cluster_category_summary ===")
    logger.info(f"Input: {input_path}")
    logger.info(f"Outdir: {outdir}")

    try:
        # 7.1 Read input
        df = read_table(input_path, args.sheet, logger)
        if df.empty:
            logger.error("Input table is empty.")
            return 1

        # 7.2 Auto-guess columns if not provided
        cluster_col = args.cluster_col
        category_col = args.category_col
        if cluster_col is None:
            cluster_col = guess_column(df, ["cluster_id", "cluster", "clusterid", "Cluster_ID"], logger)
        if category_col is None:
            category_col = guess_column(df, ["Inclusion_criterion", "category", "annotation", "Inclusion"], logger)
        if cluster_col is None or category_col is None:
            logger.error("Could not determine cluster or category columns automatically. Please pass --cluster-col and --category-col.")
            return 1
        logger.info(f"Using cluster column: '{cluster_col}'")
        logger.info(f"Using category column: '{category_col}'")

        # 7.3 Normalize strings
        df_clean = normalize_strings_df(df, logger)

        # 7.4 Build presence matrix
        pres = build_presence_matrix(df_clean, cluster_col, category_col, logger)

        # 7.5 Optionally filter low-support categories
        if args.min_support > 0:
            logger.info(f"Filtering categories with support < {args.min_support}")
            keep = pres.sum(axis=0) >= args.min_support
            pres = pres.loc[:, keep]
            logger.info(f"Presence matrix after filter: {pres.shape}")

        # 7.6 Build summaries
        cat_summary = build_category_summary(pres, logger)
        cluster_summary = build_cluster_counts(pres, logger)

        # 7.7 Save outputs
        base = input_path.stem
        cat_csv = outdir / f"{base}_category_summary_{ts}.csv"
        cluster_csv = outdir / f"{base}_cluster_summary_{ts}.csv"
        matrix_csv = outdir / f"{base}_cluster_category_matrix_{ts}.csv"
        save_csv(cat_summary, cat_csv, logger)
        save_csv(cluster_summary, cluster_csv, logger)
        # save presence matrix with clusters as first column
        pres_export = pres.copy()
        pres_export.insert(0, "cluster", pres_export.index)
        save_csv(pres_export.reset_index(drop=True), matrix_csv, logger)

        # optional top-N
        if args.top_n and args.top_n > 0:
            topn = cat_summary.head(args.top_n)
            topn_csv = outdir / f"{base}_category_top{args.top_n}_{ts}.csv"
            save_csv(topn, topn_csv, logger)
            logger.info(f"Top-{args.top_n} categories saved: {topn_csv.name}")

        # Save Excel workbook with sheets
        workbook = outdir / f"{base}_category_summary_workbook_{ts}.xlsx"
        sheets = {
            "category_summary": cat_summary,
            "cluster_summary": cluster_summary,
            "cluster_category_matrix": pres_export.reset_index(drop=True)
        }
        save_workbook(sheets, workbook, logger)

        # 7.8 Final summary
        logger.info("=== SUMMARY ===")
        logger.info(f"Category summary rows: {len(cat_summary)}")
        logger.info(f"Cluster summary rows: {len(cluster_summary)}")
        logger.info(f"Presence matrix shape: {pres.shape}")
        logger.info(f"Saved files: {cat_csv.name}, {cluster_csv.name}, {matrix_csv.name}, {workbook.name}")
        logger.info("=== END: cluster_category_summary ===")
        return 0

    except Exception as e:
        logger.exception(f"Fatal error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
