#!/usr/bin/env python3
# ====================================================================================
# Script: phylo_profile_matrix_qc.py
# Purpose: Fast, memory‑aware QC/inspection for large phylogenetic profile matrices
#          (LPP/NPP). Produces a concise summary, problem lists, and optional
#          downcasted Parquet output. Designed for portfolio publication.
# Author: Shalev Yaacov  |  Maintainer: (you)
# Created: 2025-10-28
#
# What this script does (macro):
#   1) Loads a wide gene × species matrix (CSV/TSV/Parquet/Feather) with a gene ID column.
#   2) Validates structure (gene column, shape, dtypes), detects NAs/zeros/constant rows/cols,
#      duplicates, and computes basic statistics (ranges, percentiles, variance).
#   3) Writes a timestamped QC bundle into output_data/: summary.csv, dtypes.csv,
#      problem_rows.csv / problem_cols.csv, preview.csv, and a run log.
#   4) Optionally writes a memory‑efficient Parquet copy (--export-parquet) with numeric downcast.
#
# Notes:
#   • Input paths are controlled via CLI flags; defaults are sensible and easy to change.
#   • Logging goes to console and to a timestamped .log file.
# ====================================================================================

# ==========================================
# 0) Imports and Constants
# ==========================================
import argparse
import logging
import os
import sys
import json
from datetime import datetime
from typing import List, Optional

import numpy as np
import pandas as pd

# Candidate column names for gene identifiers (case‑insensitive match handled later)
CANDIDATE_GENE_COLS = [
    "gene", "genes", "symbol", "gene_symbol", "hgnc_symbol",
    "Gene", "Genes", "SYMBOL", "GeneSymbol", "HGNC_symbol"
]


# ==========================================
# 1) Utilities: logger, io, helpers
# ==========================================

def make_timestamp() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def setup_logger(log_path: str) -> logging.Logger:
    logger = logging.getLogger("LPPInspect")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    # Console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter("%(asctime)s | %(levelname)s | %(message)s"))

    # File handler
    fh = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    fh.setLevel(logging.INFO)
    fh.setFormatter(logging.Formatter("%(asctime)s | %(levelname)s | %(message)s"))

    logger.addHandler(ch)
    logger.addHandler(fh)
    return logger


def infer_gene_col(df: pd.DataFrame, user_col: Optional[str]) -> str:
    """Return the gene column name. If user_col is provided, validate it; otherwise infer.
    Raises ValueError if not found."""
    if user_col is not None:
        if user_col in df.columns:
            return user_col
        # try case-insensitive match
        for c in df.columns:
            if c.lower() == user_col.lower():
                return c
        raise ValueError(f"Requested gene column '{user_col}' not found in file.")

    # infer by candidates (case-insensitive)
    lower_map = {c.lower(): c for c in df.columns}
    for cand in CANDIDATE_GENE_COLS:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    raise ValueError(
        "Could not infer gene column. Specify with --gene-col. Candidates tried: "
        + ", ".join(CANDIDATE_GENE_COLS)
    )


def read_matrix(path: str, delimiter: Optional[str], logger: logging.Logger) -> pd.DataFrame:
    """Read CSV/TSV/Parquet/Feather. If delimiter is None, infer by extension."""
    ext = os.path.splitext(path)[1].lower()
    logger.info(f"[LOAD] Reading matrix from: {path}")

    if ext in (".parquet",):
        df = pd.read_parquet(path)
    elif ext in (".feather", ".ft"):
        df = pd.read_feather(path)
    else:
        # CSV/TSV (or unknown → default to CSV)
        if delimiter is None:
            delimiter = "\t" if ext in (".tsv", ".tab") else ","
        df = pd.read_csv(path, sep=delimiter, low_memory=False)

    logger.info(f"[LOAD] Done. Shape: {df.shape[0]} rows × {df.shape[1]} cols")
    return df


def numeric_downcast(df: pd.DataFrame) -> pd.DataFrame:
    """Downcast numeric columns to reduce memory footprint without changing values."""
    out = df.copy()
    for c in out.select_dtypes(include=[np.number]).columns:
        col = out[c]
        if pd.api.types.is_float_dtype(col):
            out[c] = pd.to_numeric(col, downcast="float")
        else:
            out[c] = pd.to_numeric(col, downcast="integer")
    return out


# ==========================================
# 2) Argument Parsing
# ==========================================

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "QC/inspection for LPP/NPP matrices: validates structure, computes basic stats, "
            "identifies data issues, and writes a timestamped QC bundle."
        )
    )
    p.add_argument("--input", required=True, help="Path to matrix file (CSV/TSV/Parquet/Feather)")
    p.add_argument("--gene-col", default=None, help="Gene identifier column name (optional)")
    p.add_argument("--delimiter", default=None, help="CSV/TSV delimiter override (optional)")
    p.add_argument("--matrix-name", default="LPP", help="Label for reports (e.g., 'LPP' or 'NPP')")
    p.add_argument("--outdir", default="output_data", help="Directory for outputs")
    p.add_argument("--export-parquet", action="store_true", help="Also write a downcasted Parquet file")
    p.add_argument("--preview-rows", type=int, default=10, help="Number of rows to include in preview.csv")
    return p


# ==========================================
# 3) Main QC Routine
# ==========================================

def run_qc(args: argparse.Namespace, logger: logging.Logger) -> int:
    logger.info("[START] LPP/NPP matrix inspection")

    # 3.1 Load
    df = read_matrix(args.input, args.delimiter, logger)

    # 3.2 Detect gene column and set index
    gene_col = infer_gene_col(df, args.gene_col)
    logger.info(f"[STRUCT] Using gene column: '{gene_col}'")

    # Move gene column to first position for clarity
    if df.columns[0] != gene_col:
        cols = [gene_col] + [c for c in df.columns if c != gene_col]
        df = df[cols]

    # 3.3 Basic structure
    n_rows, n_cols = df.shape
    num_cols = df.select_dtypes(include=[np.number]).shape[1]
    non_num_cols = n_cols - num_cols
    logger.info(f"[STRUCT] Rows: {n_rows} | Cols: {n_cols} | Numeric cols: {num_cols} | Non‑numeric cols: {non_num_cols}")

    # 3.4 Duplicates and missing gene IDs
    dup_genes = df[gene_col].astype(str).duplicated(keep=False)
    n_dup = int(dup_genes.sum())
    n_missing_gene_ids = int(df[gene_col].isna().sum())
    logger.info(f"[STRUCT] Duplicate gene IDs: {n_dup} | Missing gene IDs: {n_missing_gene_ids}")

    # 3.5 Numeric block statistics
    numeric_df = df.drop(columns=[gene_col], errors="ignore").select_dtypes(include=[np.number])
    if numeric_df.empty:
        raise ValueError("No numeric columns detected after removing the gene column. Abort.")

    # NA, zero, const detection
    row_na_cnt = numeric_df.isna().sum(axis=1)
    col_na_cnt = numeric_df.isna().sum(axis=0)
    row_zero_frac = (numeric_df == 0).sum(axis=1) / numeric_df.shape[1]
    col_zero_frac = (numeric_df == 0).sum(axis=0) / numeric_df.shape[0]

    row_var = numeric_df.var(axis=1, ddof=0)
    col_var = numeric_df.var(axis=0, ddof=0)

    const_rows = row_var.index[row_var == 0].tolist()
    const_cols = col_var.index[col_var == 0].tolist()

    # value ranges & percentiles (on the full numeric block, ignoring NAs)
    valid_vals = numeric_df.values[np.isfinite(numeric_df.values)]
    vmin = float(np.nanmin(valid_vals)) if valid_vals.size else np.nan
    vmax = float(np.nanmax(valid_vals)) if valid_vals.size else np.nan
    p01, p50, p99 = [
        float(np.nanpercentile(valid_vals, q)) if valid_vals.size else np.nan
        for q in (1, 50, 99)
    ]

    # 3.6 Build summaries
    qc_summary = {
        "matrix_name": args.matrix_name,
        "input_path": os.path.abspath(args.input),
        "rows": n_rows,
        "cols": n_cols,
        "numeric_cols": int(num_cols),
        "non_numeric_cols": int(non_num_cols),
        "duplicate_gene_ids": int(n_dup),
        "missing_gene_ids": int(n_missing_gene_ids),
        "const_rows": int(len(const_rows)),
        "const_cols": int(len(const_cols)),
        "any_na": bool(numeric_df.isna().values.any()),
        "min": vmin,
        "p01": p01,
        "median": p50,
        "p99": p99,
        "max": vmax,
    }

    # dtypes overview (for transparency)
    dtype_report = pd.DataFrame({
        "column": df.columns,
        "dtype": [str(t) for t in df.dtypes]
    })

    # Problem lists
    problem_rows = pd.DataFrame({
        gene_col: df[gene_col],
        "is_duplicate_gene": dup_genes,
        "row_na_count": row_na_cnt,
        "row_zero_fraction": row_zero_frac,
        "row_variance": row_var,
        "is_constant_row": row_var == 0,
    })

    problem_cols = pd.DataFrame({
        "column": numeric_df.columns,
        "col_na_count": col_na_cnt.values,
        "col_zero_fraction": col_zero_frac.values,
        "col_variance": col_var.values,
        "is_constant_col": (col_var.values == 0),
    })

    # 3.7 Save outputs
    ts = make_timestamp()
    ensure_dir(args.outdir)
    base = f"{args.matrix_name}_qc_{ts}"

    summary_csv = os.path.join(args.outdir, f"{base}_summary.csv")
    dtypes_csv = os.path.join(args.outdir, f"{base}_dtypes.csv")
    prow_csv    = os.path.join(args.outdir, f"{base}_problem_rows.csv")
    pcol_csv    = os.path.join(args.outdir, f"{base}_problem_cols.csv")
    preview_csv = os.path.join(args.outdir, f"{base}_preview.csv")
    summary_json = os.path.join(args.outdir, f"{base}_summary.json")

    logger.info("[SAVE] Writing summary/dtypes/problem lists/preview…")
    pd.DataFrame([qc_summary]).to_csv(summary_csv, index=False)
    dtype_report.to_csv(dtypes_csv, index=False)
    problem_rows.to_csv(prow_csv, index=False)
    problem_cols.to_csv(pcol_csv, index=False)
    df.head(max(1, args.preview_rows)).to_csv(preview_csv, index=False)

    with open(summary_json, "w", encoding="utf-8") as f:
        json.dump(qc_summary, f, indent=2)

    logger.info(f"[SAVE] Done. Files created under: {os.path.abspath(args.outdir)}")

    # 3.8 Optional Parquet export (downcast numeric)
    if args.export_parquet:
        parquet_path = os.path.join(args.outdir, f"{base}_downcast.parquet")
        logger.info("[EXPORT] Creating downcasted Parquet copy…")
        df_dc = df.copy()
        df_dc[df_dc.columns.difference([gene_col])] = numeric_downcast(
            df_dc.drop(columns=[gene_col], errors="ignore")
        )
        df_dc.to_parquet(parquet_path, index=False)
        logger.info(f"[EXPORT] Parquet saved: {parquet_path}")

    logger.info("[END] QC completed successfully ✅")
    return 0


# ==========================================
# 4) Entrypoint
# ==========================================

def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    ensure_dir(args.outdir)
    ts = make_timestamp()
    log_file = os.path.join(args.outdir, f"{args.matrix_name}_qc_{ts}.log")
    logger = setup_logger(log_file)

    try:
        return run_qc(args, logger)
    except Exception as e:
        logger.exception(f"[FATAL] QC failed: {e}")
        return 2


if __name__ == "__main__":
    sys.exit(main())
