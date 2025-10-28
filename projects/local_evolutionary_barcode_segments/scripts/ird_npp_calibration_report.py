#!/usr/bin/env python3
# ====================================================================================
# Script: ird_npp_calibration_report.py
# Purpose: Publication-ready statistics and visualization-config helper for wide
#          phylogenetic profile matrices (NPP/LPP). Computes distribution stats,
#          tail/NA/zero metrics, and suggests z-score clipping thresholds for heatmaps.
# Author: Shalev Yaacov  |  refactored for publication
# Created: 2025-10-28
#
# Macro (what it does):
#   1) Loads a wide gene×species matrix (CSV/TSV/Parquet/Feather) and detects the gene column.
#   2) Computes summary statistics over the numeric block (min/max/percentiles/mean/std/skew/kurtosis),
#      NA counts, zero fractions, and constant row/column counts.
#   3) Suggests z-score clipping (z_clip) by tail mass (e.g., 0.5% each tail) with optional rounding step.
#   4) Writes timestamped outputs to output_data/: summary.csv/json, percentiles.csv, per-row/col
#      problem lists, and a run log. Optionally emits small previews.
#
# Notes:
#   • Use --tail-mass to control z_clip suggestion; defaults are conservative for visualization.
# ====================================================================================

# ==========================================
# 0) Imports & constants
# ==========================================
import argparse
import logging
import os
import sys
import json
from datetime import datetime
from typing import Optional, Sequence, Tuple, List

import numpy as np
import pandas as pd

CANDIDATE_GENE_COLS = [
    "gene", "genes", "symbol", "gene_symbol", "hgnc_symbol",
    "Gene", "Genes", "SYMBOL", "GeneSymbol", "HGNC_symbol"
]


# ==========================================
# 1) Utilities: time, io, logging, detection
# ==========================================

def ts() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def setup_logger(logfile: str) -> logging.Logger:
    logger = logging.getLogger("NPPStats")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    fmt = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)
    logger.addHandler(ch)

    fh = logging.FileHandler(logfile, mode="w", encoding="utf-8")
    fh.setLevel(logging.INFO)
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    return logger


def detect_col(df: pd.DataFrame, user_col: Optional[str], candidates: Sequence[str]) -> str:
    if user_col:
        if user_col in df.columns:
            return user_col
        for c in df.columns:
            if c.lower() == user_col.lower():
                return c
        raise ValueError(f"Requested gene column '{user_col}' not found.")
    lower_map = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    return df.columns[0]


def load_matrix(path: str, delimiter: Optional[str]) -> pd.DataFrame:
    ext = os.path.splitext(path)[1].lower()
    if ext == ".parquet":
        return pd.read_parquet(path)
    if ext in (".feather", ".ft"):
        return pd.read_feather(path)
    if ext in (".xlsx", ".xls"):
        return pd.read_excel(path)
    if delimiter is None:
        delimiter = "\t" if ext in (".tsv", ".tab", ".txt") else ","
    return pd.read_csv(path, sep=delimiter, low_memory=False)


# ==========================================
# 2) Core: stats & z_clip suggestion
# ==========================================

def nanpercentile(a: np.ndarray, q: float) -> float:
    return float(np.nanpercentile(a, q))


def suggest_zclip(values: np.ndarray, tail_mass: float, round_to: float) -> Tuple[float, float]:
    """Suggest symmetric z-clip around 0 by tail mass per side.
    Example: tail_mass=0.5 → clip at [p0.5, p99.5]. Round up to the nearest step.
    """
    if values.size == 0:
        return (0.0, 0.0)
    lo_q = tail_mass
    hi_q = 100.0 - tail_mass
    lo = np.nanpercentile(values, lo_q)
    hi = np.nanpercentile(values, hi_q)
    # enforce symmetry by taking the max absolute magnitude
    m = max(abs(lo), abs(hi))
    if round_to and round_to > 0:
        m = round_to * np.ceil(m / round_to)
    return (-float(m), float(m))


def compute_stats(numeric_df: pd.DataFrame) -> Tuple[dict, pd.DataFrame, pd.DataFrame]:
    # Flatten valid values for global summaries
    vals = numeric_df.values
    finite_mask = np.isfinite(vals)
    valid = vals[finite_mask]

    if valid.size == 0:
        raise ValueError("No finite numeric values found in the matrix.")

    # Global stats
    g_min = float(np.nanmin(valid))
    g_max = float(np.nanmax(valid))
    p01  = nanpercentile(valid, 1)
    p50  = nanpercentile(valid, 50)
    p99  = nanpercentile(valid, 99)
    g_mean = float(np.nanmean(valid))
    g_std  = float(np.nanstd(valid))

    # Fisher's definition of skew/kurtosis (using pandas for robustness)
    s = pd.Series(valid)
    g_skew = float(s.skew())
    g_kurt = float(s.kurt())

    # NA/zero/constant counts (rows/cols)
    row_na = numeric_df.isna().sum(axis=1)
    col_na = numeric_df.isna().sum(axis=0)
    row_zero_frac = (numeric_df == 0).sum(axis=1) / numeric_df.shape[1]
    col_zero_frac = (numeric_df == 0).sum(axis=0) / numeric_df.shape[0]
    row_var = numeric_df.var(axis=1, ddof=0)
    col_var = numeric_df.var(axis=0, ddof=0)

    problem_rows = pd.DataFrame({
        "row_index": numeric_df.index,
        "row_na_count": row_na,
        "row_zero_fraction": row_zero_frac,
        "row_variance": row_var,
        "is_constant_row": row_var == 0,
    })
    problem_cols = pd.DataFrame({
        "column": numeric_df.columns,
        "col_na_count": col_na.values,
        "col_zero_fraction": col_zero_frac.values,
        "col_variance": col_var.values,
        "is_constant_col": (col_var.values == 0),
    })

    summary = {
        "rows": int(numeric_df.shape[0]),
        "cols": int(numeric_df.shape[1]),
        "min": g_min,
        "p01": p01,
        "median": p50,
        "p99": p99,
        "max": g_max,
        "mean": g_mean,
        "std": g_std,
        "skew": g_skew,
        "kurtosis": g_kurt,
        "rows_constant": int((row_var == 0).sum()),
        "cols_constant": int((col_var == 0).sum()),
        "any_na": bool(np.isnan(vals).any()),
        "any_zero": bool((vals == 0).any()),
    }
    return summary, problem_rows, problem_cols


# ==========================================
# 3) CLI & main
# ==========================================

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Statistics and z-clip suggestion for NPP/LPP matrices; outputs CSV/JSON reports."
        )
    )
    # Inputs
    p.add_argument("--input", required=True, help="Path to matrix (CSV/TSV/Parquet/Feather/XLSX)")
    p.add_argument("--gene-col", default=None, help="Gene column name in the matrix (optional)")
    p.add_argument("--delimiter", default=None, help="Delimiter for CSV/TSV input")
    p.add_argument("--matrix-name", default="NPP", help="Label for reports (e.g., 'NPP' or 'LPP')")

    # Percentiles & z-clip config
    p.add_argument("--percentiles", default="1,5,25,50,75,95,99", help="Comma-separated percentiles to compute")
    p.add_argument("--tail-mass", type=float, default=0.5, help="Tail mass %% per side for z_clip suggestion (default: 0.5)")
    p.add_argument("--round-step", type=float, default=0.25, help="Round |z_clip| up to nearest step (default: 0.25)")

    # Outputs
    p.add_argument("--outdir", default="output_data", help="Output directory")
    p.add_argument("--preview-rows", type=int, default=10, help="Rows to print into preview.csv")
    return p


def main() -> int:
    args = build_parser().parse_args()

    ensure_dir(args.outdir)
    run_ts = ts()
    base = f"{args.matrix_name}_stats_{run_ts}"
    logfile = os.path.join(args.outdir, f"{base}.log")
    logger = setup_logger(logfile)

    logger.info("[START] Matrix statistics & z-clip suggestion")

    # 3.1 Load
    df = load_matrix(args.input, args.delimiter)
    logger.info(f"[LOAD] Input shape: {df.shape[0]}×{df.shape[1]}")

    # 3.2 Detect gene column and peel numeric block
    gene_col = detect_col(df, args.gene_col, CANDIDATE_GENE_COLS)
    if gene_col not in df.columns:
        raise ValueError("Gene column not found after detection.")
    cols = [gene_col] + [c for c in df.columns if c != gene_col]
    df = df[cols]
    numeric_df = df.drop(columns=[gene_col], errors="ignore").select_dtypes(include=[np.number])
    if numeric_df.empty:
        raise ValueError("No numeric columns found in the matrix.")

    # 3.3 Compute core stats
    summary, problem_rows, problem_cols = compute_stats(numeric_df)

    # 3.4 Percentiles table (user-defined)
    try:
        qs = [float(q.strip()) for q in args.percentiles.split(",") if q.strip() != ""]
    except Exception:
        qs = [1.0, 5.0, 25.0, 50.0, 75.0, 95.0, 99.0]
    vals = numeric_df.values
    finite_vals = vals[np.isfinite(vals)]
    pct_rows = []
    for q in qs:
        pct_rows.append({"percentile": q, "value": float(np.nanpercentile(finite_vals, q))})
    pct_df = pd.DataFrame(pct_rows)

    # 3.5 z_clip suggestion by tail mass
    zlo, zhi = suggest_zclip(finite_vals, tail_mass=float(args.tail_mass), round_to=float(args.round_step))
    summary.update({"suggested_zclip_low": zlo, "suggested_zclip_high": zhi, "tail_mass_per_side_percent": float(args.tail_mass), "round_step": float(args.round_step)})

    # 3.6 Save outputs
    out_summary_csv = os.path.join(args.outdir, f"{base}_summary.csv")
    out_summary_json= os.path.join(args.outdir, f"{base}_summary.json")
    out_pct_csv     = os.path.join(args.outdir, f"{base}_percentiles.csv")
    out_prow_csv    = os.path.join(args.outdir, f"{base}_problem_rows.csv")
    out_pcol_csv    = os.path.join(args.outdir, f"{base}_problem_cols.csv")
    out_preview_csv = os.path.join(args.outdir, f"{base}_preview.csv")

    logger.info("[SAVE] Writing CSV/JSON artifacts…")
    pd.DataFrame([summary]).to_csv(out_summary_csv, index=False)
    with open(out_summary_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    pct_df.to_csv(out_pct_csv, index=False)
    problem_rows.to_csv(out_prow_csv, index=False)
    problem_cols.to_csv(out_pcol_csv, index=False)
    df.head(max(1, args.preview_rows)).to_csv(out_preview_csv, index=False)

    logger.info(
        f"[END] Done. Suggested z_clip = [{zlo:.2f}, {zhi:.2f}] (tail_mass={args.tail_mass}% per side). ✅"
    )
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:
        print(f"[FATAL] {e}")
        sys.exit(2)
