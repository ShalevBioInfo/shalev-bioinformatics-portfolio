#!/usr/bin/env python3
# ====================================================================================
# Script: report_gene_set_matrix_coverage.py
# Purpose: Publication-ready overlap report between an IRD gene list and a gene×species
#          profile matrix (NPP/LPP). Produces clean, timestamped outputs under
#          output_data/: three CSVs (present/missing/only_in_matrix), a 3-sheet Excel,
#          and a JSON/CSV summary, with console+file logging throughout.
# Author: Shalev Yaacov  |  refactored for publication
# Created: 2025-10-28
#
# Macro (what it does):
#   1) Loads a wide matrix (CSV/TSV/Parquet/Feather) and detects the gene column.
#   2) Loads an IRD gene list (CSV/TSV/XLSX) and detects the gene column (default 'Gene').
#   3) Computes set overlaps: present_in_both, missing_from_matrix, only_in_matrix.
#   4) Writes timestamped artifacts under output_data/ and a concise run log.
#
# Notes:
#   • Input paths and column names are CLI-controlled. Defaults are sensible.
#   • All outputs include a run timestamp for easy provenance.
# ====================================================================================

# ==========================================
# 0) Imports and constants
# ==========================================
import argparse
import logging
import os
import sys
from datetime import datetime
from typing import Optional, Sequence

import pandas as pd

# Candidate column names for gene identifiers (case-insensitive)
CANDIDATE_GENE_COLS = [
    "gene", "genes", "symbol", "gene_symbol", "hgnc_symbol",
    "Gene", "Genes", "SYMBOL", "GeneSymbol", "HGNC_symbol"
]


# ==========================================
# 1) Utilities: paths, time, logging
# ==========================================

def ts() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def setup_logger(logfile: str) -> logging.Logger:
    logger = logging.getLogger("IRDGeneOverlap")
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


# ==========================================
# 2) IO helpers and column detection
# ==========================================

def detect_col(df: pd.DataFrame, user_col: Optional[str], candidates: Sequence[str]) -> str:
    """Return the gene column. If user_col is provided, validate or match case-insensitively; otherwise infer."""
    if user_col:
        if user_col in df.columns:
            return user_col
        for c in df.columns:
            if c.lower() == user_col.lower():
                return c
        raise ValueError(f"Requested gene column '{user_col}' not found in input.")

    lower_map = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    # Fall back to first column as a pragmatic default
    return df.columns[0]


def load_matrix(path: str, delimiter: Optional[str], logger: logging.Logger) -> pd.DataFrame:
    """Load a genes×species matrix; do not set index yet."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Matrix file not found: {path}")

    ext = os.path.splitext(path)[1].lower()
    logger.info(f"[LOAD] Matrix: {path}")

    if ext == ".parquet":
        df = pd.read_parquet(path)
    elif ext in (".feather", ".ft"):
        df = pd.read_feather(path)
    else:
        # CSV/TSV
        if delimiter is None:
            delimiter = "\t" if ext in (".tsv", ".tab", ".txt") else ","
        df = pd.read_csv(path, sep=delimiter, low_memory=False)

    logger.info(f"[LOAD] Matrix shape: {df.shape[0]} rows × {df.shape[1]} cols")
    return df


def load_gene_table(path: str, delimiter: Optional[str], logger: logging.Logger) -> pd.DataFrame:
    """Load a table containing at least one gene-identifier column."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Genes file not found: {path}")

    ext = os.path.splitext(path)[1].lower()
    logger.info(f"[LOAD] Genes table: {path}")

    if ext in (".xlsx", ".xls"):
        gdf = pd.read_excel(path)
    else:
        if delimiter is None:
            delimiter = "\t" if ext in (".tsv", ".tab", ".txt") else ","
        gdf = pd.read_csv(path, sep=delimiter)

    logger.info(f"[LOAD] Genes table shape: {gdf.shape[0]} rows × {gdf.shape[1]} cols")
    return gdf


# ==========================================
# 3) Core logic
# ==========================================

def compute_overlap(matrix_df: pd.DataFrame, matrix_gene_col: str, genes: pd.Series) -> dict:
    """Compute present/missing/only_in_matrix lists and summary counts."""
    # Normalize gene strings (strip whitespace)
    matrix_genes = matrix_df[matrix_gene_col].astype(str).str.strip()
    user_genes   = genes.astype(str).str.strip()

    set_matrix = set(matrix_genes)
    set_genes  = set(user_genes)

    present = sorted(set_genes & set_matrix)
    missing = sorted(set_genes - set_matrix)
    only_in_matrix = sorted(set_matrix - set_genes)

    summary = {
        "user_list_size": len(set_genes),
        "matrix_index_size": len(set_matrix),
        "present_in_both": len(present),
        "missing_from_matrix": len(missing),
        "only_in_matrix": len(only_in_matrix),
        "present_rate": (len(present) / max(1, len(set_genes)))
    }

    return {
        "present": present,
        "missing": missing,
        "only_in_matrix": only_in_matrix,
        "summary": summary,
    }


# ==========================================
# 4) CLI and main flow
# ==========================================

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Overlap report between an IRD gene list and a phylogenetic profile matrix (NPP/LPP)."
        )
    )
    # Inputs
    p.add_argument("--matrix", required=True, help="Path to matrix (CSV/TSV/Parquet/Feather)")
    p.add_argument("--genes", required=True, help="Path to gene list table (CSV/TSV/XLSX)")
    p.add_argument("--matrix-gene-col", default=None, help="Gene column name in the matrix (optional)")
    p.add_argument("--genes-col", default="Gene", help="Gene column name in the genes table (default: 'Gene')")
    p.add_argument("--matrix-delim", default=None, help="Delimiter override for CSV/TSV matrix")
    p.add_argument("--genes-delim", default=None, help="Delimiter override for CSV/TSV genes table")
    p.add_argument("--label", default="IRD", help="Label for outputs (e.g., 'IRD' or 'RetNet')")

    # Outputs
    p.add_argument("--outdir", default="output_data", help="Directory for all outputs")
    p.add_argument("--excel", action="store_true", help="Also write a 3-sheet Excel file")

    return p


def main() -> int:
    args = build_parser().parse_args()

    ensure_dir(args.outdir)
    run_ts = ts()
    base = f"{args.label}_overlap_{run_ts}"
    logfile = os.path.join(args.outdir, f"{base}.log")
    logger = setup_logger(logfile)

    logger.info("[START] IRD gene overlap report")

    # 4.1 Load inputs
    mtx = load_matrix(args.matrix, args.matrix_delim, logger)
    genes_tbl = load_gene_table(args.genes, args.genes_delim, logger)

    # 4.2 Detect columns
    matrix_gene_col = detect_col(mtx, args.matrix_gene_col, CANDIDATE_GENE_COLS)
    genes_col = detect_col(genes_tbl, args.genes_col, CANDIDATE_GENE_COLS)
    logger.info(f"[STRUCT] Matrix gene column: '{matrix_gene_col}' | Genes table column: '{genes_col}'")

    # 4.3 Compute overlap
    result = compute_overlap(mtx, matrix_gene_col, genes_tbl[genes_col])
    summary = result["summary"]

    logger.info(
        f"[STATS] user={summary['user_list_size']} | matrix={summary['matrix_index_size']} | "
        f"present={summary['present_in_both']} | missing={summary['missing_from_matrix']} | "
        f"only_in_matrix={summary['only_in_matrix']} | present_rate={summary['present_rate']:.3f}"
    )

    # 4.4 Write outputs
    present_df = pd.DataFrame({"Gene": result["present"]})
    missing_df = pd.DataFrame({"Gene": result["missing"]})
    only_df    = pd.DataFrame({"Gene": result["only_in_matrix"]})

    present_csv = os.path.join(args.outdir, f"{base}_present.csv")
    missing_csv = os.path.join(args.outdir, f"{base}_missing.csv")
    only_csv    = os.path.join(args.outdir, f"{base}_only_in_matrix.csv")
    summary_csv = os.path.join(args.outdir, f"{base}_summary.csv")
    summary_json= os.path.join(args.outdir, f"{base}_summary.json")

    logger.info("[SAVE] Writing CSV/JSON outputs…")
    present_df.to_csv(present_csv, index=False)
    missing_df.to_csv(missing_csv, index=False)
    only_df.to_csv(only_csv, index=False)
    pd.DataFrame([summary]).to_csv(summary_csv, index=False)
    pd.DataFrame([summary]).to_json(summary_json, orient="records", indent=2)

    if args.excel:
        xlsx_path = os.path.join(args.outdir, f"{base}_report.xlsx")
        logger.info("[SAVE] Writing 3-sheet Excel report…")
        try:
            with pd.ExcelWriter(xlsx_path, engine="openpyxl") as xw:
                present_df.to_excel(xw, sheet_name="Present_in_matrix", index=False)
                missing_df.to_excel(xw, sheet_name="Missing_from_matrix", index=False)
                only_df.to_excel(xw, sheet_name="Only_in_matrix", index=False)
        except Exception:
            # Fallback to xlsxwriter if openpyxl is missing
            with pd.ExcelWriter(xlsx_path, engine="xlsxwriter") as xw:
                present_df.to_excel(xw, sheet_name="Present_in_matrix", index=False)
                missing_df.to_excel(xw, sheet_name="Missing_from_matrix", index=False)
                only_df.to_excel(xw, sheet_name="Only_in_matrix", index=False)
        logger.info(f"[SAVE] Excel report: {xlsx_path}")

    logger.info("[END] Overlap report completed successfully ✅")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:
        print(f"[FATAL] {e}")
        sys.exit(2)
