#!/usr/bin/env python3
# ====================================================================================
# Script: classify_gene_matrix_overlap.py
# Purpose: Publication‑ready overlap report between an IRD gene list and a profile
#          matrix (NPP/LPP), using a SYMBOL↔ENSG mapping to classify matches by
#          SYMBOL, by ENSG, or by BOTH, and optionally export a subset matrix.
# Author: Shalev Yaacov  |  refactored for publication
# Created: 2025-10-28
#
# Macro (what it does):
#   1) Loads matrix (CSV/TSV/Parquet/Feather) and detects the gene column.
#   2) Loads mapping with two columns (symbol, ensembl_gene_id) and normalizes case.
#   3) Loads IRD (or any) gene list and tries to match each gene via mapping to the
#      matrix index, classifying per‑gene as: SYMBOL / ENSG / BOTH / NONE.
#   4) Writes timestamped outputs under output_data/: per‑class CSVs, a summary CSV/JSON,
#      optional 4‑sheet Excel, and (optional) a subset matrix containing only matched genes.
#   5) Logs to console and to a .log file with the same timestamp suffix.
#
# Notes:
#   • Mapping is expected to have columns named (or detected as) 'symbol' and 'ensembl_gene_id'.
#   • The matrix gene column can be symbol or ENSG; detection is automatic, case‑insensitive.
# ====================================================================================

# ==========================================
# 0) Imports & constants
# ==========================================
import argparse
import logging
import os
import sys
from datetime import datetime
from typing import Optional, Sequence, Tuple

import pandas as pd

CANDIDATE_GENE_COLS = [
    "gene", "genes", "symbol", "gene_symbol", "hgnc_symbol",
    "Gene", "Genes", "SYMBOL", "GeneSymbol", "HGNC_symbol"
]
MAP_SYMBOL_COLS = ["symbol", "hgnc_symbol", "approved_symbol", "Symbol", "SYMBOL"]
MAP_ENSG_COLS   = ["ensembl_gene_id", "Ensembl Gene ID", "ENSG", "ensembl"]


# ==========================================
# 1) Utilities: time, io, logging, detection
# ==========================================

def ts() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def setup_logger(logfile: str) -> logging.Logger:
    logger = logging.getLogger("OverlapWithMapping")
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
        raise ValueError(f"Requested column '{user_col}' not found.")
    lower_map = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    # Fallback: first column
    return df.columns[0]


def load_table(path: str, delimiter: Optional[str] = None) -> pd.DataFrame:
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
# 2) Core logic
# ==========================================

def normalize_series(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.upper()


def classify_matches(
    matrix_df: pd.DataFrame,
    matrix_gene_col: str,
    mapping_df: pd.DataFrame,
    map_symbol_col: str,
    map_ensg_col: str,
    genes_df: pd.DataFrame,
    genes_col: str,
):
    # Normalize matrix index strings
    mtx_gene = normalize_series(matrix_df[matrix_gene_col])
    set_mtx = set(mtx_gene)

    # Normalize mapping
    mdf = mapping_df[[map_symbol_col, map_ensg_col]].dropna().copy()
    mdf.columns = ["symbol", "ensembl_gene_id"]
    mdf["symbol"] = normalize_series(mdf["symbol"])  # uppercase symbols
    mdf["ensembl_gene_id"] = normalize_series(mdf["ensembl_gene_id"])  # uppercase ENSG

    # Normalize user genes
    g = normalize_series(genes_df[genes_col]).rename("Gene")

    # Join mapping to user list (left join) to get candidate ENSG for each input symbol
    user_map = pd.merge(pd.DataFrame({"Gene": g}), mdf, how="left", left_on="Gene", right_on="symbol")

    # Precompute presence sets
    def is_in_matrix(values: pd.Series) -> pd.Series:
        return values.isin(mtx_gene)  # matrix values are uppercase too

    # Classification flags
    by_symbol = is_in_matrix(user_map["symbol"]).fillna(False)
    by_ensg   = is_in_matrix(user_map["ensembl_gene_id"]).fillna(False)

    # BOTH if both true, NONE if neither
    cls = (
        by_symbol.astype(int) + by_ensg.astype(int)
    ).map({2: "BOTH", 1: None, 0: "NONE"})

    # Replace 1 with specific type
    cls = cls.where(~(by_symbol ^ by_ensg), other=("SYMBOL" if by_symbol.any() else "ENSG"))

    # The above line sets one global label; we need per-row. Fix with vectorized choice:
    cls = (
        (by_symbol & by_ensg).map({True: "BOTH", False: None})
        .fillna((by_symbol & ~by_ensg).map({True: "SYMBOL", False: None}))
        .fillna((~by_symbol & by_ensg).map({True: "ENSG", False: None}))
        .fillna("NONE")
    )

    out = user_map[["Gene", "symbol", "ensembl_gene_id"]].copy()
    out["match_class"] = cls.values
    out["in_matrix_symbol"] = by_symbol.values
    out["in_matrix_ensg"] = by_ensg.values

    # Unique gene sets per class
    classes = {
        "SYMBOL": out.loc[out.match_class == "SYMBOL", "Gene"].drop_duplicates().tolist(),
        "ENSG":   out.loc[out.match_class == "ENSG", "Gene"].drop_duplicates().tolist(),
        "BOTH":   out.loc[out.match_class == "BOTH", "Gene"].drop_duplicates().tolist(),
        "NONE":   out.loc[out.match_class == "NONE", "Gene"].drop_duplicates().tolist(),
    }

    # Also return the augmented table for Excel/debug
    return out, classes


def subset_matrix(
    matrix_df: pd.DataFrame,
    matrix_gene_col: str,
    matches_tbl: pd.DataFrame,
) -> pd.DataFrame:
    # Collect all matrix keys matched by symbol or ensg
    keys = set(normalize_series(matches_tbl.loc[matches_tbl["in_matrix_symbol"], "symbol"]))
    keys |= set(normalize_series(matches_tbl.loc[matches_tbl["in_matrix_ensg"], "ensembl_gene_id"]))

    mtx_key_norm = normalize_series(matrix_df[matrix_gene_col])
    mask = mtx_key_norm.isin(keys)
    return matrix_df.loc[mask].copy()


# ==========================================
# 3) CLI & main
# ==========================================

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Overlap between IRD genes and a matrix using SYMBOL↔ENSG mapping; classify and optionally subset the matrix."
        )
    )
    # Inputs
    p.add_argument("--matrix", required=True, help="Path to matrix (CSV/TSV/Parquet/Feather)")
    p.add_argument("--mapping", required=True, help="Path to mapping with columns [symbol, ensembl_gene_id]")
    p.add_argument("--genes", required=True, help="Path to IRD (or other) gene list (CSV/TSV/XLSX)")

    # Column overrides
    p.add_argument("--matrix-gene-col", default=None, help="Gene column name in matrix (optional)")
    p.add_argument("--mapping-symbol-col", default=None, help="Symbol column in mapping (optional)")
    p.add_argument("--mapping-ensg-col", default=None, help="ENSG column in mapping (optional)")
    p.add_argument("--genes-col", default="Gene", help="Gene column in genes list (default: 'Gene')")

    # Delimiter overrides
    p.add_argument("--matrix-delim", default=None, help="Delimiter for CSV/TSV matrix")
    p.add_argument("--mapping-delim", default=None, help="Delimiter for CSV/TSV mapping")
    p.add_argument("--genes-delim", default=None, help="Delimiter for CSV/TSV genes")

    # Outputs
    p.add_argument("--label", default="IRD", help="Label prefix for outputs")
    p.add_argument("--outdir", default="output_data", help="Output directory")
    p.add_argument("--excel", action="store_true", help="Write 4‑sheet Excel (SYMBOL/ENSG/BOTH/NONE + augmented)")
    p.add_argument("--subset-matrix", action="store_true", help="Export subset matrix of matched rows")
    return p


def main() -> int:
    args = build_parser().parse_args()

    ensure_dir(args.outdir)
    run_ts = ts()
    base = f"{args.label}_overlap_mapping_{run_ts}"
    logfile = os.path.join(args.outdir, f"{base}.log")
    logger = setup_logger(logfile)

    logger.info("[START] Overlap with mapping")

    # 3.1 Load inputs
    mtx = load_table(args.matrix if args.matrix else "")
    mp  = load_table(args.mapping if args.mapping else "")
    gl  = load_table(args.genes if args.genes else "")

    # 3.2 Detect columns
    mtx_gene_col = detect_col(mtx, args.matrix_gene_col, CANDIDATE_GENE_COLS)
    map_sym_col  = detect_col(mp,  args.mapping_symbol_col, MAP_SYMBOL_COLS)
    map_ens_col  = detect_col(mp,  args.mapping_ensg_col,   MAP_ENSG_COLS)
    genes_col    = detect_col(gl,  args.genes_col,          CANDIDATE_GENE_COLS)
    logger.info(f"[STRUCT] Matrix gene col='{mtx_gene_col}' | Mapping: symbol='{map_sym_col}', ensg='{map_ens_col}' | Genes col='{genes_col}'")

    # 3.3 Classify matches
    aug_tbl, classes = classify_matches(
        matrix_df=mtx,
        matrix_gene_col=mtx_gene_col,
        mapping_df=mp,
        map_symbol_col=map_sym_col,
        map_ensg_col=map_ens_col,
        genes_df=gl,
        genes_col=genes_col,
    )

    n_user = gl[genes_col].dropna().astype(str).str.upper().nunique()
    stats = {
        "n_user_genes": int(n_user),
        "n_symbol": len(classes["SYMBOL"]),
        "n_ensg": len(classes["ENSG"]),
        "n_both": len(classes["BOTH"]),
        "n_none": len(classes["NONE"]),
        "present_rate": ( (len(classes["SYMBOL"]) + len(classes["ENSG"]) + len(classes["BOTH"]) ) / max(1, n_user) )
    }
    logger.info(
        f"[STATS] user={stats['n_user_genes']} | SYMBOL={stats['n_symbol']} | ENSG={stats['n_ensg']} | BOTH={stats['n_both']} | NONE={stats['n_none']} | present_rate={stats['present_rate']:.3f}"
    )

    # 3.4 Save per‑class CSVs and augmented table
    aug_csv   = os.path.join(args.outdir, f"{base}_augmented_table.csv")
    sym_csv   = os.path.join(args.outdir, f"{base}_SYMBOL.csv")
    ensg_csv  = os.path.join(args.outdir, f"{base}_ENSG.csv")
    both_csv  = os.path.join(args.outdir, f"{base}_BOTH.csv")
    none_csv  = os.path.join(args.outdir, f"{base}_NONE.csv")
    summary_csv = os.path.join(args.outdir, f"{base}_summary.csv")
    summary_json= os.path.join(args.outdir, f"{base}_summary.json")

    aug_tbl.to_csv(aug_csv, index=False)
    pd.DataFrame({"Gene": classes["SYMBOL"]}).to_csv(sym_csv, index=False)
    pd.DataFrame({"Gene": classes["ENSG"]}).to_csv(ensg_csv, index=False)
    pd.DataFrame({"Gene": classes["BOTH"]}).to_csv(both_csv, index=False)
    pd.DataFrame({"Gene": classes["NONE"]}).to_csv(none_csv, index=False)
    pd.DataFrame([stats]).to_csv(summary_csv, index=False)
    pd.DataFrame([stats]).to_json(summary_json, orient="records", indent=2)

    # 3.5 Optional Excel
    if args.excel:
        xlsx = os.path.join(args.outdir, f"{base}_report.xlsx")
        logger.info("[SAVE] Writing Excel report…")
        try:
            with pd.ExcelWriter(xlsx, engine="openpyxl") as xw:
                aug_tbl.to_excel(xw, sheet_name="Augmented", index=False)
                pd.DataFrame({"Gene": classes["SYMBOL"]}).to_excel(xw, sheet_name="SYMBOL", index=False)
                pd.DataFrame({"Gene": classes["ENSG"]}).to_excel(xw, sheet_name="ENSG", index=False)
                pd.DataFrame({"Gene": classes["BOTH"]}).to_excel(xw, sheet_name="BOTH", index=False)
                pd.DataFrame({"Gene": classes["NONE"]}).to_excel(xw, sheet_name="NONE", index=False)
        except Exception:
            with pd.ExcelWriter(xlsx, engine="xlsxwriter") as xw:
                aug_tbl.to_excel(xw, sheet_name="Augmented", index=False)
                pd.DataFrame({"Gene": classes["SYMBOL"]}).to_excel(xw, sheet_name="SYMBOL", index=False)
                pd.DataFrame({"Gene": classes["ENSG"]}).to_excel(xw, sheet_name="ENSG", index=False)
                pd.DataFrame({"Gene": classes["BOTH"]}).to_excel(xw, sheet_name="BOTH", index=False)
                pd.DataFrame({"Gene": classes["NONE"]}).to_excel(xw, sheet_name="NONE", index=False)
        logger.info(f"[SAVE] Excel: {xlsx}")

    # 3.6 Optional subset matrix
    if args.subset_matrix:
        sub = subset_matrix(mtx, mtx_gene_col, aug_tbl)
        sub_path = os.path.join(args.outdir, f"{base}_subset_matrix.csv")
        sub.to_csv(sub_path, index=False)
        logger.info(f"[SAVE] Subset matrix: {sub_path}")

    logger.info("[END] Overlap with mapping completed ✅")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:
        print(f"[FATAL] {e}")
        sys.exit(2)
