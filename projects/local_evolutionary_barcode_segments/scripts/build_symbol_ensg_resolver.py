#!/usr/bin/env python3
# ====================================================================================
# Script: build_symbol_ensg_resolver.py
# Purpose: Build a clean SYMBOL ↔ ENSG mapping from BioMart (optionally enriched with
#          HGNC synonyms) and generate a coverage report for a user gene list.
# Author: Shalev Yaacov  |  refactored for publication
# Created: 2025-10-28
#
# What this script does (macro):
#   1) Loads a BioMart export (CSV/TSV/XLSX/Parquet/Feather) and detects relevant columns.
#   2) Optionally merges HGNC synonyms (CSV/TSV/XLSX) and explodes alias/previous symbols.
#   3) Normalizes and deduplicates to produce a clean two-column mapping: symbol, ensembl_gene_id.
#   4) (Optional) Loads a user gene list and computes coverage (present/missing, multi-maps).
#   5) Writes timestamped outputs to output_data/: id_map.csv, summary.json/csv, problem lists
#      and an optional 3–5 sheet Excel report, plus a run log.
#
# Notes:
#   • Inputs/columns are configurable via CLI flags. Logging to console and to a .log file.
# ====================================================================================

# ==========================================
# 0) Imports & constants
# ==========================================
import argparse
import logging
import os
import sys
from datetime import datetime
from typing import Optional, Sequence, Tuple, List

import numpy as np
import pandas as pd

CANDIDATE_SYMBOL_COLS = [
    "symbol", "gene", "gene_name", "hgnc_symbol", "approved_symbol",
    "HGNC symbol", "HGNC_symbol", "Gene", "Symbol", "SYMBOL"
]
CANDIDATE_ENSG_COLS = [
    "ensembl_gene_id", "ensembl gene id", "ensembl_gene_id(s)",
    "ensembl_gene_id(s) [ensembl]", "ensembl", "ENSG"
]
CANDIDATE_ALIAS_COLS = [
    "alias_symbol", "alias symbols", "synonyms", "synonym", "aliases"
]
CANDIDATE_PREV_COLS = [
    "prev_symbol", "previous symbols", "previous_symbol"
]

# Delimiters used by various exports for multiple synonyms
MULTI_SEPARATORS = [";", ",", "|"]


# ==========================================
# 1) Utilities: time, io, logging, column detection
# ==========================================

def ts() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def setup_logger(logfile: str) -> logging.Logger:
    logger = logging.getLogger("IDMapBuilder")
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


def detect_col(df: pd.DataFrame, user_col: Optional[str], candidates: Sequence[str]) -> Optional[str]:
    """Detect a column name with case-insensitive matching; return None if not found."""
    if user_col:
        if user_col in df.columns:
            return user_col
        for c in df.columns:
            if c.lower() == user_col.lower():
                return c
        return None
    lower_map = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    return None


def load_table(path: str, delimiter: Optional[str] = None) -> pd.DataFrame:
    ext = os.path.splitext(path)[1].lower()
    if ext == ".parquet":
        return pd.read_parquet(path)
    if ext in (".feather", ".ft"):
        return pd.read_feather(path)
    if ext in (".xlsx", ".xls"):
        return pd.read_excel(path)
    # CSV/TSV
    if delimiter is None:
        delimiter = "\t" if ext in (".tsv", ".tab", ".txt") else ","
    return pd.read_csv(path, sep=delimiter, low_memory=False)


def explode_multi_series(s: pd.Series) -> pd.Series:
    """Explode strings by common separators into one symbol per row."""
    out = s.astype(str)
    for sep in MULTI_SEPARATORS:
        out = out.str.replace(f"{sep}", ";", regex=False)
    out = out.str.split(";")
    return out.explode().astype(str).str.strip()


# ==========================================
# 2) Build mapping from BioMart (+ optional HGNC)
# ==========================================

def build_id_map(
    biomart_df: pd.DataFrame,
    symbol_col: Optional[str] = None,
    ensg_col: Optional[str] = None,
    alias_df: Optional[pd.DataFrame] = None,
    alias_symbol_col: Optional[str] = None,
    alias_ensg_col: Optional[str] = None,
    alias_alias_col: Optional[str] = None,
    alias_prev_col: Optional[str] = None,
    logger: Optional[logging.Logger] = None,
) -> pd.DataFrame:
    # 2.1 Detect columns in BioMart
    sym = detect_col(biomart_df, symbol_col, CANDIDATE_SYMBOL_COLS)
    ensg = detect_col(biomart_df, ensg_col, CANDIDATE_ENSG_COLS)
    if sym is None or ensg is None:
        raise ValueError("Could not detect symbol/Ensembl columns in the BioMart table. Use --biomart-symbol/--biomart-ensg.")

    if logger: logger.info(f"[Biomart] Using columns: symbol='{sym}', ensembl='{ensg}'")

    base = biomart_df[[sym, ensg]].copy()
    base.columns = ["symbol", "ensembl_gene_id"]

    # 2.2 Optional HGNC alias/previous symbols
    extra = pd.DataFrame(columns=["symbol", "ensembl_gene_id"])  # empty default
    if alias_df is not None and not alias_df.empty:
        a_sym = detect_col(alias_df, alias_symbol_col, CANDIDATE_SYMBOL_COLS) or sym
        a_ens = detect_col(alias_df, alias_ensg_col, CANDIDATE_ENSG_COLS) or ensg
        a_alias = detect_col(alias_df, alias_alias_col, CANDIDATE_ALIAS_COLS)
        a_prev  = detect_col(alias_df, alias_prev_col,  CANDIDATE_PREV_COLS)
        if logger: logger.info(
            f"[HGNC] Using columns: symbol='{a_sym}', ensembl='{a_ens}', alias='{a_alias}', previous='{a_prev}'")

        cols_to_keep = {a_sym: "symbol", a_ens: "ensembl_gene_id"}
        extra_blocks: List[pd.DataFrame] = []

        # aliases
        if a_alias and a_ens in alias_df.columns:
            tmp = alias_df[[a_alias, a_ens]].dropna()
            tmp = tmp.assign(symbol=explode_multi_series(tmp[a_alias])).drop(columns=[a_alias])
            tmp = tmp.rename(columns={a_ens: "ensembl_gene_id"})
            extra_blocks.append(tmp[["symbol", "ensembl_gene_id"]])

        # previous symbols
        if a_prev and a_ens in alias_df.columns:
            tmp = alias_df[[a_prev, a_ens]].dropna()
            tmp = tmp.assign(symbol=explode_multi_series(tmp[a_prev])).drop(columns=[a_prev])
            tmp = tmp.rename(columns={a_ens: "ensembl_gene_id"})
            extra_blocks.append(tmp[["symbol", "ensembl_gene_id"]])

        # approved symbols in HGNC file (if present)
        if a_sym in alias_df.columns and a_ens in alias_df.columns:
            tmp = alias_df[[a_sym, a_ens]].dropna()
            tmp = tmp.rename(columns={a_sym: "symbol", a_ens: "ensembl_gene_id"})
            extra_blocks.append(tmp[["symbol", "ensembl_gene_id"]])

        if extra_blocks:
            extra = pd.concat(extra_blocks, ignore_index=True)

    # 2.3 Concatenate and clean
    mapping = pd.concat([base, extra], ignore_index=True)
    mapping["symbol"] = mapping["symbol"].astype(str).str.strip()
    mapping["ensembl_gene_id"] = mapping["ensembl_gene_id"].astype(str).str.strip()

    # Drop empty
    mapping = mapping[(mapping["symbol"] != "") & (mapping["ensembl_gene_id"] != "")]

    # Normalize: HGNC symbols are typically uppercase; ENSG uppercase
    mapping["symbol"] = mapping["symbol"].str.upper()
    mapping["ensembl_gene_id"] = mapping["ensembl_gene_id"].str.upper()

    # Remove obvious non-ENSG ids if present (keep strings starting with ENSG)
    mask_ensg = mapping["ensembl_gene_id"].str.startswith("ENSG", na=False)
    if mask_ensg.any():
        mapping = mapping[mask_ensg]

    # Drop duplicates (exact pairs)
    mapping = mapping.drop_duplicates(ignore_index=True)

    # 2.4 Flag multi-maps for the report stage
    # Note: do not drop; just keep info for coverage
    return mapping


# ==========================================
# 3) Coverage report (optional)
# ==========================================

def compute_coverage(mapping: pd.DataFrame, genes_df: pd.DataFrame, genes_col: Optional[str]) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    # Detect gene column in user list
    gcol = detect_col(genes_df, genes_col, CANDIDATE_SYMBOL_COLS) or (genes_df.columns[0])
    genes = genes_df[gcol].astype(str).str.strip().str.upper()

    set_genes = set(genes)
    set_map_symbols = set(mapping["symbol"])

    present = sorted(set_genes & set_map_symbols)
    missing = sorted(set_genes - set_map_symbols)

    present_df = pd.DataFrame({"Gene": present})
    missing_df = pd.DataFrame({"Gene": missing})

    # multi-map stats
    dup_symbol = mapping.groupby("symbol")["ensembl_gene_id"].nunique().reset_index(name="n_ensembl")
    dup_ensg   = mapping.groupby("ensembl_gene_id")["symbol"].nunique().reset_index(name="n_symbol")

    summary = {
        "n_input_genes": len(set_genes),
        "n_mapping_pairs": int(len(mapping)),
        "n_present": len(present),
        "n_missing": len(missing),
        "n_symbols_multi_map": int((dup_symbol["n_ensembl"]>1).sum()),
        "n_ensg_multi_map": int((dup_ensg["n_symbol"]>1).sum()),
        "present_rate": (len(present) / max(1, len(set_genes)))
    }

    return present_df, missing_df, dup_symbol, dup_ensg, summary


# ==========================================
# 4) CLI & main
# ==========================================

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Build SYMBOL↔ENSG mapping from BioMart (+HGNC synonyms) and compute coverage for a user list."
        )
    )
    # Inputs
    p.add_argument("--biomart", required=True, help="Path to BioMart export (CSV/TSV/XLSX/Parquet/Feather)")
    p.add_argument("--hgnc", default=None, help="Path to optional HGNC synonyms file (CSV/TSV/XLSX)")
    p.add_argument("--genes", default=None, help="Path to user gene list for coverage (CSV/TSV/XLSX)")

    # Column overrides
    p.add_argument("--biomart-symbol", default=None, help="Override symbol column name in BioMart table")
    p.add_argument("--biomart-ensg", default=None, help="Override Ensembl column name in BioMart table")
    p.add_argument("--hgnc-symbol", default=None, help="Override symbol column in HGNC file")
    p.add_argument("--hgnc-ensg", default=None, help="Override Ensembl column in HGNC file")
    p.add_argument("--hgnc-alias", default=None, help="Override alias/synonyms column in HGNC file")
    p.add_argument("--hgnc-prev", default=None, help="Override previous symbols column in HGNC file")
    p.add_argument("--genes-col", default=None, help="Override gene column in user list")

    # Delimiter overrides for CSV/TSV
    p.add_argument("--biomart-delim", default=None, help="Delimiter for BioMart CSV/TSV")
    p.add_argument("--hgnc-delim", default=None, help="Delimiter for HGNC CSV/TSV")
    p.add_argument("--genes-delim", default=None, help="Delimiter for user list CSV/TSV")

    # Outputs
    p.add_argument("--label", default="IDMAP", help="Label prefix for output files")
    p.add_argument("--outdir", default="output_data", help="Output directory")
    p.add_argument("--excel", action="store_true", help="Write a multi-sheet Excel report as well")
    return p


def main() -> int:
    args = build_parser().parse_args()

    ensure_dir(args.outdir)
    run_ts = ts()
    base = f"{args.label}_idmap_{run_ts}"
    logfile = os.path.join(args.outdir, f"{base}.log")
    logger = setup_logger(logfile)

    logger.info("[START] Build SYMBOL↔ENSG mapping")

    # 4.1 Load tables
    if not os.path.exists(args.biomart):
        logger.error(f"BioMart file not found: {args.biomart}")
        return 2
    biomart_df = load_table(args.biomart, args.biomart_delim)
    logger.info(f"[LOAD] BioMart shape: {biomart_df.shape}")

    alias_df = None
    if args.hgnc:
        if not os.path.exists(args.hgnc):
            logger.error(f"HGNC file not found: {args.hgnc}")
            return 2
        alias_df = load_table(args.hgnc, args.hgnc_delim)
        logger.info(f"[LOAD] HGNC shape: {alias_df.shape}")

    # 4.2 Build mapping
    mapping = build_id_map(
        biomart_df,
        symbol_col=args.biomart_symbol,
        ensg_col=args.biomart_ensg,
        alias_df=alias_df,
        alias_symbol_col=args.hgnc_symbol,
        alias_ensg_col=args.hgnc_ensg,
        alias_alias_col=args.hgnc_alias,
        alias_prev_col=args.hgnc_prev,
        logger=logger,
    )
    logger.info(f"[MAP] Clean mapping rows: {len(mapping)}")

    # 4.3 Save mapping
    map_csv = os.path.join(args.outdir, f"{base}_id_map.csv")
    mapping.to_csv(map_csv, index=False)
    logger.info(f"[SAVE] Mapping CSV: {map_csv}")

    # 4.4 Optional coverage
    present_df = missing_df = dup_symbol = dup_ensg = None
    summary = {}
    if args.genes:
        if not os.path.exists(args.genes):
            logger.error(f"Genes file not found: {args.genes}")
            return 2
        genes_df = load_table(args.genes, args.genes_delim)
        logger.info(f"[LOAD] Genes list shape: {genes_df.shape}")
        present_df, missing_df, dup_symbol, dup_ensg, summary = compute_coverage(mapping, genes_df, args.genes_col)

        # Save coverage artifacts
        present_csv = os.path.join(args.outdir, f"{base}_present.csv")
        missing_csv = os.path.join(args.outdir, f"{base}_missing.csv")
        dup_sym_csv = os.path.join(args.outdir, f"{base}_symbols_multi_map.csv")
        dup_ensg_csv= os.path.join(args.outdir, f"{base}_ensg_multi_map.csv")
        summary_csv = os.path.join(args.outdir, f"{base}_summary.csv")
        summary_json= os.path.join(args.outdir, f"{base}_summary.json")

        present_df.to_csv(present_csv, index=False)
        missing_df.to_csv(missing_csv, index=False)
        dup_symbol.to_csv(dup_sym_csv, index=False)
        dup_ensg.to_csv(dup_ensg_csv, index=False)
        pd.DataFrame([summary]).to_csv(summary_csv, index=False)
        pd.DataFrame([summary]).to_json(summary_json, orient="records", indent=2)
        logger.info("[SAVE] Coverage CSV/JSON artifacts written")

        if args.excel:
            xlsx_path = os.path.join(args.outdir, f"{base}_report.xlsx")
            logger.info("[SAVE] Writing Excel report…")
            try:
                with pd.ExcelWriter(xlsx_path, engine="openpyxl") as xw:
                    mapping.to_excel(xw, sheet_name="ID_Map", index=False)
                    present_df.to_excel(xw, sheet_name="Present", index=False)
                    missing_df.to_excel(xw, sheet_name="Missing", index=False)
                    dup_symbol.to_excel(xw, sheet_name="Symbols_multi_map", index=False)
                    dup_ensg.to_excel(xw, sheet_name="ENSG_multi_map", index=False)
            except Exception:
                with pd.ExcelWriter(xlsx_path, engine="xlsxwriter") as xw:
                    mapping.to_excel(xw, sheet_name="ID_Map", index=False)
                    present_df.to_excel(xw, sheet_name="Present", index=False)
                    missing_df.to_excel(xw, sheet_name="Missing", index=False)
                    dup_symbol.to_excel(xw, sheet_name="Symbols_multi_map", index=False)
                    dup_ensg.to_excel(xw, sheet_name="ENSG_multi_map", index=False)
            logger.info(f"[SAVE] Excel: {xlsx_path}")

    # 4.5 Final log
    if summary:
        logger.info(
            f"[STATS] input_genes={summary['n_input_genes']} present={summary['n_present']} missing={summary['n_missing']} "
            f"multi_symbol={summary['n_symbols_multi_map']} multi_ensg={summary['n_ensg_multi_map']} present_rate={summary['present_rate']:.3f}"
        )
    logger.info("[END] ID map build completed ✅")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:
        print(f"[FATAL] {e}")
        sys.exit(2)
