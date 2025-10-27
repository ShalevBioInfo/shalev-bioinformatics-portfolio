#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# File: extract_cilia_annotations.py
# Purpose: Ingest and perform an initial filter for "cilia"-related rows
#          from a master table (CSV/XLSX). Produces cleaned master file and
#          filtered subset files (CSV + XLSX) in output_data/ with timestamps.
# Created: 2025-10-27
# Author: Shalev Yaacov (refactored for publication)
#
# Short description:
#   - Reads a CSV or Excel file containing cluster/gene annotations.
#   - Normalizes string columns and searches for one or more search terms
#     (default: "cilia") across either all columns or a user-specified subset.
#   - Exports: cleaned master (strings normalized), filtered subset (rows
#     matching any search term) as CSV and XLSX, plus a run log.
#
# Example usage:
#   python step1_extract_cilia.py \
#       --input data/all_343_clusters.xlsx \
#       --sheet "343 clusters - ordered by clust" \
#       --terms cilia,cilial \
#       --outdir output_data \
#       --columns Cluster_ID,Description
#
# -----------------------------------------------------------------------------

import argparse
import logging
import sys
from pathlib import Path
from datetime import datetime
from typing import List

import pandas as pd
import numpy as np

# -------------------- 1. Utilities & Logging setup ---------------------------
def timestamp_now() -> str:
    """Return compact timestamp for filenames."""
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def ensure_dir(path: Path):
    """Create directory if it does not exist."""
    path.mkdir(parents=True, exist_ok=True)


def setup_logger(log_path: Path):
    """
    Setup logger that prints to stdout and writes to a logfile (INFO).
    Returns the configured logger.
    """
    logger = logging.getLogger("step1_extract_cilia")
    logger.setLevel(logging.INFO)
    # Clear existing handlers to avoid duplicates in long sessions
    if logger.hasHandlers():
        logger.handlers.clear()

    fmt = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S")

    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.INFO)
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    fh = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    fh.setLevel(logging.INFO)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    return logger


# -------------------- 2. Argument parsing -----------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Step 1: Extract rows that mention 'cilia' (or other terms) from a master table."
    )
    p.add_argument("--input", "-i", required=True, help="Path to input CSV or XLSX table.")
    p.add_argument("--sheet", "-s", default=None, help="If input is Excel, the sheet name (or index).")
    p.add_argument("--outdir", "-o", default="output_data", help="Directory for outputs and logs.")
    p.add_argument("--terms", "-t", default="cilia",
                   help="Comma-separated search terms (default: 'cilia'). Terms are treated as plain substrings unless --regex is set.")
    p.add_argument("--columns", "-c", default=None,
                   help="Comma-separated list of columns to search. If omitted, search all columns.")
    p.add_argument("--regex", action="store_true",
                   help="Treat terms as regular expressions (case-insensitive).")
    p.add_argument("--min-rows", type=int, default=0,
                   help="If >0, require the input to have at least this many rows (sanity check).")
    return p.parse_args()


# -------------------- 3. I/O helpers & normalization ------------------------
def read_table(input_path: Path, sheet, logger) -> pd.DataFrame:
    """Read CSV/TSV/XLSX into pandas DataFrame with safe defaults."""
    logger.info(f"Reading input: {input_path}")
    if not input_path.exists():
        logger.error(f"Input file not found: {input_path}")
        raise FileNotFoundError(input_path)

    suffix = input_path.suffix.lower()
    if suffix in (".csv", ".tsv", ".txt"):
        sep = "," if suffix == ".csv" else "\t"
        df = pd.read_csv(input_path, sep=sep, dtype=str, encoding="utf-8", low_memory=False)
    elif suffix in (".xls", ".xlsx"):
        df = pd.read_excel(input_path, sheet_name=sheet, dtype=str)
    else:
        # fallback to CSV
        df = pd.read_csv(input_path, dtype=str, encoding="utf-8", low_memory=False)

    # Trim column names
    df.columns = [str(c).strip() for c in df.columns]
    logger.info(f"Loaded table with shape: {df.shape}")
    return df


def normalize_strings_df(df: pd.DataFrame, logger) -> pd.DataFrame:
    """
    Return a copy of df where all object/string columns are stripped of
    leading/trailing whitespace and None/NaN are replaced by empty string.
    """
    logger.info("Normalizing string columns (strip whitespace, fill NaN -> '').")
    df_copy = df.copy()
    obj_cols = df_copy.select_dtypes(include=["object"]).columns
    for col in obj_cols:
        # Convert to str and strip
        df_copy[col] = df_copy[col].fillna("").astype(str).str.strip()
    return df_copy


# -------------------- 4. Core filter logic ---------------------------------
def build_search_mask(df: pd.DataFrame, terms: List[str], columns: List[str], use_regex: bool, logger) -> pd.Series:
    """
    Build boolean mask of rows where any of the terms appears (in any of the given columns).
    - terms: list of strings (already split).
    - columns: columns to search (subset of df.columns).
    - use_regex: whether to treat terms as regex (case-insensitive).
    Returns pandas Series of booleans indexed like df.
    """
    logger.info(f"Building search mask for terms={terms} (regex={use_regex}).")
    if columns is None:
        cols_to_search = df.columns.tolist()
    else:
        # keep only existing columns
        cols_to_search = [c for c in columns if c in df.columns]
        missing = [c for c in columns if c not in df.columns]
        if missing:
            logger.warning(f"Requested columns not found and will be ignored: {missing}")
        if not cols_to_search:
            logger.error("No valid columns to search after filtering; aborting.")
            raise ValueError("No valid columns to search.")

    # Build mask per column then combine with OR
    mask_total = pd.Series(False, index=df.index)
    for col in cols_to_search:
        col_series = df[col].astype(str).fillna("")
        if use_regex:
            # combine multiple regex terms into one big regex OR
            pattern = "|".join(f"(?:{t})" for t in terms)
            try:
                mask_col = col_series.str.contains(pattern, case=False, na=False, regex=True)
            except Exception as e:
                logger.error(f"Regex error for column '{col}': {e}")
                raise
        else:
            # plain substring search across terms
            mask_col = pd.Series(False, index=df.index)
            for t in terms:
                mask_col = mask_col | col_series.str.contains(t, case=False, na=False)
        mask_total = mask_total | mask_col
    logger.info(f"Search mask created: {mask_total.sum()} matching rows out of {len(df)}.")
    return mask_total


# -------------------- 5. Save outputs --------------------------------------
def save_outputs(clean_df: pd.DataFrame, filtered_df: pd.DataFrame, outdir: Path, base_name: str, ts: str, logger):
    """
    Save cleaned master and filtered subset to CSV and XLSX with timestamped names.
    Returns paths to saved files.
    """
    ensure_dir(outdir)
    cleaned_csv = outdir / f"{base_name}_cleaned_{ts}.csv"
    filtered_csv = outdir / f"{base_name}_filtered_{ts}.csv"
    filtered_xlsx = outdir / f"{base_name}_filtered_{ts}.xlsx"

    logger.info(f"Saving cleaned master CSV: {cleaned_csv.name}")
    clean_df.to_csv(cleaned_csv, index=False, encoding="utf-8")

    logger.info(f"Saving filtered CSV: {filtered_csv.name}")
    filtered_df.to_csv(filtered_csv, index=False, encoding="utf-8")

    # Save XLSX for easier viewing (if pandas supports)
    try:
        logger.info(f"Saving filtered XLSX: {filtered_xlsx.name}")
        filtered_df.to_excel(filtered_xlsx, index=False)
    except Exception as e:
        logger.warning(f"Could not save XLSX file ({e}); CSV is available.")

    return cleaned_csv, filtered_csv, filtered_xlsx


# -------------------- 6. Main ------------------------------------------------
def main():
    args = parse_args()
    input_path = Path(args.input)
    outdir = Path(args.outdir)
    ensure_dir(outdir)
    ts = timestamp_now()
    log_path = outdir / f"step1_extract_cilia_{ts}.log"
    logger = setup_logger(log_path)

    logger.info("=== START: step1_extract_cilia ===")
    logger.info(f"Input: {input_path}")
    logger.info(f"Outdir: {outdir}")

    try:
        df = read_table(input_path, args.sheet, logger)

        if args.min_rows > 0 and len(df) < args.min_rows:
            logger.error(f"Input has {len(df)} rows which is less than required min_rows={args.min_rows}. Aborting.")
            raise SystemExit(2)

        # Normalize strings
        df_clean = normalize_strings_df(df, logger)

        # Columns to search
        columns = None
        if args.columns:
            columns = [c.strip() for c in args.columns.split(",") if c.strip()]

        # Terms list
        terms = [t.strip() for t in args.terms.split(",") if t.strip()]
        if not terms:
            logger.error("No search terms provided.")
            raise ValueError("No search terms provided.")

        # Build mask and filter
        mask = build_search_mask(df_clean, terms=terms, columns=columns, use_regex=args.regex, logger=logger)
        filtered = df_clean[mask].copy()

        # Save outputs
        base_name = input_path.stem
        cleaned_csv, filtered_csv, filtered_xlsx = save_outputs(df_clean, filtered, outdir, base_name, ts, logger)

        # Final summary
        logger.info("=== SUMMARY ===")
        logger.info(f"Total input rows: {len(df_clean)}")
        logger.info(f"Filtered rows matching terms ({', '.join(terms)}): {len(filtered)}")
        logger.info(f"Cleaned master: {cleaned_csv.name}")
        logger.info(f"Filtered CSV: {filtered_csv.name}")
        logger.info(f"Filtered XLSX: {filtered_xlsx.name} (if saved)")
        logger.info("=== END ===")
        return 0

    except Exception as e:
        logger.exception(f"Fatal error: {e}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
