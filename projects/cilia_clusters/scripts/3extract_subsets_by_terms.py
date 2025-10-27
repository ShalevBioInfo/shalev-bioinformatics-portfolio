#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# File: extract_subsets_by_terms.py
# Purpose: Search a master table for rows matching provided terms (plain substrings
#          or regex), extract per-term subsets, and save timestamped CSV/XLSX outputs
#          and a small summary report. Supports searching specific columns, regex,
#          and optional explosion of list-columns (e.g., "GENE1;GENE2").
# Created: 2025-10-27
# Author: Shalev Yaacov (refactored for GitHub publication)
#
# Short description:
#  - Read input CSV/XLSX table.
#  - Optionally explode one column that contains delimited lists into long form.
#  - For each search term (provided via CLI or file), find rows where term appears
#    in any of the selected columns (substring or regex, case-insensitive).
#  - Save per-term CSV and a combined Excel workbook (one sheet per term), plus a
#    summary CSV with counts and output paths.
#
# Example usage:
#  python groups_finder.py --input data/master_table.csv \
#      --outdir output_data --terms cilia,cilial --columns Cluster_ID,Description \
#      --explode-col cluster_genes --explode-sep "[,;|\\s]+" --regex False
#
# -----------------------------------------------------------------------------

# 1. Imports & environment ----------------------------------------------------
import argparse
import logging
import sys
from pathlib import Path
from datetime import datetime
from typing import List, Optional

import pandas as pd
import numpy as np
import re

# 2. Utilities: timestamp, logging, dirs -------------------------------------
def timestamp_now() -> str:
    """2.1 Return a compact timestamp for filenames."""
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def ensure_dir(path: Path):
    """2.2 Create directory if it doesn't exist."""
    path.mkdir(parents=True, exist_ok=True)


def setup_logger(log_path: Path):
    """2.3 Setup a logger that writes to stdout and to a file (INFO)."""
    logger = logging.getLogger("groups_finder")
    logger.setLevel(logging.INFO)
    if logger.hasHandlers():
        logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s", "%Y-%m-%d %H:%M:%S")
    # console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)
    # file handler
    fh = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    fh.setLevel(logging.INFO)
    fh.setFormatter(fmt)
    logger.addHandler(ch)
    logger.addHandler(fh)
    return logger


# 3. I/O helpers and normalization -------------------------------------------
def read_table(path: Path, sheet: Optional[str], logger) -> pd.DataFrame:
    """3.1 Read CSV/TSV/XLSX safely and trim column names."""
    logger.info(f"Reading input file: {path}")
    if not path.exists():
        logger.error(f"Input file not found: {path}")
        raise FileNotFoundError(path)
    sfx = path.suffix.lower()
    if sfx in {".csv", ".tsv", ".txt"}:
        sep = "," if sfx == ".csv" else "\t"
        df = pd.read_csv(path, sep=sep, dtype=str, encoding="utf-8", low_memory=False)
    elif sfx in {".xls", ".xlsx"}:
        df = pd.read_excel(path, sheet_name=sheet, dtype=str)
    else:
        df = pd.read_csv(path, dtype=str, encoding="utf-8", low_memory=False)
    df.columns = [str(c).strip() for c in df.columns]
    logger.info(f"Loaded table shape: {df.shape}")
    return df


def normalize_string_series(s: pd.Series) -> pd.Series:
    """3.2 Normalize strings: fillna, convert to str, strip whitespace."""
    return s.fillna("").astype(str).str.strip()


# 4. Explode helper ----------------------------------------------------------
def explode_column(df: pd.DataFrame, col: str, sep_pattern: str, logger) -> pd.DataFrame:
    """
    4.1 Explode a column that contains delimited lists into long form.
    - sep_pattern: regex for separators (e.g., r'[;,|\s]+').
    Returns new DataFrame with the exploded column named as original col.
    """
    logger.info(f"Exploding column '{col}' using separator pattern '{sep_pattern}'.")
    if col not in df.columns:
        logger.error(f"Explode column '{col}' not found in DataFrame.")
        raise KeyError(col)
    # Ensure column is string
    df = df.copy()
    df[col] = df[col].fillna("").astype(str)
    # Split into lists
    splitted = df[col].str.split(sep_pattern, expand=False)
    # Create long form by repeating index and exploding
    long_rows = []
    for idx, row in df.iterrows():
        items = splitted.at[idx] if idx in splitted.index else []
        if not items:
            # keep original row but with empty value to preserve context if desired
            long_rows.append((idx, ""))
        else:
            for it in items:
                it = it.strip()
                if it:
                    long_rows.append((idx, it))
    # Build exploded df by merging back with original other columns
    if not long_rows:
        logger.warning("Explode produced no items; returning original DataFrame.")
        return df
    idxs, values = zip(*long_rows)
    exploded_col = pd.Series(values, index=range(len(values)))
    # take original rows repeated
    df_rep = df.loc[list(idxs)].reset_index(drop=True)
    df_rep[col] = exploded_col.values
    logger.info(f"Exploded to {len(df_rep)} rows (from {len(df)} original rows).")
    return df_rep


# 5. Search mask builder -----------------------------------------------------
def build_search_mask(df: pd.DataFrame, terms: List[str], columns: List[str], use_regex: bool, logger) -> pd.Series:
    """
    5.1 Build boolean mask for rows matching any of the terms in any of the specified columns.
    - terms: list of search terms (strings).
    - columns: list of columns to search (subset of df.columns).
    - use_regex: treat terms as regex patterns (combined with OR) if True; otherwise substring search.
    Returns a boolean Series aligned with df.index.
    """
    logger.info(f"Building search mask: {len(terms)} terms, regex={use_regex}, columns={columns if columns else 'ALL'}")
    if not terms:
        logger.error("No search terms provided.")
        raise ValueError("No terms provided.")

    # ensure columns exist
    if columns:
        cols = [c for c in columns if c in df.columns]
        missing = [c for c in columns if c not in df.columns]
        if missing:
            logger.warning(f"Requested columns not found and will be ignored: {missing}")
        if not cols:
            logger.error("None of the requested columns exist in the DataFrame.")
            raise ValueError("No valid columns to search.")
    else:
        cols = df.columns.tolist()

    # Pre-normalize columns to strings
    for c in cols:
        df[c] = normalize_string_series(df[c])

    mask_total = pd.Series(False, index=df.index)

    if use_regex:
        # combine terms into one regex OR pattern (case-insensitive)
        # protect if user provided patterns already containing groups; simply join by |
        pattern = "|".join(f"(?:{t})" for t in terms)
        try:
            for c in cols:
                mask_col = df[c].str.contains(pattern, case=False, na=False, regex=True)
                mask_total = mask_total | mask_col
        except re.error as e:
            logger.error(f"Regex compilation error: {e}")
            raise
    else:
        # substring search: loop terms and OR accumulate (case-insensitive)
        for t in terms:
            if t == "":
                continue
            for c in cols:
                mask_col = df[c].str.contains(re.escape(t), case=False, na=False, regex=True)
                mask_total = mask_total | mask_col

    logger.info(f"Search mask: {mask_total.sum()} matching rows found.")
    return mask_total


# 6. Save outputs ------------------------------------------------------------
def save_term_outputs(term: str, term_df: pd.DataFrame, outdir: Path, base_name: str, ts: str, logger):
    """
    6.1 Save term-specific dataframe to CSV and also append to an Excel workbook.
    Returns (csv_path, sheet_name).
    """
    safe_term = re.sub(r"[^\w\-]+", "_", term)[:100]  # safe sheet/file name
    csv_path = outdir / f"{base_name}_group_{safe_term}_{ts}.csv"
    logger.info(f"Saving CSV for term '{term}': {csv_path.name}")
    term_df.to_csv(csv_path, index=False, encoding="utf-8")
    sheet_name = safe_term[:31]  # Excel sheet name limit
    return csv_path, sheet_name


# 7. CLI ---------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Find groups (row subsets) in a master table by search terms.")
    p.add_argument("--input", "-i", required=True, help="Path to input CSV/XLSX master table.")
    p.add_argument("--sheet", default=None, help="If input is Excel, sheet name or index.")
    p.add_argument("--outdir", "-o", default="output_data", help="Output directory.")
    p.add_argument("--terms", "-t", default=None,
                   help="Comma-separated terms (e.g. cilia,transport). Ignored if --terms-file is used alone.")
    p.add_argument("--terms-file", "-f", default=None,
                   help="Path to a plain text file with one term per line (comments starting with # ignored).")
    p.add_argument("--columns", "-c", default=None,
                   help="Comma-separated list of columns to search. If omitted, search all columns.")
    p.add_argument("--regex", action="store_true", help="Treat terms as regular expressions (case-insensitive).")
    p.add_argument("--explode-col", default=None,
                   help="Optional column name to explode (contains delimited lists, e.g. genes).")
    p.add_argument("--explode-sep", default=r"[,;|\s]+",
                   help="Regex pattern for separators when exploding (default: '[,;|\\s]+').")
    p.add_argument("--min-matches", type=int, default=0,
                   help="If >0, only save term outputs with at least this many matches.")
    return p.parse_args()


# 8. Main --------------------------------------------------------------------
def main():
    args = parse_args()
    input_path = Path(args.input)
    outdir = Path(args.outdir)
    ensure_dir(outdir)
    ts = timestamp_now()
    log_path = outdir / f"groups_finder_{ts}.log"
    logger = setup_logger(log_path)

    logger.info("=== START: groups_finder ===")
    logger.info(f"Input: {input_path}")
    logger.info(f"Outdir: {outdir}")
    logger.info(f"Regex mode: {args.regex}")

    try:
        # 8.1 Read input table
        df = read_table(input_path, args.sheet, logger)
        original_shape = df.shape

        # 8.2 Optionally explode a list-column into long form
        if args.explode_col:
            logger.info(f"Explode column requested: {args.explode_col}")
            df = explode_column(df, args.explode_col, args.explode_sep, logger)
            logger.info(f"Data shape after explode: {df.shape}")

        # 8.3 Build terms list (from CLI and/or file)
        terms: List[str] = []
        if args.terms_file:
            tf = Path(args.terms_file)
            if not tf.exists():
                logger.error(f"Terms file not found: {tf}")
                raise FileNotFoundError(tf)
            with tf.open("r", encoding="utf-8") as fh:
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    terms.append(line)
        if args.terms:
            terms_cli = [t.strip() for t in args.terms.split(",") if t.strip()]
            terms.extend(terms_cli)
        # Deduplicate while preserving order
        seen = set()
        terms = [x for x in terms if not (x in seen or seen.add(x))]
        if not terms:
            logger.error("No terms provided via --terms or --terms-file.")
            raise ValueError("No search terms supplied.")

        # 8.4 Determine columns to search
        columns = None
        if args.columns:
            columns = [c.strip() for c in args.columns.split(",") if c.strip()]
        logger.info(f"Searching columns: {columns if columns else 'ALL'}")

        # 8.5 Prepare for saving combined workbook
        base_name = input_path.stem
        workbook_path = outdir / f"{base_name}_groups_{ts}.xlsx"
        excel_writer = pd.ExcelWriter(workbook_path, engine="openpyxl")

        # 8.6 Iterate terms, build mask, save outputs, collect summary
        summary_rows = []
        sheets_written = []
        for term in terms:
            logger.info(f"--- Processing term: '{term}' ---")
            mask = build_search_mask(df, [term], columns, use_regex=args.regex, logger=logger)
            matched = df.loc[mask].copy()
            match_count = matched.shape[0]
            logger.info(f"Term '{term}': {match_count} rows matched.")

            # Possibly filter by min_matches
            if match_count < args.min_matches:
                logger.info(f"Skipping saving for term '{term}' (matches {match_count} < min_matches {args.min_matches}).")
                summary_rows.append({
                    "term": term,
                    "matched_rows": match_count,
                    "csv": "",
                    "sheet": "",
                })
                continue

            # Save CSV and add sheet to workbook
            csv_path, sheet_name = save_term_outputs(term, matched, outdir, base_name, ts, logger)
            try:
                matched.to_excel(excel_writer, sheet_name=sheet_name, index=False)
                sheets_written.append(sheet_name)
            except Exception as e:
                logger.warning(f"Could not write sheet '{sheet_name}' to workbook: {e}")

            summary_rows.append({
                "term": term,
                "matched_rows": match_count,
                "csv": str(csv_path.name),
                "sheet": sheet_name,
            })

        # finalize workbook
        try:
            excel_writer.save()
            logger.info(f"Combined workbook saved: {workbook_path.name} (sheets: {len(sheets_written)})")
        except Exception as e:
            logger.warning(f"Could not save combined workbook: {e}")

        # 8.7 Save summary CSV
        summary_df = pd.DataFrame(summary_rows)
        summary_csv = outdir / f"{base_name}_groups_summary_{ts}.csv"
        summary_df.to_csv(summary_csv, index=False, encoding="utf-8")
        logger.info(f"Summary saved: {summary_csv.name}")

        # 8.8 Final log and exit
        logger.info("=== SUMMARY ===")
        logger.info(f"Original shape: {original_shape}")
        logger.info(f"Terms processed: {len(terms)}")
        logger.info(f"Total saved sheets: {len(sheets_written)}")
        logger.info("=== END: groups_finder ===")
        return 0

    except Exception as e:
        logger.exception(f"Fatal error in groups_finder: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
