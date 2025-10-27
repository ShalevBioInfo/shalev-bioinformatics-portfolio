#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# File: group_overlap_analysis.py
# Purpose: Compute overlaps between gene sets (groups) and produce overlap
#          matrices (counts and/or Jaccard) + summary CSVs for downstream
#          visualization and ranking of novel candidates.
# Created: 2025-10-27
# Author: Shalev Yaacov (refactored for GitHub publication)
#
# Short description:
#  - Input modes:
#     1) master table input with group column containing delimited gene lists
#        (--input + --group-col + --genes-col), OR
#     2) directory of one-file-per-group CSV/TSV (each file is list of genes)
#        (--groups-dir)
#  - Normalizes gene symbols (strip -> upper).
#  - Builds pairwise overlap matrix (counts of shared genes) and optional
#    Jaccard similarity matrix.
#  - Filters groups by minimum size if requested.
#  - Writes outputs: groups_list CSV, overlap_counts CSV, overlap_jaccard CSV
#    (if requested), and a small summary CSV. All outputs are timestamped and
#    stored in --outdir.
#
# Usage examples:
#  Mode A (master table):
#    python novel_overlap_analysis.py --input data/master_table.csv \
#      --group-col Inclusion_criterion --genes-col genes_list --sep '[;,|\\s]+' \
#      --method jaccard --outdir output_data
#
#  Mode B (groups dir):
#    python novel_overlap_analysis.py --groups-dir data/groups_csvs/ \
#      --method count --outdir output_data
#
# -----------------------------------------------------------------------------

# 1. Imports & environment ----------------------------------------------------
import argparse
import logging
import sys
from pathlib import Path
from datetime import datetime
from typing import Dict, Set, List, Tuple, Optional
import re

import pandas as pd
import numpy as np

# 2. Utilities: timestamp, dirs, logging -------------------------------------
def timestamp_now() -> str:
    """2.1 Compact timestamp for filenames."""
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def ensure_dir(path: Path):
    """2.2 Create directory if missing."""
    path.mkdir(parents=True, exist_ok=True)


def setup_logger(log_path: Path):
    """2.3 Setup logger to console and file (INFO)."""
    logger = logging.getLogger("novel_overlap")
    logger.setLevel(logging.INFO)
    if logger.hasHandlers():
        logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s", "%Y-%m-%d %H:%M:%S")
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)
    fh = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    fh.setLevel(logging.INFO)
    fh.setFormatter(fmt)
    logger.addHandler(ch)
    logger.addHandler(fh)
    return logger


# 3. I/O helpers --------------------------------------------------------------
def read_table(path: Path, sheet: Optional[str], logger) -> pd.DataFrame:
    """3.1 Read CSV/TSV/XLSX into DataFrame and trim column names."""
    logger.info(f"Reading table: {path}")
    if not path.exists():
        logger.error(f"Input not found: {path}")
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
    logger.info(f"Loaded table with shape {df.shape}")
    return df


# 4. Normalization & parsing --------------------------------------------------
def normalize_gene_name(g: str) -> str:
    """4.1 Normalize single gene name: strip and uppercase."""
    if pd.isna(g):
        return ""
    return str(g).strip().upper()


def split_gene_list(s: str, sep_pattern: str) -> List[str]:
    """4.2 Split a delimited gene list string by regex separators and normalize."""
    if s is None:
        return []
    if not isinstance(s, str):
        s = str(s)
    parts = re.split(sep_pattern, s)
    parts = [p.strip() for p in parts if p and p.strip() != ""]
    return parts


def build_groups_from_master(df: pd.DataFrame, group_col: str, genes_col: str,
                             sep_pattern: str, logger) -> Dict[str, Set[str]]:
    """
    4.3 From master table produce dict[group_name] -> set(gene_norm).
    If multiple rows share the same group value, gene lists are unioned.
    """
    logger.info("Building groups from master table.")
    if group_col not in df.columns:
        logger.error(f"Group column '{group_col}' not found in table.")
        raise KeyError(group_col)
    if genes_col not in df.columns:
        logger.error(f"Genes column '{genes_col}' not found in table.")
        raise KeyError(genes_col)

    groups: Dict[str, Set[str]] = {}
    for idx, row in df[[group_col, genes_col]].fillna("").iterrows():
        group = str(row[group_col]).strip()
        gene_list_raw = str(row[genes_col])
        if group == "" or gene_list_raw == "":
            continue
        items = split_gene_list(gene_list_raw, sep_pattern)
        normed = {normalize_gene_name(x) for x in items if normalize_gene_name(x) != ""}
        if not normed:
            continue
        if group not in groups:
            groups[group] = set()
        groups[group].update(normed)
    logger.info(f"Built {len(groups)} groups from table.")
    return groups


def build_groups_from_dir(dir_path: Path, logger) -> Dict[str, Set[str]]:
    """
    4.4 Read a directory of files; each file contains a list of genes (one per row
    or a single column). Returns dict[file_stem] -> set(gene_norm).
    """
    logger.info(f"Building groups from directory: {dir_path}")
    if not dir_path.exists():
        logger.error(f"Groups directory not found: {dir_path}")
        raise FileNotFoundError(dir_path)
    groups: Dict[str, Set[str]] = {}
    for p in sorted(dir_path.iterdir()):
        if p.is_file() and p.suffix.lower() in {".csv", ".tsv", ".txt", ".xlsx", ".xls"}:
            try:
                if p.suffix.lower() in {".csv", ".tsv", ".txt"}:
                    sep = "," if p.suffix.lower() == ".csv" else "\t"
                    df = pd.read_csv(p, sep=sep, dtype=str, header=None, encoding="utf-8", low_memory=False)
                    # assume first column contains gene names
                    col_series = df.iloc[:, 0].astype(str)
                else:
                    df = pd.read_excel(p, dtype=str)
                    # choose first non-empty column
                    col_series = df.iloc[:, 0].astype(str)
                gene_set = {normalize_gene_name(x) for x in col_series.fillna("").str.strip().tolist() if str(x).strip() != ""}
                groups[p.stem] = gene_set
                logger.info(f"Read group '{p.stem}' with {len(gene_set)} genes.")
            except Exception as e:
                logger.warning(f"Failed to read {p}: {e}")
    logger.info(f"Built {len(groups)} groups from directory.")
    return groups


# 5. Overlap computations -----------------------------------------------------
def compute_overlap_matrices(groups: Dict[str, Set[str]], method: str, logger) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    5.1 Compute two matrices:
      - counts_df: pairwise counts of shared genes between groups (symmetric)
      - jaccard_df: pairwise Jaccard similarities (intersection / union) if method includes 'jaccard'
    Returns (counts_df, jaccard_df) where jaccard_df may be None if not computed.
    """
    logger.info("Computing overlap matrices.")
    names = list(groups.keys())
    n = len(names)
    counts = np.zeros((n, n), dtype=int)
    # Precompute sizes
    sizes = [len(groups[name]) for name in names]
    for i in range(n):
        gi = groups[names[i]]
        for j in range(i, n):
            gj = groups[names[j]]
            inter = gi & gj
            c = len(inter)
            counts[i, j] = c
            counts[j, i] = c
    counts_df = pd.DataFrame(counts, index=names, columns=names)

    jaccard_df = None
    if method == "jaccard":
        logger.info("Computing Jaccard similarity matrix.")
        jacc = np.zeros((n, n), dtype=float)
        for i in range(n):
            gi = groups[names[i]]
            for j in range(i, n):
                gj = groups[names[j]]
                union = gi | gj
                if len(union) == 0:
                    val = 0.0
                else:
                    val = len(gi & gj) / len(union)
                jacc[i, j] = val
                jacc[j, i] = val
        jaccard_df = pd.DataFrame(jacc, index=names, columns=names)

    logger.info("Overlap matrices computed.")
    return counts_df, jaccard_df


# 6. Save outputs -------------------------------------------------------------
def save_groups_list(groups: Dict[str, Set[str]], outdir: Path, base_name: str, ts: str, logger) -> Path:
    """6.1 Save per-group gene lists as a CSV for reference."""
    ensure_dir(outdir)
    groups_df = pd.DataFrame({
        "group": list(groups.keys()),
        "size": [len(groups[g]) for g in groups.keys()],
        "genes": [";".join(sorted(groups[g])) for g in groups.keys()]
    })
    path = outdir / f"{base_name}_groups_list_{ts}.csv"
    logger.info(f"Saving groups list CSV: {path.name}")
    groups_df.to_csv(path, index=False, encoding="utf-8")
    return path


def save_matrices(counts_df: pd.DataFrame, jaccard_df: Optional[pd.DataFrame],
                  outdir: Path, base_name: str, ts: str, logger) -> Tuple[Path, Optional[Path]]:
    """6.2 Save counts and jaccard matrices as CSVs."""
    counts_path = outdir / f"{base_name}_overlap_counts_{ts}.csv"
    logger.info(f"Saving counts matrix CSV: {counts_path.name}")
    counts_df.to_csv(counts_path, index=True, encoding="utf-8")
    jacc_path = None
    if jaccard_df is not None:
        jacc_path = outdir / f"{base_name}_overlap_jaccard_{ts}.csv"
        logger.info(f"Saving Jaccard matrix CSV: {jacc_path.name}")
        jaccard_df.to_csv(jacc_path, index=True, encoding="utf-8")
    return counts_path, jacc_path


def save_summary(groups: Dict[str, Set[str]], counts_df: pd.DataFrame, outdir: Path, base_name: str, ts: str, logger) -> Path:
    """6.3 Save a simple summary CSV with group sizes and top overlaps."""
    rows = []
    for g in counts_df.index:
        sizes = len(groups[g])
        # find top 5 overlaps excluding self
        row = counts_df.loc[g].drop(labels=g).sort_values(ascending=False)
        top_overlap = ";".join([f"{idx}:{val}" for idx, val in row.head(5).items() if val > 0])
        rows.append({"group": g, "size": sizes, "top_overlaps": top_overlap})
    summary_df = pd.DataFrame(rows)
    summary_path = outdir / f"{base_name}_overlap_summary_{ts}.csv"
    logger.info(f"Saving summary CSV: {summary_path.name}")
    summary_df.to_csv(summary_path, index=False, encoding="utf-8")
    return summary_path


# 7. CLI ---------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Compute overlaps between gene groups (counts and/or Jaccard).")
    p.add_argument("--input", "-i", default=None, help="Master table CSV/XLSX with group column containing gene lists.")
    p.add_argument("--sheet", default=None, help="Sheet name/index if input is Excel.")
    p.add_argument("--group-col", default="group", help="Column name containing the group label in master table.")
    p.add_argument("--genes-col", default="genes", help="Column name containing gene lists in master table.")
    p.add_argument("--groups-dir", "-g", default=None, help="Directory with one file per group (each file: gene list).")
    p.add_argument("--sep", default=r"[;,\|\s]+", help="Regex separators for genes in list (default: '[;,,|\\s]+').")
    p.add_argument("--method", default="count", choices=["count", "jaccard"], help="Overlap metric to compute.")
    p.add_argument("--outdir", "-o", default="output_data", help="Output directory.")
    p.add_argument("--min-size", type=int, default=0, help="Minimum group size to keep (filter small groups).")
    return p.parse_args()


# 8. Main --------------------------------------------------------------------
def main():
    args = parse_args()
    ts = timestamp_now()
    outdir = Path(args.outdir)
    ensure_dir(outdir)
    log_path = outdir / f"novel_overlap_{ts}.log"
    logger = setup_logger(log_path)
    logger.info("=== START: novel_overlap_analysis ===")

    try:
        # 8.1 Build groups dict depending on mode
        groups: Dict[str, Set[str]] = {}
        if args.groups_dir:
            groups = build_groups_from_dir(Path(args.groups_dir), logger)
        elif args.input:
            df = read_table(Path(args.input), args.sheet, logger)
            groups = build_groups_from_master(df, args.group_col, args.genes_col, args.sep, logger)
        else:
            logger.error("Either --input or --groups-dir must be provided.")
            raise ValueError("Missing input source.")

        # 8.2 Filter by min-size
        if args.min_size > 0:
            before_n = len(groups)
            groups = {k: v for k, v in groups.items() if len(v) >= args.min_size}
            after_n = len(groups)
            logger.info(f"Filtered groups by min_size={args.min_size}: {before_n} -> {after_n}")

        if not groups:
            logger.error("No groups available after processing/filtering. Exiting.")
            return 1

        base_name = (Path(args.input).stem if args.input else Path(args.groups_dir).name)
        # 8.3 Save groups list CSV
        groups_list_path = save_groups_list(groups, outdir, base_name, ts, logger)

        # 8.4 Compute overlap matrices
        counts_df, jaccard_df = compute_overlap_matrices(groups, method=args.method, logger=logger)

        # 8.5 Save matrices and summary
        counts_path, jacc_path = save_matrices(counts_df, jaccard_df, outdir, base_name, ts, logger)
        summary_path = save_summary(groups, counts_df, outdir, base_name, ts, logger)

        # 8.6 Final summary log
        logger.info("=== SUMMARY ===")
        logger.info(f"Groups list: {groups_list_path.name}")
        logger.info(f"Counts matrix: {counts_path.name}")
        if jacc_path:
            logger.info(f"Jaccard matrix: {jacc_path.name}")
        logger.info(f"Summary: {summary_path.name}")
        logger.info("=== END: novel_overlap_analysis ===")
        return 0

    except Exception as e:
        logger.exception(f"Fatal error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
