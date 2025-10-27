#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# File: merge_clusters_genes_with_reference.py
# Purpose: Merge cluster gene lists (clusters file) with an external reference
#          table (CiliaCarta) by gene name. Produces exploded merged table
#          (one row per cluster:gene), a cluster-level summary (matched genes,
#          counts, unmatched genes), and an Excel workbook with provenance.
# Created: 2025-10-27
# Author: Shalev Yaacov (refactored for GitHub publication)
#
# Short description:
#  - Reads clusters table (CSV/XLSX) and CiliaCarta table (CSV/XLSX).
#  - Normalizes gene name strings and explodes cluster gene lists into one
#    gene per row (configurable separators).
#  - Joins exploded cluster genes with CiliaCarta on normalized gene name
#    (case-insensitive exact match).
#  - Aggregates back to cluster-level with counts and lists of matched/unmatched genes.
#  - Saves outputs (exploded merged CSV, cluster summary CSV, combined Excel workbook)
#    into output_data/ with timestamped filenames and run log.
#
# Notes:
#  - Matching is exact on normalized names (strip + upper). If you need fuzzy
#    synonyms mapping, provide a mapping beforehand and join on mapped names.
#  - Dependencies: pandas, numpy, openpyxl (for XLSX writing). Keep environment minimal.
#
# -----------------------------------------------------------------------------

# 1) Imports and environment ---------------------------------------------------
import argparse
import logging
import sys
from pathlib import Path
from datetime import datetime
from typing import List, Optional
import re

import pandas as pd
import numpy as np

# 2) Utilities: timestamp, logging, dirs -------------------------------------
def timestamp_now() -> str:
    """Return compact timestamp for filenames."""
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def ensure_dir(p: Path):
    """Create directory if missing."""
    p.mkdir(parents=True, exist_ok=True)


def setup_logger(log_path: Path):
    """2.1 Setup logger to stdout and file (INFO)."""
    logger = logging.getLogger("merge_clusters")
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


# 3) I/O helpers --------------------------------------------------------------
def read_table(path: Path, sheet: Optional[str], logger) -> pd.DataFrame:
    """3.1 Read CSV/TSV/XLSX safely and return DataFrame with trimmed column names."""
    logger.info(f"Reading input file: {path}")
    if not path.exists():
        logger.error(f"Input not found: {path}")
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
    logger.info(f"Loaded table shape: {df.shape}")
    return df


# 4) Normalization & explode logic -------------------------------------------
def normalize_gene_name(s: Optional[str]) -> str:
    """4.1 Normalize single gene name to canonical comparison form (strip + upper)."""
    if s is None:
        return ""
    return str(s).strip().upper()


def split_genes(s: str, sep_pattern: str) -> List[str]:
    """4.2 Split a gene-list string into individual gene names using regex sep."""
    if not isinstance(s, str) or s.strip() == "":
        return []
    parts = re.split(sep_pattern, s)
    # remove empties, strip spaces
    return [p.strip() for p in parts if p and p.strip()]


def explode_cluster_genes(df: pd.DataFrame, genes_col: str, cluster_id_col: str, sep_pattern: str, logger) -> pd.DataFrame:
    """
    4.3 Explode cluster gene lists into long form:
      Input: df with one row per cluster, and a 'genes' column (e.g. "GENE1;GENE2")
      Output: DataFrame with columns [cluster_id_col, gene_raw, gene_norm, source_row_index]
    """
    logger.info("Exploding cluster gene lists into one gene per row.")
    rows = []
    for idx, row in df[[cluster_id_col, genes_col]].fillna("").iterrows():
        cluster_id = str(row[cluster_id_col]).strip()
        gene_list_raw = row[genes_col]
        if gene_list_raw is None or str(gene_list_raw).strip() == "":
            continue
        genes = split_genes(str(gene_list_raw), sep_pattern)
        for g in genes:
            g_norm = normalize_gene_name(g)
            if g_norm == "":
                continue
            rows.append({cluster_id_col: cluster_id, "gene_raw": g, "gene_norm": g_norm, "_source_index": idx})
    long = pd.DataFrame(rows)
    logger.info(f"Exploded into {len(long)} cluster:gene rows (from {df.shape[0]} cluster rows).")
    return long


# 5) Merge logic -------------------------------------------------------------
def merge_with_reference(exploded_df: pd.DataFrame, ref_df: pd.DataFrame, ref_gene_col: str, logger) -> pd.DataFrame:
    """
    5.1 Join exploded cluster genes with the reference table on normalized gene name.
       Adds ref columns prefixed with 'ref_' to the exploded_df.
    """
    logger.info("Preparing reference table for join (normalize gene column).")
    ref = ref_df.copy()
    ref_cols = [c for c in ref.columns]
    # Add normalized gene column to reference
    ref["_ref_gene_norm"] = ref[ref_gene_col].fillna("").astype(str).apply(normalize_gene_name)
    # Build join keys
    left = exploded_df.copy()
    left["_gene_norm"] = left["gene_norm"]  # already normalized
    # Join
    logger.info("Merging exploded cluster genes with reference by normalized gene name.")
    merged = left.merge(ref, left_on="_gene_norm", right_on="_ref_gene_norm", how="left", suffixes=("", "_ref"))
    # Rename reference columns to 'ref_<col>'
    ref_renames = {}
    for c in ref_cols:
        ref_renames[c] = f"ref_{c}"
    merged = merged.rename(columns=ref_renames)
    logger.info(f"Merged exploded records: {merged.shape[0]} rows.")
    return merged


# 6) Aggregate cluster-level summary -----------------------------------------
def aggregate_cluster_summary(merged_exploded: pd.DataFrame, cluster_id_col: str, logger) -> pd.DataFrame:
    """
    6.1 Aggregate merged exploded table to cluster-level summary:
       - matched_genes_list: semicolon-separated matched gene names (raw)
       - matched_count: number matched
       - unmatched_genes_list, unmatched_count
       - matched_ref_info: semicolon-separated 'GENE|ref_col1:val|ref_col2:val' (if ref exists)
    """
    logger.info("Aggregating cluster-level summary.")
    def agg_for_group(gdf):
        raw_genes = gdf["gene_raw"].astype(str).tolist()
        matched_mask = gdf["ref_Associated Gene Name"].notna() if "ref_Associated Gene Name" in gdf.columns else gdf.filter(like="ref_").any(axis=1)
        matched_genes = gdf.loc[matched_mask, "gene_raw"].astype(str).tolist()
        unmatched_genes = [x for x in raw_genes if x not in matched_genes]
        matched_count = len(matched_genes)
        unmatched_count = len(unmatched_genes)
        # Build a compact matched_ref_info if reference columns exist
        ref_cols = [c for c in gdf.columns if c.startswith("ref_")]
        matched_ref_info = []
        if ref_cols:
            for _, row in gdf.loc[matched_mask, :].iterrows():
                gene = row["gene_raw"]
                pieces = [gene]
                for rc in ref_cols:
                    val = row.get(rc)
                    if pd.notna(val) and str(val).strip() != "":
                        pieces.append(f"{rc[4:]}:{str(val)}")
                matched_ref_info.append("|".join(pieces))
        return pd.Series({
            "matched_genes_list": ";".join(matched_genes),
            "matched_count": matched_count,
            "unmatched_genes_list": ";".join(unmatched_genes),
            "unmatched_count": unmatched_count,
            "matched_ref_info": ";".join(matched_ref_info)
        })

    grouped = merged_exploded.groupby(cluster_id_col).apply(agg_for_group).reset_index()
    logger.info(f"Aggregated cluster summary for {len(grouped)} clusters.")
    return grouped


# 7) Save outputs (CSV + XLSX workbook) --------------------------------------
def save_results(exploded_merged: pd.DataFrame, cluster_summary: pd.DataFrame,
                 outdir: Path, base_name: str, ts: str, logger):
    """7.1 Save exploded merged CSV, cluster summary CSV, and XLSX workbook with sheets."""
    ensure_dir(outdir)
    exploded_csv = outdir / f"{base_name}_exploded_merged_{ts}.csv"
    summary_csv = outdir / f"{base_name}_cluster_summary_{ts}.csv"
    workbook = outdir / f"{base_name}_merged_{ts}.xlsx"

    logger.info(f"Saving exploded merged CSV: {exploded_csv.name}")
    exploded_merged.to_csv(exploded_csv, index=False, encoding="utf-8")

    logger.info(f"Saving cluster summary CSV: {summary_csv.name}")
    cluster_summary.to_csv(summary_csv, index=False, encoding="utf-8")

    # Save XLSX workbook with two sheets
    try:
        logger.info(f"Saving combined Excel workbook: {workbook.name}")
        with pd.ExcelWriter(workbook, engine="openpyxl") as writer:
            cluster_summary.to_excel(writer, sheet_name="cluster_summary", index=False)
            exploded_merged.to_excel(writer, sheet_name="exploded_merged", index=False)
        logger.info("Workbook saved.")
    except Exception as e:
        logger.warning(f"Could not save workbook ({e}); CSV outputs are available.")

    return exploded_csv, summary_csv, workbook


# 8) CLI ---------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Merge cluster gene lists with CiliaCarta reference by gene name.")
    p.add_argument("--clusters", "-c", required=True, help="Path to clusters CSV/XLSX (filtered for cilia).")
    p.add_argument("--clusters-sheet", default=None, help="If clusters input is Excel, sheet name/index.")
    p.add_argument("--ciliacarta", "-r", required=True, help="Path to CiliaCarta CSV/XLSX reference table.")
    p.add_argument("--ciliacarta-sheet", default=None, help="If CiliaCarta input is Excel, sheet name/index.")
    p.add_argument("--outdir", "-o", default="output_data", help="Output directory.")
    p.add_argument("--cluster-id-col", default="Cluster_ID", help="Column name in clusters table that identifies cluster.")
    p.add_argument("--cluster-genes-col", default="cluster_genes", help="Column name with genes list in clusters table.")
    p.add_argument("--ciliacarta-gene-col", default="Associated Gene Name", help="Gene name column in CiliaCarta.")
    p.add_argument("--gene-sep", default=r"[;,|\s]+", help="Regex of separators for splitting gene lists (default: '[;,|\\s]+').")
    p.add_argument("--min-matches", type=int, default=0, help="Optionally require clusters to have at least this many matched genes (filter summary).")
    return p.parse_args()


# 9) Main flow ---------------------------------------------------------------
def main():
    args = parse_args()
    clusters_path = Path(args.clusters)
    ref_path = Path(args.ciliacarta)
    outdir = Path(args.outdir)
    ensure_dir(outdir)
    ts = timestamp_now()
    log_path = outdir / f"merge_clusters_{ts}.log"
    logger = setup_logger(log_path)

    logger.info("=== START: merge_clusters ===")
    logger.info(f"Clusters input: {clusters_path}")
    logger.info(f"CiliaCarta ref: {ref_path}")
    logger.info(f"Outdir: {outdir}")

    try:
        # 1. Read inputs
        logger.info("=== 1. READ INPUTS ===")
        df_clusters = read_table(clusters_path, args.clusters_sheet, logger)
        df_ref = read_table(ref_path, args.ciliacarta_sheet, logger)

        # 2. Validate columns
        logger.info("=== 2. VALIDATE COLUMNS ===")
        if args.cluster_id_col not in df_clusters.columns:
            logger.error(f"Cluster ID column '{args.cluster_id_col}' not found in clusters input.")
            raise KeyError(args.cluster_id_col)
        if args.cluster_genes_col not in df_clusters.columns:
            logger.error(f"Cluster genes column '{args.cluster_genes_col}' not found in clusters input.")
            raise KeyError(args.cluster_genes_col)
        if args.ciliacarta_gene_col not in df_ref.columns:
            logger.error(f"CiliaCarta gene column '{args.ciliacarta_gene_col}' not found in reference input.")
            raise KeyError(args.ciliacarta_gene_col)

        # 3. Normalize reference gene names
        logger.info("=== 3. NORMALIZE REFERENCE ===")
        df_ref_clean = df_ref.copy()
        df_ref_clean["_ref_gene_norm"] = df_ref_clean[args.ciliacarta_gene_col].fillna("").astype(str).apply(normalize_gene_name)

        # 4. Explode cluster genes
        logger.info("=== 4. EXPLODE CLUSTER GENES ===")
        exploded = explode_cluster_genes(df_clusters, genes_col=args.cluster_genes_col,
                                         cluster_id_col=args.cluster_id_col,
                                         sep_pattern=args.gene_sep, logger=logger)
        if exploded.empty:
            logger.warning("No exploded rows produced; exiting.")
            return 0

        # 5. Merge exploded with reference
        logger.info("=== 5. MERGE WITH REFERENCE ===")
        merged = merge_with_reference(exploded, df_ref_clean, ref_gene_col=args.ciliacarta_gene_col, logger=logger)

        # 6. Aggregate to cluster summary
        logger.info("=== 6. AGGREGATE TO CLUSTER LEVEL ===")
        cluster_summary = aggregate_cluster_summary(merged, cluster_id_col=args.cluster_id_col, logger=logger)

        # 7. Optionally filter clusters with few matches
        if args.min_matches > 0:
            pre_count = cluster_summary.shape[0]
            cluster_summary = cluster_summary[cluster_summary["matched_count"] >= args.min_matches].reset_index(drop=True)
            post_count = cluster_summary.shape[0]
            logger.info(f"Filtered clusters by min_matches={args.min_matches}: {pre_count} -> {post_count}")

        # 8. Save outputs
        logger.info("=== 7. SAVE OUTPUTS ===")
        base_name = clusters_path.stem
        exploded_csv, summary_csv, workbook = save_results(merged, cluster_summary, outdir, base_name, ts, logger)

        # 9. Provenance & summary
        logger.info("=== 8. PROVENANCE & SUMMARY ===")
        logger.info(f"Exploded merged CSV: {exploded_csv.name}")
        logger.info(f"Cluster summary CSV: {summary_csv.name}")
        logger.info(f"Workbook: {workbook.name}")
        logger.info("=== END: merge_clusters ===")
        return 0

    except Exception as e:
        logger.exception(f"Fatal error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
