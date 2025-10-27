#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# File: cluster_evidence_viz.py
# Purpose: Create reproducible evidence-summary tables and publication-ready
#          visualizations (heatmap and bar plots) that summarize categorical
#          evidence/annotation counts per cluster (e.g., "Inclusion_criterion")
# Created: 2025-10-27 (refactor)
# Author: Shalev Yaacov (refactored for GitHub publication)
#
# Short description:
#  - Reads a CSV/XLSX master table with at least a cluster identifier column and
#    a categorical evidence/annotation column.
#  - Computes cluster x category presence/count matrix (and optional normalized
#    frequencies), saves matrix as CSV, and produces a heatmap PNG and a barplot
#    PNG (top clusters by total evidence).
#  - Fully parameterized via CLI, logs to console and a timestamped logfile,
#    and writes all outputs to an output directory with timestamped filenames.
#
# Example:
#  python cluster_evidence_viz.py \
#    --input data/filtered_cilia_clusters.csv \
#    --cluster-col Cluster_ID --evidence-col Inclusion_criterion \
#    --outdir output_data --target "Novel cilia-associated candidate" \
#    --top-clusters 40 --method jaccard
#
# -----------------------------------------------------------------------------

# 1) Imports & environment ----------------------------------------------------
import argparse
import logging
import sys
from pathlib import Path
from datetime import datetime
from typing import Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# seaborn optional (nicer defaults); fallback to matplotlib if not available
try:
    import seaborn as sns  # type: ignore
    SEABORN_AVAILABLE = True
except Exception:
    SEABORN_AVAILABLE = False

# 2) Utilities: timestamp, dirs, logging -------------------------------------
def now_ts() -> str:
    """Return compact timestamp string for filenames."""
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def ensure_dir(p: Path):
    """Create directory if it doesn't exist."""
    p.mkdir(parents=True, exist_ok=True)


def setup_logger(log_path: Path):
    """
    Setup logger that prints to stdout and writes to a file. Returns logger.
    """
    logger = logging.getLogger("cluster_evidence_viz")
    logger.setLevel(logging.INFO)
    # clear old handlers to avoid duplicates in interactive runs
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


# 3) I/O helpers --------------------------------------------------------------
def read_table(path: Path, sheet: Optional[str], logger: logging.Logger) -> pd.DataFrame:
    """Read CSV/TSV/XLSX into a pandas DataFrame with trimmed column names."""
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
        # fallback to csv parsing
        df = pd.read_csv(path, dtype=str, encoding="utf-8", low_memory=False)
    # normalize column names
    df.columns = [str(c).strip() for c in df.columns]
    logger.info(f"Loaded table with shape {df.shape}")
    return df


# 4) Data processing ---------------------------------------------------------
def normalize_series(s: pd.Series) -> pd.Series:
    """Normalize string series: fillna, str, strip."""
    return s.fillna("").astype(str).str.strip()


def build_presence_count_matrix(df: pd.DataFrame, cluster_col: str, evidence_col: str,
                                logger: logging.Logger) -> pd.DataFrame:
    """
    Build a cluster x evidence-category count/presence matrix.
    Returns a DataFrame indexed by cluster, columns are categories (counts: int).
    """
    logger.info("Building presence/count matrix.")
    temp = df[[cluster_col, evidence_col]].copy()
    temp[cluster_col] = normalize_series(temp[cluster_col])
    temp[evidence_col] = normalize_series(temp[evidence_col])

    # drop empty rows
    temp = temp[(temp[cluster_col] != "") & (temp[evidence_col] != "")]
    if temp.empty:
        logger.error("No valid rows after cleaning cluster/evidence columns.")
        raise ValueError("No valid rows for building matrix.")

    # create boolean presence per cluster-category and count (if duplicate gene entries exist)
    # Use drop_duplicates to count presence (1) per cluster-category instead of raw duplicates
    pres = (temp.drop_duplicates([cluster_col, evidence_col])
            .assign(pres=1)
            .pivot_table(index=cluster_col, columns=evidence_col, values='pres', fill_value=0))

    # ensure integer dtype
    pres = pres.astype(int)
    logger.info(f"Presence matrix shape: {pres.shape} (clusters x categories).")
    return pres


def compute_similarity_matrix(presence_df: pd.DataFrame, method: str, logger: logging.Logger):
    """
    Compute derived matrices from presence:
     - method == 'raw' -> return presence counts (C)
     - method == 'jaccard' -> return jaccard matrix between categories
     - method == 'cooccurrence' -> return C^T C (category x category counts)
    """
    logger.info(f"Computing similarity matrix using method='{method}'.")
    if method == "raw":
        # return presence_df (clusters x categories)
        return presence_df
    # For category-category matrices, compute counts: C^T * C
    cat_counts = presence_df.T.dot(presence_df).astype(int)
    if method == "cooccurrence":
        return cat_counts
    if method == "jaccard":
        # compute jaccard per pair
        diag = np.diag(cat_counts).astype(float)
        union = diag[:, None] + diag[None, :] - cat_counts.values
        with np.errstate(divide='ignore', invalid='ignore'):
            jacc = cat_counts.values.astype(float) / union
            jacc[np.isnan(jacc)] = 0.0
        jacc_df = pd.DataFrame(jacc, index=cat_counts.index, columns=cat_counts.columns)
        return jacc_df
    raise ValueError(f"Unknown method: {method}")


# 5) Plotting helpers --------------------------------------------------------
def plot_heatmap(matrix: pd.DataFrame, outpath: Path, title: str, logger: logging.Logger,
                 figsize=(10, 8), cmap="vlag", annot=False, vmax: Optional[float] = None,
                 vmin: Optional[float] = None, rotate_xticks: bool = True):
    """Plot and save heatmap (seaborn preferred)."""
    logger.info(f"Plotting heatmap to {outpath.name} (shape {matrix.shape}).")
    plt.close("all")
    fig, ax = plt.subplots(figsize=figsize)
    if SEABORN_AVAILABLE:
        sns.set(style="white")
        sns.heatmap(matrix, ax=ax, cmap=cmap, annot=annot, vmax=vmax, vmin=vmin,
                    cbar_kws={"shrink": 0.8})
    else:
        im = ax.imshow(matrix.values, aspect="auto", origin="lower", cmap=cmap, vmax=vmax, vmin=vmin)
        plt.colorbar(im, ax=ax, shrink=0.8)

    # ticks handling: reduce fontsize / rotate if many labels
    ncols = len(matrix.columns)
    nrows = len(matrix.index)
    if rotate_xticks:
        if ncols > 60:
            ax.set_xticks([])
        else:
            ax.set_xticks(range(ncols))
            ax.set_xticklabels(matrix.columns, rotation=90, fontsize=8 if ncols > 30 else 10)
    else:
        ax.set_xticks(range(ncols))
        ax.set_xticklabels(matrix.columns, rotation=90, fontsize=10)

    if nrows > 200:
        ax.set_yticks([])
    else:
        ax.set_yticks(range(nrows))
        ax.set_yticklabels(matrix.index, rotation=0, fontsize=8 if nrows > 30 else 10)

    ax.set_title(title, fontsize=12)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)
    logger.info("Heatmap saved.")


def plot_top_clusters_bar(presence_df: pd.DataFrame, outpath: Path, top_n: int,
                          title: str, logger: logging.Logger, figsize=(10, 6)):
    """Plot stacked bar of top N clusters by total evidence (counts per category)."""
    logger.info(f"Plotting top-{top_n} clusters barplot to {outpath.name}.")
    totals = presence_df.sum(axis=1).sort_values(ascending=False)
    top_idx = totals.head(top_n).index.tolist()
    if not top_idx:
        logger.warning("No clusters to plot in bar chart (empty selection).")
        return
    df_top = presence_df.loc[top_idx]
    # stacked horizontal bar (categories as colored segments)
    plt.close("all")
    fig, ax = plt.subplots(figsize=figsize)
    df_top.plot(kind="bar", stacked=True, ax=ax, width=0.8)
    ax.set_ylabel("Count (presence)")
    ax.set_xlabel("Cluster")
    ax.set_title(title)
    plt.xticks(rotation=45, ha="right", fontsize=8)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)
    logger.info("Barplot saved.")


# 6) CLI ---------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Cluster evidence visualizer: build CSV/heatmap from a master table.")
    p.add_argument("--input", "-i", required=True, help="Path to input CSV/XLSX master table.")
    p.add_argument("--sheet", default=None, help="Sheet name/index if input is Excel.")
    p.add_argument("--outdir", "-o", default="output_data", help="Output directory.")
    p.add_argument("--cluster-col", default=None, help="Column name for cluster identifier (auto-guess if omitted).")
    p.add_argument("--evidence-col", default=None, help="Column name for evidence/category (auto-guess if omitted).")
    p.add_argument("--method", default="cooccurrence", choices=["raw", "cooccurrence", "jaccard"],
                   help="Which matrix to compute: raw (cluster x category), cooccurrence (category counts), or jaccard.")
    p.add_argument("--target", default=None, help="Optional category name to highlight in heatmap.")
    p.add_argument("--top-clusters", type=int, default=40, help="Number of top clusters to show in barplot.")
    p.add_argument("--figsize", default="12,8", help="Figure size as 'width,height' in inches.")
    p.add_argument("--cmap", default="vlag", help="Colormap for heatmap.")
    p.add_argument("--annot", action="store_true", help="Annotate heatmap cells with values.")
    p.add_argument("--max-cols", type=int, default=120, help="If more columns than this, x-tick labels are suppressed.")
    return p.parse_args()


# 7) Main --------------------------------------------------------------------
def main():
    args = parse_args()
    input_path = Path(args.input)
    outdir = Path(args.outdir)
    ensure_dir(outdir)
    ts = now_ts()
    log_path = outdir / f"cluster_evidence_viz_{ts}.log"
    logger = setup_logger(log_path)

    logger.info("=== START: cluster_evidence_viz ===")
    logger.info(f"Input: {input_path}")
    logger.info(f"Outdir: {outdir}")
    logger.info(f"Plot method: {args.method}")

    try:
        # 1) Read input
        df = read_table(input_path, args.sheet, logger)
        if df.empty:
            logger.error("Input table is empty. Exiting.")
            return 1

        # 2) Guess columns if not provided (simple heuristics)
        cluster_col = args.cluster_col
        evidence_col = args.evidence_col
        if cluster_col is None:
            for cand in ["cluster_id", "cluster", "clusterid", "Cluster_ID", "cluster name"]:
                if cand in df.columns:
                    cluster_col = cand
                    break
        if evidence_col is None:
            for cand in ["Inclusion_criterion", "inclusion_criterion", "Inclusion criterion", "evidence", "category", "Annotation"]:
                if cand in df.columns:
                    evidence_col = cand
                    break
        if cluster_col is None or evidence_col is None:
            logger.error("Could not auto-detect cluster or evidence column. Please pass --cluster-col and --evidence-col explicitly.")
            return 1
        logger.info(f"Using cluster column: '{cluster_col}'")
        logger.info(f"Using evidence column: '{evidence_col}'")

        # 3) Build presence/count matrix
        presence = build_presence_count_matrix(df, cluster_col, evidence_col, logger)
        base = input_path.stem

        # 4) Save presence matrix CSV
        matrix_csv = outdir / f"{base}_cluster_category_matrix_{ts}.csv"
        logger.info(f"Saving presence matrix CSV: {matrix_csv.name}")
        # export with cluster as first column
        pres_export = presence.copy()
        pres_export.insert(0, cluster_col, pres_export.index)
        pres_export.reset_index(drop=True).to_csv(matrix_csv, index=False, encoding="utf-8")

        # 5) Compute requested matrix for plotting
        matrix_for_plot = compute_similarity_matrix(presence, method=args.method, logger=logger)

        # 6) Optionally focus on target (if method produces category x category matrix)
        if args.target and args.method in ("cooccurrence", "jaccard"):
            if args.target not in matrix_for_plot.index:
                logger.warning(f"Target '{args.target}' not found in matrix rows; continuing without highlight.")
            else:
                logger.info(f"Target '{args.target}' found and will be highlighted in the heatmap title.")

        # 7) Limit columns/rows for plotting if very large (take top by variance or counts)
        # convert figsize
        try:
            w, h = [float(x) for x in args.figsize.split(",")]
        except Exception:
            w, h = 12.0, 8.0

        # If plotting category x category and many categories, optionally reduce to top-k by diagonal counts
        plot_matrix = matrix_for_plot
        if args.method in ("cooccurrence", "jaccard"):
            ncats = plot_matrix.shape[0]
            if ncats > args.max_cols:
                # pick top categories by self-count (diagonal) or row sums
                try:
                    diag = np.diag(plot_matrix) if args.method == "cooccurrence" else plot_matrix.sum(axis=1)
                except Exception:
                    diag = plot_matrix.sum(axis=1)
                top_k = args.max_cols
                top_idx = np.argsort(-diag)[:top_k]
                cols_keep = plot_matrix.index[top_idx]
                plot_matrix = plot_matrix.loc[cols_keep, cols_keep]
                logger.info(f"Reducing categories from {ncats} -> {plot_matrix.shape[0]} for plotting (top by support).")

        # 8) Heatmap
        heatmap_png = outdir / f"{base}_evidence_heatmap_{args.method}_{ts}.png"
        title = f"Evidence heatmap ({args.method}) — {base}"
        if args.target:
            title += f" — target: {args.target}"
        plot_heatmap(plot_matrix, heatmap_png, title, logger,
                     figsize=(w, h), cmap=args.cmap, annot=args.annot,
                     vmax=None, vmin=None, rotate_xticks=(plot_matrix.shape[1] <= args.max_cols))

        # 9) Barplot for top clusters
        barplot_png = outdir / f"{base}_top{args.top_clusters}_clusters_bar_{ts}.png"
        plot_top_clusters_bar(presence, barplot_png, top_n=args.top_clusters,
                              title=f"Top {args.top_clusters} clusters by evidence count — {base}",
                              logger=logger, figsize=(min(20, w), min(10, h)))

        # 10) Save a small run summary CSV
        summary_csv = outdir / f"{base}_evidence_summary_{ts}.csv"
        logger.info(f"Saving summary CSV: {summary_csv.name}")
        summary_df = pd.DataFrame({
            "cluster": presence.index,
            "total_evidence_counts": presence.sum(axis=1).astype(int),
            "num_categories": (presence > 0).sum(axis=1).astype(int)
        }).sort_values("total_evidence_counts", ascending=False).reset_index(drop=True)
        summary_df.to_csv(summary_csv, index=False, encoding="utf-8")

        logger.info("=== SUMMARY ===")
        logger.info(f"Saved matrix CSV: {matrix_csv.name}")
        logger.info(f"Saved heatmap PNG: {heatmap_png.name}")
        logger.info(f"Saved barplot PNG: {barplot_png.name}")
        logger.info(f"Saved summary CSV: {summary_csv.name}")
        logger.info("=== END: cluster_evidence_viz ===")
        return 0

    except Exception as e:
        logger.exception(f"Fatal error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
