#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# File: co_occurrence_heatmap_cilia.py
# Purpose: Build a publication-ready co-occurrence matrix and heatmap of
#          categorical annotations (e.g., "Inclusion_criterion") across
#          clusters (e.g., phylogenetic clusters / gene clusters).
# Created: 2025-10-27
# Author: Shalev Yaacov (adapted for GitHub publication)
#
# Short description:
#   - Reads a CSV/XLSX table with at least two columns: cluster identifier
#     and a categorical annotation per row (for example 'Cluster_ID'
#     and 'Inclusion_criterion').
#   - Builds a cluster x category presence matrix (boolean: is category
#     present in cluster?). From that it computes a category x category
#     co-occurrence matrix (counts of clusters where both categories appear).
#   - Supports outputting raw co-occurrence counts or Jaccard similarity.
#   - Produces a labeled heatmap PNG and a CSV sidecar (co-occurrence matrix).
#   - Logs all steps to console and to a timestamped run log in output_data/.
#
# Key design decisions:
#   - Robust validation of input columns (user can override default guesses).
#   - Fast matrix computation uses boolean pivot + matrix multiplication.
#   - Minimal external dependencies: pandas, numpy, matplotlib, seaborn (optional).
#
# Usage example:
#   python co_occurrence_heatmap_cilia.py --input data/example.csv \
#       --outdir output_data --cluster-col Cluster_ID \
#       --category-col Inclusion_criterion --method jaccard \
#       --min-clusters 2 --target "Novel cilia-associated candidate"
#
# -----------------------------------------------------------------------------

# 1. Imports and environment --------------------------------------------------
import argparse
import logging
import sys
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# seaborn is optional but provides nicer heatmaps; fallback to matplotlib if missing
try:
    import seaborn as sns  # type: ignore
    SEABORN_AVAILABLE = True
except Exception:
    SEABORN_AVAILABLE = False

# 2. Helper functions ---------------------------------------------------------
def setup_logger(log_path: Path):
    """2.1 Setup logging to console and to a file (INFO level)."""
    logger = logging.getLogger("coocc")
    logger.setLevel(logging.INFO)
    # Avoid duplicate handlers in interactive runs
    if logger.hasHandlers():
        logger.handlers.clear()

    fmt = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S")

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)
    logger.addHandler(ch)

    fh = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    fh.setLevel(logging.INFO)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    return logger


def timestamp_now() -> str:
    """2.2 Return a compact timestamp string for filenames."""
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def guess_column(df: pd.DataFrame, candidates: list, logger=None):
    """2.3 Find the first column name in df that matches any candidate (case-insensitive)."""
    cols_lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols_lower:
            return cols_lower[cand.lower()]
    # fuzzy: try partial match
    for cand in candidates:
        for col in df.columns:
            if cand.lower() in col.lower():
                if logger:
                    logger.info(f"Guessed column '{col}' for candidate '{cand}'.")
                return col
    return None


def ensure_dir(path: Path):
    """2.4 Create directory if it doesn't exist."""
    path.mkdir(parents=True, exist_ok=True)


# 3. Core processing functions ------------------------------------------------
def read_table(path: Path, logger):
    """3.1 Read CSV/XLSX into pandas DataFrame with safe defaults."""
    logger.info(f"Reading input file: {path}")
    if not path.exists():
        logger.error(f"Input file not found: {path}")
        raise FileNotFoundError(path)
    suffix = path.suffix.lower()
    if suffix in [".csv", ".tsv", ".txt"]:
        sep = "," if suffix == ".csv" else "\t"
        df = pd.read_csv(path, sep=sep, dtype=str, encoding="utf-8", low_memory=False)
    elif suffix in [".xls", ".xlsx"]:
        df = pd.read_excel(path, dtype=str)
    else:
        # try CSV read as default
        df = pd.read_csv(path, dtype=str, encoding="utf-8", low_memory=False)
    df.columns = [c.strip() for c in df.columns]
    logger.info(f"Loaded table with shape {df.shape}")
    return df


def normalize_strings(series: pd.Series) -> pd.Series:
    """3.2 Normalize string values: strip whitespace and convert NaN to empty string."""
    return series.fillna("").astype(str).str.strip()


def build_presence_matrix(df: pd.DataFrame, cluster_col: str, category_col: str, logger):
    """
    3.3 Build presence matrix (clusters x category) boolean:
        - Rows: unique clusters
        - Columns: unique category labels
        - Value = 1 if category appears in cluster, else 0
    """
    logger.info("Building presence matrix (cluster x category).")
    temp = df[[cluster_col, category_col]].copy()
    temp[cluster_col] = normalize_strings(temp[cluster_col])
    temp[category_col] = normalize_strings(temp[category_col])

    # Remove empty categories and clusters
    temp = temp[(temp[cluster_col] != "") & (temp[category_col] != "")]
    if temp.empty:
        logger.error("No valid rows found after cleaning cluster/category columns.")
        raise ValueError("No valid rows in input for cluster/category processing.")

    # Create boolean presence via groupby size and pivot
    pres = (temp
            .drop_duplicates([cluster_col, category_col])
            .assign(pres=1)
            .pivot_table(index=cluster_col, columns=category_col, values='pres', fill_value=0))

    logger.info(f"Presence matrix shape: {pres.shape} (clusters x categories).")
    return pres


def compute_cooccurrence(presence_df: pd.DataFrame, method: str, logger):
    """
    3.4 Compute category x category co-occurrence matrix.
      - method == 'raw' returns counts of clusters where both categories are present.
      - method == 'jaccard' returns intersection / union.
    """
    logger.info(f"Computing co-occurrence matrix using method='{method}'.")
    # Convert to int (0/1)
    mat = presence_df.astype(int)
    # cooccurrence counts: C^T * C -> categories x categories
    counts = mat.T.dot(mat)
    counts = counts.astype(int)

    if method == "raw":
        logger.info("Returning raw co-occurrence counts.")
        return counts

    if method == "jaccard":
        # Intersection = counts, Union = counts_i + counts_j - counts_ij
        diag = np.diag(counts).astype(float)
        union = diag[:, None] + diag[None, :] - counts.values
        with np.errstate(divide='ignore', invalid='ignore'):
            jacc = counts.values.astype(float) / union
            jacc[np.isnan(jacc)] = 0.0
        cooc = pd.DataFrame(jacc, index=counts.index, columns=counts.columns)
        logger.info("Returning Jaccard similarity matrix.")
        return cooc

    logger.error(f"Unknown method '{method}'. Supported: raw, jaccard.")
    raise ValueError("Unknown co-occurrence method.")


def filter_by_min_clusters(counts_df: pd.DataFrame, min_clusters: int, logger):
    """3.5 Drop categories that appear in fewer than min_clusters (based on diagonal counts)."""
    if min_clusters <= 1:
        return counts_df
    logger.info(f"Filtering categories with fewer than {min_clusters} supporting clusters.")
    counts_diag = np.diag(counts_df)
    keep = counts_diag >= min_clusters
    kept_names = counts_df.index[keep].tolist()
    dropped = counts_df.index[~keep].tolist()
    logger.info(f"Keeping {len(kept_names)} categories, dropping {len(dropped)} categories.")
    return counts_df.loc[kept_names, kept_names]


def save_matrix(df: pd.DataFrame, outpath: Path, logger):
    """3.6 Save DataFrame to CSV with timestamp and nice formatting."""
    logger.info(f"Saving matrix CSV to: {outpath}")
    df.to_csv(outpath, index=True, encoding="utf-8")


def plot_heatmap(matrix: pd.DataFrame, outpath: Path, title: str, logger,
                 figsize=(10, 8), cmap="vlag", annot=False, vmax=None, vmin=None, target=None):
    """3.7 Plot and save heatmap (seaborn if available, otherwise matplotlib imshow)."""
    logger.info(f"Plotting heatmap to: {outpath} (size={figsize})")
    plt.close("all")
    fig, ax = plt.subplots(figsize=figsize)

    if SEABORN_AVAILABLE:
        sns.set(style="white")
        sns.heatmap(matrix, ax=ax, cmap=cmap, square=False, cbar_kws={"shrink": 0.8},
                    annot=annot, vmax=vmax, vmin=vmin)
    else:
        im = ax.imshow(matrix.values, aspect="auto", origin="lower", cmap=cmap, vmax=vmax, vmin=vmin)
        plt.colorbar(im, ax=ax, shrink=0.8)

    ax.set_xticks(np.arange(len(matrix.columns)))
    ax.set_xticklabels(matrix.columns, rotation=90, fontsize=8)
    ax.set_yticks(np.arange(len(matrix.index)))
    ax.set_yticklabels(matrix.index, rotation=0, fontsize=8)
    ax.set_title(title, fontsize=12)

    # Optionally highlight target row/column
    if target and target in matrix.index:
        try:
            idx = list(matrix.index).index(target)
            # Draw rectangle around target row/col
            from matplotlib.patches import Rectangle
            ax.add_patch(Rectangle((-0.5, idx - 0.5), len(matrix.columns), 1, fill=False, lw=2, edgecolor="black"))
            ax.add_patch(Rectangle((idx - 0.5, -0.5), 1, len(matrix.index), fill=False, lw=2, edgecolor="black"))
            logger.info(f"Highlighted target category '{target}' in heatmap.")
        except Exception:
            logger.warning("Failed to highlight target; continuing without highlight.")

    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)
    logger.info("Heatmap saved.")


# 4. CLI and main -------------------------------------------------------------
def parse_args():
    """4.1 Parse command-line arguments."""
    p = argparse.ArgumentParser(
        description="Co-occurrence heatmap of categories across clusters (publication-ready).")
    p.add_argument("--input", "-i", required=True, help="Path to input CSV/XLSX table.")
    p.add_argument("--outdir", "-o", default="output_data", help="Output directory.")
    p.add_argument("--cluster-col", default=None,
                   help="Column name for cluster identifier (auto-guessed if omitted).")
    p.add_argument("--category-col", default=None,
                   help="Column name for category / annotation (auto-guessed if omitted).")
    p.add_argument("--method", default="jaccard", choices=["raw", "jaccard"],
                   help="Co-occurrence aggregation: raw counts or jaccard similarity.")
    p.add_argument("--min-clusters", type=int, default=1,
                   help="Minimum number of clusters a category must appear in to be kept.")
    p.add_argument("--target", default=None,
                   help="Optional category name to highlight in heatmap.")
    p.add_argument("--cmap", default="vlag", help="Matplotlib/seaborn colormap.")
    p.add_argument("--figsize", default="10,8",
                   help="Figure size as 'width,height' in inches (e.g. 12,10).")
    return p.parse_args()


def main():
    """4.2 Main execution flow (orchestrates numbered steps and logging)."""
    args = parse_args()

    # Prepare output directory and log
    outdir = Path(args.outdir)
    ensure_dir(outdir)
    ts = timestamp_now()
    log_path = outdir / f"cooccurrence_run_{ts}.log"
    logger = setup_logger(log_path)
    logger.info("=== 1. START: co_occurrence_heatmap_cilia ===")

    input_path = Path(args.input)
    try:
        # Step 1: Read
        logger.info("=== 2. READ INPUT ===")
        df = read_table(input_path, logger)

        # Step 2: Detect columns if not provided
        logger.info("=== 3. DETECT / VALIDATE COLUMNS ===")
        cluster_col = args.cluster_col
        category_col = args.category_col

        if cluster_col is None:
            cluster_col = guess_column(df, ["cluster_id", "cluster", "clusterid", "cluster_id", "Cluster_ID"], logger)
        if category_col is None:
            category_col = guess_column(df, ["Inclusion_criterion", "inclusion_criterion",
                                             "category", "Category", "Annotation"], logger)

        if cluster_col is None or category_col is None:
            logger.error("Could not detect required columns automatically.")
            logger.error("Please pass --cluster-col and --category-col explicitly.")
            raise ValueError("Missing required columns.")

        logger.info(f"Using cluster column: '{cluster_col}'")
        logger.info(f"Using category column: '{category_col}'")

        # Step 3: Build presence matrix
        logger.info("=== 4. BUILD PRESENCE MATRIX ===")
        presence = build_presence_matrix(df, cluster_col, category_col, logger)

        # Step 4: Compute co-occurrence counts / similarity
        logger.info("=== 5. COMPUTE COOCCURRENCE ===")
        coocc = compute_cooccurrence(presence, method=args.method, logger=logger)

        # Step 5: Optionally filter low-support categories
        logger.info("=== 6. FILTER CATEGORIES ===")
        if args.min_clusters > 1:
            # filter based on raw counts diagonal; if method == jaccard, we need raw counts
            raw_counts = presence.T.dot(presence).astype(int)
            filtered_raw = filter_by_min_clusters(raw_counts, args.min_clusters, logger)
            keep_names = filtered_raw.index.tolist()
            coocc = coocc.loc[keep_names, keep_names]
            presence = presence.loc[:, keep_names] if not presence.empty else presence

        # Step 6: Save matrix CSV
        logger.info("=== 7. SAVE MATRIX CSV ===")
        matrix_csv = outdir / f"cooccurrence_matrix_{args.method}_{ts}.csv"
        save_matrix(coocc, matrix_csv, logger)

        # Step 7: Plot heatmap
        logger.info("=== 8. PLOT HEATMAP ===")
        w, h = (float(s) for s in args.figsize.split(","))
        heatmap_png = outdir / f"cooccurrence_heatmap_{args.method}_{ts}.png"
        title = f"Co-occurrence ({args.method}) — {input_path.name}"
        plot_heatmap(coocc, heatmap_png, title, logger,
                     figsize=(w, h), cmap=args.cmap, annot=False, target=args.target)

        # Step 8: Summary log
        logger.info("=== 9. SUMMARY ===")
        logger.info(f"Input: {input_path}")
        logger.info(f"Output dir: {outdir.resolve()}")
        logger.info(f"Saved matrix CSV: {matrix_csv.name}")
        logger.info(f"Saved heatmap PNG: {heatmap_png.name}")
        logger.info("=== 10. END ===")

    except Exception as e:
        logger.exception(f"Fatal error during run: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
