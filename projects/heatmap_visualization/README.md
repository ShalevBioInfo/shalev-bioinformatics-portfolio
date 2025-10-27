# README — Heatmap scripts (LPP & NPP)

**Purpose (short):** Two R scripts to produce publication-ready phylogenetic profile heatmaps for a user-provided gene list.

* `gene_list_to_lpp_heatmap.R` — Local Phylogenetic Profiles (LPP) heatmaps.
* `gene_list_to_npp_heatmap.R` — z-scored Normalized Phylogenetic Profiles (NPP) heatmaps.

Each script:

* reads a gene list and a profile matrix (LPP or NPP),
* aligns species columns to a clades mapping (fixed phylogenetic order),
* determines row order (clustering / input order / explicit file),
* renders a ComplexHeatmap (PNG + PDF) and writes sidecar files (sorted matrix, row order, session info, run params).

---

## Files & locations

Place these files under the repository (recommended paths):

```
projects/heatmap_visualization/
├─ scripts/
│  ├─ gene_list_to_lpp_heatmap.R
│  └─ gene_list_to_npp_heatmap.R
├─ input/                # small demo inputs (synthetic) — recommended
│  ├─ demo_genes.csv
│  ├─ demo_lpp.tsv
│  ├─ demo_npp.tsv
│  └─ species_clades.csv
└─ README.md             # this file
```

---

## Requirements

* R 4.0+ recommended.
* CRAN packages: `optparse`, `readr`, `stringr`, `dplyr`, `glue`, `data.table`.
* Bioconductor packages: `ComplexHeatmap`, `circlize`.
* Recommended reproducible environment: `renv` or Docker image with R + Bioconductor.
* On some Windows setups you may need Cairo support for high-quality PNGs; see *Notes*.

---

## Input formats (examples)

### 1) Gene list

Either provide `--genes "GENE1,GENE2,..."` or a CSV (`--csv path`) where the first column or a `gene` column lists gene symbols.

Example `demo_genes.csv`:

```csv
gene
ABCA4
RHO
USH2A
```

### 2) Profile matrix (LPP / NPP)

Tab/CSV with first column = gene names, remaining columns = species identifiers (taxid or scientific_name).

Example (TSV):

```
gene    9606    10090   3702
ABCA4   0.95    0.85    0.02
RHO     0.99    0.88    0.01
```

### 3) Clades mapping

CSV with columns `scientific_name,taxid,clade`. Species order in this file defines the fixed column order used by the scripts.

Example:

```csv
scientific_name,taxid,clade
Homo_sapiens,9606,Mammalia
Mus_musculus,10090,Mammalia
Danio_rerio,7955,Actinopteri
```

---

## How to run — examples

### LPP heatmap

```bash
# genes provided inline
Rscript scripts/gene_list_to_lpp_heatmap.R \
  --genes "ABCA4,RHO,USH2A" \
  --lpp input/demo_lpp.tsv \
  --clades input/species_clades.csv \
  --outdir heatmap_runs \
  --out_prefix lpp_heatmap \
  --min_height_in 12

# or genes from CSV
Rscript scripts/gene_list_to_lpp_heatmap.R \
  --csv input/demo_genes.csv \
  --lpp input/demo_lpp.tsv \
  --clades input/species_clades.csv
```

### NPP heatmap

```bash
# cluster rows (default), symmetric clipping to ±4
Rscript scripts/gene_list_to_npp_heatmap.R \
  --genes "PNPT1,TUB,USHA1" \
  --npp input/demo_npp.tsv \
  --clades input/species_clades.csv \
  --row_order cluster \
  --z_clip 4 \
  --outdir heatmap_runs \
  --out_prefix npp_heatmap

# keep input row order
Rscript scripts/gene_list_to_npp_heatmap.R \
  --csv input/demo_genes.csv \
  --npp input/demo_npp.tsv \
  --clades input/species_clades.csv \
  --row_order input
```

---

## Main CLI options (summary)

Common options

* `--genes` : comma-separated gene symbols (e.g., `"ABCA4,RHO"`).
* `--csv` : path to gene CSV (first column or `gene` header).
* `--clades` : clades mapping CSV (required).
* `--outdir` : base output dir (default `heatmap_runs`).
* `--out_prefix` : file name prefix for outputs.
* `--min_height_in` : minimum PNG height in inches (defaults kept intentionally high to avoid cropping).

LPP-specific

* `--lpp` : path to LPP matrix file.

NPP-specific

* `--npp` : path to NPP matrix file.
* `--drop_first_rowcol` : `TRUE/FALSE` — drop first row & column if meta-headers exist (default `TRUE`).
* `--z_clip` : numeric — symmetric clip for color scale (e.g., `3` means [-3,3]); if not provided, 1%/99% quantiles are used.
* `--row_order` : `'cluster'` (default), `'input'`, or path to file listing desired row order.

---

## Outputs (per run)

Outputs are written under an auto-created run folder (e.g., `heatmap_runs/run_YYYYMMDD_HHMMSS/`):

* `<out_prefix>.png` — raster figure (high quality).
* `<out_prefix>.pdf` — vector figure.
* `<out_prefix>_sorted.csv` — matrix sorted by final row order (gene first column).
* `<out_prefix>_row_order.txt` — final row order (one gene per line).
* `session_info.txt` — R session info for reproducibility.
* `run_params.txt` — copy of run parameters and key metadata.

---

## Key design choices & rationale (short)

* **Clustering:** rows clustered with `distance = 1 - Pearson` and `average` linkage — robust for continuous profiles; correlation computed with `pairwise.complete.obs` to handle missing values.
* **Columns:** species columns are *not* clustered — they are shown in fixed phylogenetic order from the clades mapping to preserve biological meaning.
* **Color scale:** NPP uses a symmetric diverging scale centered at 0 (z-scores). Use `--z_clip` for a fixed dynamic range or defaults to quantile-based clipping to avoid extreme outliers dominating the palette.
* **Large-minimum height:** default `min_height_in` is high to avoid row-label/figure cropping on tall heatmaps; configurable for compact figures.

---

## Reproducibility & best practices

* Prefer using `renv` or a Docker image to lock package versions. The scripts may install missing Bioconductor packages via `BiocManager` but this is not a replacement for environment management.
* Each run writes `session_info.txt` and `run_params.txt` — keep these with your figures to enable exact reproduction of results.
* For Windows users: Cairo raster support may be required for `type="cairo"` PNG; if unavailable, set `type` parameter or run on Linux.

---

## Troubleshooting (quick)

* **Script errors reading matrix:** check separators (tab vs comma), ensure first column is gene names and header exists.
* **No genes found:** check capitalization and gene symbol conventions; scripts perform case-insensitive rescue but ensure consistent identifiers.
* **Clade mapping misses:** update `species_clades.csv` so that taxid or scientific_name matches column names in the matrix.
* **Huge PNG files:** reduce `--min_height_in` (e.g., 6 or 12) to generate smaller images for screen preview.

---

## Licensing, authorship, contact

* **Author:** Shalev Yaacov
* **License:** MIT (see repository `LICENSE`)
* **Contact:** add preferred email or GitHub handle in the top-level README or `METADATA.md`.

---

## Suggested git commit message (when adding scripts)

```
analysis: add heatmap scripts (gene_list_to_lpp_heatmap.R, gene_list_to_npp_heatmap.R) + README
```

---
