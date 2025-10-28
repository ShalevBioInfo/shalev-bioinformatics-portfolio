# README — Heatmap Generation Scripts (LPP & NPP)

**Purpose:** A collection of R scripts to produce publication-ready phylogenetic profile heatmaps (LPP and NPP) for user-provided gene lists, demonstrating various `ComplexHeatmap` functionalities.

**Scripts Included:**

1.  `gene_list_to_lpp_heatmap.R`: Generates a standard LPP heatmap for a single gene list with hierarchical clustering.
2.  `gene_list_to_npp_heatmap.R`: Generates a standard z-scored NPP heatmap for a single gene list, with options for row order and color clipping.
3.  `lpp_multi_cluster_heatmap_with_inclusion.R`: Generates an advanced LPP heatmap displaying multiple clusters (from a pre-annotated input file) with row gaps, per-cluster clustering, and evidence annotations (`Inclusion_criterion`).

Each script generally:

* Reads gene lists and profile matrices (LPP or NPP).
* Aligns species columns to a clades mapping file (maintaining fixed phylogenetic order).
* Determines row order (clustering, input order, or pre-defined).
* Renders a `ComplexHeatmap` figure (PNG, optionally PDF) and writes sidecar files for reproducibility (e.g., sorted matrix, row order, session info, run parameters).

---

## Files & Locations

Recommended repository structure for this project:

```

projects/heatmap\_visualization/
├─ scripts/
│  ├─ gene\_list\_to\_lpp\_heatmap.R
│  ├─ gene\_list\_to\_npp\_heatmap.R
│  └─ lpp\_multi\_cluster\_heatmap\_with\_inclusion.R \# Added script
├─ input/                 \# Small demo inputs (synthetic)
│  ├─ demo\_genes.csv
│  ├─ demo\_clusters\_genes\_inclusion.csv \# Input for the 3rd script
│  ├─ demo\_lpp.tsv
│  ├─ demo\_npp.tsv
│  └─ species\_clades.csv
├─ results/               \# Example outputs
│  └─ example\_multi\_cluster\_heatmap.png \# Added example output
└─ README.md              \# This file

````

---

## Requirements

* R 4.0+ recommended.
* CRAN packages: `optparse`, `readr`, `stringr`, `dplyr`, `glue`, `data.table`, `tidyr`, `tools`, `randomcoloR`.
* Bioconductor packages: `ComplexHeatmap`, `circlize`.
* **Note:** All scripts include automatic dependency checking and installation via `BiocManager`.
* Recommended reproducible environment: `renv` or Docker image with R + Bioconductor.
* On some Windows setups, you may need Cairo support for high-quality PNGs (`type="cairo"`).

---

## Input Formats (Examples)

### 1) Gene list (for scripts 1 & 2)

Either `--genes "GENE1,GENE2,..."` or a CSV (`--csv path`) with gene symbols.

Example `demo_genes.csv`:

```csv
gene
ABCA4
RHO
USH2A
````

### 2\) Cluster Details File (for script 3)

CSV with **one gene per row**, including cluster and evidence info. Required columns: `cluster_id`, `cluster_genes`, `Inclusion_criterion`.

Example `demo_clusters_genes_inclusion.csv`:

```csv
cluster_id,cluster_genes,Inclusion_criterion
949,MT-ND5,Literature
949,MT-CYB,Literature
1257,TTC30B,CiliaCarta: Gold Standard
1257,CLUAP1,CiliaCarta: Predicted
```

### 3\) Profile Matrix (LPP / NPP)

Tab/CSV/RDS with first column = gene names, remaining columns = species identifiers (taxid or scientific\_name).

Example (TSV):

```
gene    9606    10090   3702
ABCA4   0.95    0.85    0.02
RHO     0.99    0.88    0.01
```

### 4\) Clades Mapping

CSV with columns `scientific_name`, `taxid`, `clade`. Species order defines fixed column order.

Example:

```csv
scientific_name,taxid,clade
Homo_sapiens,9606,Mammalia
Mus_musculus,10090,Mammalia
Danio_rerio,7955,Actinopteri
```

-----

## How to Run — Examples

*(Note: Assumes scripts are in `scripts/` and demo data in `input/` relative to where you run the command)*

### Standard LPP Heatmap (`gene_list_to_lpp_heatmap.R`)

```bash
Rscript scripts/gene_list_to_lpp_heatmap.R \
  --genes "ABCA4,RHO,USH2A" \
  --lpp input/demo_lpp.tsv \
  --clades input/species_clades.csv \
  --outdir heatmap_runs \
  --out_prefix demo_lpp \
  --min_height_in 8
```

### Standard NPP Heatmap (`gene_list_to_npp_heatmap.R`)

```bash
Rscript scripts/gene_list_to_npp_heatmap.R \
  --csv input/demo_genes.csv \
  --npp input/demo_npp.tsv \
  --clades input/species_clades.csv \
  --row_order cluster \
  --z_clip 3 \
  --outdir heatmap_runs \
  --out_prefix demo_npp
```

### Multi-Cluster Annotated LPP Heatmap (`lpp_multi_cluster_heatmap_with_inclusion.R`)

```bash
Rscript scripts/lpp_multi_cluster_heatmap_with_inclusion.R \
  --clusters input/demo_clusters_genes_inclusion.csv \
  --lpp input/demo_lpp.tsv \
  --species input/species_clades.csv \
  --selected "949,1257" \
  --outdir heatmap_runs \
  --label demo_multi_cluster \
  --seed 42
```

-----

## Main CLI Options (Summary)

*(Run script with `--help` for full details)*

**Common Options (Most Scripts):**

  * `--clusters`, `--csv`, `--genes`: Input gene list / cluster details.
  * `--lpp`, `--npp`: Input profile matrix.
  * `--species`, `--clades`: Input species-clade mapping.
  * `--outdir`: Output directory.
  * `--label`, `--out_prefix`: Prefix for output files.
  * `--seed`: Random seed for colors (where applicable).
  * `--png_width`, `--min_height_in`, `--png_height_per_gene`: Figure dimension controls.

**NPP Specific (`gene_list_to_npp_heatmap.R`):**

  * `--drop_first_rowcol`: Handle meta headers.
  * `--z_clip`: Symmetric color scale clipping.
  * `--row_order`: `'cluster'`, `'input'`, or path to file.

**Multi-Cluster Specific (`lpp_multi_cluster_heatmap_with_inclusion.R`):**

  * `--selected`: Comma-separated cluster IDs to include.
  * `--allow_explode`: Handle semicolon-separated genes in input.
  * `--method`: Clustering method within clusters (e.g., `'average'`, `'single'`).

-----

## Outputs (Per Run)

Outputs are typically written under a timestamped sub-folder within `--outdir`:

  * `<prefix>_<timestamp>.png`: The main heatmap figure.
  * (Optional PDF): Vector version of the figure.
  * `<prefix>_sorted.csv` / `<prefix>_genes_map.csv`: Data underlying the heatmap rows.
  * `<prefix>_row_order.txt`: Final order of genes in the heatmap.
  * `runinfo_<timestamp>.txt`: Parameters and metadata for the run.
  * `runlog_<timestamp>.txt`: Console log output.
  * (Optional `session_info.txt`): R session details.

-----

## Example Output

An example of the multi-cluster annotated heatmap produced by `lpp_multi_cluster_heatmap_with_inclusion.R` can be found in the `results/` folder:

  * `results/example_multi_cluster_heatmap.png`

This demonstrates the visualization of distinct clusters with internal clustering and associated evidence annotations.

projects/heatmap_visualization/results/Clusters_Dominated_by_Known_Genes.png
projects/heatmap_visualization/results/LPP_7_genes_HeatMap_2025.png
projects/heatmap_visualization/results/NPP_7_genes_HeatMap_2025.png

-----

## Key Design Choices & Rationale (Short)

  * **Clustering:** Default uses `1 - Pearson` distance and `average` or `single` linkage for robust clustering of profiles. Handles NAs via `pairwise.complete.obs`.
  * **Columns (Species):** Always kept in fixed phylogenetic order (from species map) for biological interpretability; never clustered.
  * **Color Scales:** LPP uses `white` to `blue` (0 to 1). NPP uses a diverging `blue`-`white`-`red` scale centered at 0, with clipping options.
  * **Reproducibility:** Automatic dependency handling, detailed run info logging, and use of seeds aim for maximum reproducibility.

-----

## Troubleshooting (Quick)

  * **Errors reading matrix/CSV:** Check separators (tab vs comma), headers, and ensure first column contains gene names as expected.
  * **Genes not found / Missing from LPP:** Verify gene symbols match between input list and matrix rownames (case-insensitive rescue is attempted).
  * **Clade mapping issues:** Ensure `taxid` or `scientific_name` in the species map matches LPP column headers. Update map if needed.
  * **Large PNG files / Cropping:** Adjust figure dimension parameters (`--png_width`, `--png_height_per_gene`, `--min_height_in`).

-----

## Licensing, Authorship, Contact

  * **Author:** Shalev Yaacov
  * **License:** MIT (see repository `LICENSE`)
  * **Contact:** Refer to the top-level repository README.

-----

