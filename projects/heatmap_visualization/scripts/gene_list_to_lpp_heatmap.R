#!/usr/bin/env Rscript
################################################################################
# Title:  Gene-list -> LPP Heatmap Generator (git-ready, documented)
# File:   gene_list_to_lpp_heatmap.R
# Author: Shalev Yaacov
# Date:   2025-10-21 (revised)
#
# Purpose (one line):
#   Generate publication-ready heatmaps of Local Phylogenetic Profiles (LPP)
#   for a user-specified gene list; export PNG/PDF + sorted CSV + row order.
#
# High-level NOTES (professional warnings in English):
#   - Bioconductor packages (ComplexHeatmap, circlize) are installed via BiocManager
#     if missing — this script will attempt installation when run interactively but
#     in production it's recommended to manage dependencies outside the script
#     (renv, Docker, or a requirements installer).
#   - The default minimum PNG height is intentionally large (36 inches) to avoid
#     cropping on very tall heatmaps; change --min_height_in to a smaller value
#     when generating smaller output. See CLI below.
#   - Single-gene or empty-gene inputs are handled gracefully (no clustering).
#
# Suggested commit message when adding to git:
#   analysis: add gene_list_to_lpp_heatmap.R (LPP heatmap generator, documented)
################################################################################

# ------------------------------
# OPERATIONAL INPUT GUIDE (copy-paste examples)
# ------------------------------
# 1) Gene-list CSV (simple; header 'gene' or first column used):
#    my_genes.csv
#    ----------------
#    gene
#    ABCA4
#    RHO
#    USH2A
#
# 2) LPP matrix (tab/CSV) — first column MUST be gene names, subsequent columns species
#    example (TSV or CSV):
#    ----------------
#    gene    9606    10090   3702
#    ABCA4   1       0.8     0.0
#    RHO     1       0.75    0.0
#
#    - Column names for species may be taxid or scientific_name. See 'clades' mapping.
#
# 3) Clades mapping CSV (required columns: scientific_name,taxid,clade):
#    species_clades.csv
#    ----------------
#    scientific_name,taxid,clade
#    Homo_sapiens,9606,Mammalia
#    Mus_musculus,10090,Mammalia
#
# 4) CLI usage examples:
#    Rscript gene_list_to_lpp_heatmap.R --genes "ABCA4,RHO,USH2A" --lpp path/to/lpp.tsv --clades path/to/clades.csv
#    Rscript gene_list_to_lpp_heatmap.R --csv my_genes.csv --lpp lpp.tsv --clades species_clades.csv --min_height_in 12 --png_width_in 10
#
# 5) Output (created under --outdir/run_YYYYMMDD_HHMMSS/):
#    - <out_prefix>_N_genes_YYYYMMDD.png
#    - <out_prefix>_N_genes_YYYYMMDD.pdf
#    - <out_prefix>_N_genes_YYYYMMDD_sorted.csv
#    - <out_prefix>_N_genes_YYYYMMDD_row_order.txt
#    - session_info.txt, run_params.txt
#
# ------------------------------
# Block 0: Dependencies & helper logging (NOTE: professional logging functions)
# ------------------------------
suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
})

cran_pkgs <- c("optparse","readr","stringr","dplyr","glue","data.table")
bioc_pkgs <- c("ComplexHeatmap","circlize")

for (p in cran_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
for (p in bioc_pkgs) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE)

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(stringr)
  library(dplyr)
  library(glue)
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

log_info <- function(...) cat(glue("[INFO] {format(Sys.time(), '%Y-%m-%d %H:%M:%S')} | ", ... , "\n"))
log_ok   <- function(...) cat(glue("[OK]   {format(Sys.time(), '%Y-%m-%d %H:%M:%S')} | ", ... , "\n"))
log_warn <- function(...) cat(glue("[WARN] {format(Sys.time(), '%Y-%m-%d %H:%M:%S')} | ", ... , "\n"))
log_err  <- function(...) { cat(glue("[ERR]  {format(Sys.time(), '%Y-%m-%d %H:%M:%S')} | ", ... , "\n")); quit(status = 1) }

# fallback color palette for clade annotation (modify as lab palette evolves)
all_clade_colors <- c(
  Mammalia="#E31A1C", Aves="#A6CEE3", Lepidosauria="#004529",
  Amphibia="#80CDC1", Actinopteri="#1F78B4", Platyhelminthes="#3F007D",
  Nematoda="#6A3D9A", Arthropoda="#CAB2D6", Ascomycota="#FFFF99",
  Basidiomycota="#FF7F00", Viridiplantae="#33A02C", Rhodophyta="#7F0000",
  Ochrophyta="#B2DF8A", Oomycota="#B15928", Apicomplexa="#FDBF6F",
  Amoebozoa="#EF6548", Unknown="#BBBBBB"
)

# ------------------------------
# Block 1: CLI options (USER-FACING PARAMETERS)
# NOTE: min_height_in kept intentionally high by default to match user's previous workflow;
#       adjust via CLI if you want smaller figures (common values: 6, 8, 12 inches).
# ------------------------------
option_list <- list(
  make_option("--genes", type = "character", default = NULL,
              help = "Comma-separated gene symbols, e.g. \"ABCA4,RHO,USH2A\""),
  make_option("--csv",   type = "character", default = "my_genes.csv",
              help = "CSV with gene list (uses 'gene' column or first column)"),
  make_option("--lpp",   type = "character", default = "Homo_sapiens_LPP_matrix.tsv",
              help = "Path to LPP matrix (first column = gene names; subsequent columns = species)"),
  make_option("--clades",type = "character", default = "species_clades.csv",
              help = "CSV with columns: scientific_name,taxid,clade"),
  make_option("--outdir", type = "character", default = "heatmap_runs",
              help = "Output directory base (default: heatmap_runs)"),
  make_option("--out_prefix", type = "character", default = "heatmap",
              help = "Output filename prefix (default: heatmap)"),
  make_option("--min_height_in", type = "double", default = 36,
              help = "Minimum PNG height in inches (default 36; reduce if desired)"),
  make_option("--png_width_in", type = "double", default = 12,
              help = "PNG width in inches (default 12)")
)
args <- parse_args(OptionParser(option_list = option_list))
log_info("CLI options parsed.")
log_ok(glue("Inputs: csv={args$csv} | lpp={args$lpp} | clades={args$clades} | outdir={args$outdir}"))

if (args$min_height_in > 20) {
  log_warn(glue("min_height_in is {args$min_height_in} in — unusually large; keep if intentional."))
}

# ------------------------------
# Block 2: Load/validate gene list (OPERATIONAL: supported input forms)
# ------------------------------
sanitize_genes <- function(x) { x <- trimws(x); x <- x[x != "" & !is.na(x)]; unique(x) }

get_genes_from_cli <- function(s) {
  if (is.null(s)) return(NULL)
  str_split(s, ",", simplify = TRUE) %>% as.vector() %>% sanitize_genes()
}

get_genes_from_csv <- function(path) {
  if (!file.exists(path)) log_err(glue("Gene CSV not found: {path}"))
  log_info(glue("Reading gene CSV: {path}"))
  df <- suppressWarnings(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
  if (nrow(df) == 0) log_err("Gene CSV is empty.")
  genes <- if ("gene" %in% names(df)) df$gene else df[[1]]
  sanitize_genes(genes)
}

genes <- get_genes_from_cli(args$genes)
if (is.null(genes)) {
  log_info("No --genes provided; reading --csv.")
  genes <- get_genes_from_csv(args$csv)
}
if (is.null(genes) || length(genes) == 0) log_err("No genes loaded. Use --genes or provide a valid CSV.")
if (any(duplicated(genes))) {
  log_warn(glue("Found and removed {sum(duplicated(genes))} duplicate gene(s)."))
  genes <- unique(genes)
}
log_ok(glue("Loaded {length(genes)} genes (examples: {paste(head(genes,5), collapse=', ')})"))

# ------------------------------
# Block 3: Load LPP matrix (format expectations & quick checks)
# NOTE: matrix must be numeric; first column = gene names; header present.
# ------------------------------
if (!file.exists(args$lpp)) log_err(glue("LPP file not found: {args$lpp}"))
log_info(glue("Reading LPP matrix: {args$lpp}"))
lpp_dt <- suppressWarnings(fread(args$lpp, header = TRUE, data.table = FALSE))
if (ncol(lpp_dt) < 2) log_err("LPP must have at least two columns (gene + >=1 species).")

rn <- as.character(lpp_dt[[1]])
lpp_dt[[1]] <- NULL
LPP <- as.matrix(lpp_dt)
rownames(LPP) <- rn

if (!is.numeric(LPP[1,1])) log_err("LPP content is not numeric — check file separators and decimals.")
val_range <- range(LPP, na.rm = TRUE)
if (val_range[1] < 0 || val_range[2] > 1) {
  log_warn(glue("LPP numeric range out of [0,1]: {round(val_range[1],4)}..{round(val_range[2],4)}"))
} else {
  log_ok(glue("LPP numeric range looks OK: {round(val_range[1],4)}..{round(val_range[2],4)}"))
}
na_count <- sum(is.na(LPP))
if (na_count > 0) log_warn(glue("LPP contains {na_count} NA cells; correlations will use pairwise.complete.obs."))

log_ok(glue("Loaded LPP with {nrow(LPP)} genes x {ncol(LPP)} species."))

# ------------------------------
# Block 4: Clade mapping & annotation (robust matching)
# NOTE: clades CSV must contain 'scientific_name','taxid','clade'
# ------------------------------
if (!file.exists(args$clades)) log_err(glue("Clades file not found: {args$clades}"))
stc <- suppressWarnings(fread(args$clades, header = TRUE, data.table = FALSE))
req_cols <- c("scientific_name","taxid","clade")
if (!all(req_cols %in% names(stc))) log_err("Clades file must include columns: scientific_name,taxid,clade")
stc$taxid <- as.character(stc$taxid)

n_taxid <- sum(colnames(LPP) %in% stc$taxid)
n_scien <- sum(colnames(LPP) %in% stc$scientific_name)
match_mode <- if (n_taxid >= n_scien) "taxid" else "scientific_name"
log_info(glue("Column matching mode: {match_mode} (taxid_hits={n_taxid}, scien_hits={n_scien})"))

if (match_mode == "scientific_name" && n_scien > 0) {
  mapping <- setNames(stc$taxid, stc$scientific_name)
  new_cols <- as.character(mapping[colnames(LPP)])
  repl_idx <- which(!is.na(new_cols))
  colnames(LPP)[repl_idx] <- new_cols[repl_idx]
}

clade_by_col <- stc$clade[match(colnames(LPP), stc$taxid)]
clade_by_col[is.na(clade_by_col) | clade_by_col == ""] <- "Unknown"
unknown_n <- sum(clade_by_col == "Unknown")
if (unknown_n > 0) log_warn(glue("There are {unknown_n} 'Unknown' clade annotations — consider updating clades CSV."))

log_ok("Clade mapping prepared (vector aligned to LPP columns).")

# ------------------------------
# Block 5: Subset LPP by gene list (case-insensitive rescue & maintain input order)
# ------------------------------
present_exact <- intersect(genes, rownames(LPP))
missing_now <- setdiff(genes, rownames(LPP))

if (length(missing_now)) {
  rn_lower <- tolower(rownames(LPP))
  gn_lower <- tolower(missing_now)
  hit_idx <- match(gn_lower, rn_lower, nomatch = NA_integer_)
  rescued <- !is.na(hit_idx)
  if (any(rescued)) {
    rescued_in <- missing_now[rescued]
    rescued_map <- rownames(LPP)[hit_idx[rescued]]
    log_warn(glue("Rescued {sum(rescued)} gene(s) by case-insensitive match. Example: {paste(head(paste0(rescued_in,'→',rescued_map),5), collapse=', ')}"))
    present_exact <- c(present_exact, rescued_map)
    missing_now <- setdiff(missing_now, rescued_in)
  }
}

present <- unique(present_exact)
missing <- missing_now

if (length(present) == 0) log_err("None of the provided genes were found in the LPP matrix.")
if (length(missing) > 0) log_warn(glue("Missing {length(missing)} gene(s) dropped (example: {paste(head(missing,5), collapse=', ')})"))

LPP_sub <- LPP[present, , drop = FALSE]
ord_input <- match(tolower(rownames(LPP_sub)), tolower(genes))
LPP_sub <- LPP_sub[order(ord_input, na.last = TRUE), , drop = FALSE]
log_ok(glue("Submatrix prepared: {nrow(LPP_sub)} genes x {ncol(LPP_sub)} species."))

# ------------------------------
# Block 6: Clustering (1 - Pearson, average linkage) — handle single-gene case
# ------------------------------
do_clustering <- TRUE
row_hclust <- NULL

if (nrow(LPP_sub) < 2) {
  log_warn("Less than 2 genes: clustering skipped. Output will use input order.")
  do_clustering <- FALSE
} else {
  # compute pairwise Pearson correlation across species (pairwise.complete.obs)
  row_cor <- suppressWarnings(cor(t(LPP_sub), method = "pearson", use = "pairwise.complete.obs"))
  if (any(is.na(row_cor))) log_warn("Correlation matrix contains NA values; verify rows with zero variance.")
  row_dist <- as.dist(1 - row_cor)
  row_hclust <- hclust(row_dist, method = "average")
  log_ok("Row clustering computed (1 - Pearson; average linkage).")
}

# ------------------------------
# Block 7: Build heatmap object & export (geometry, safe defaults, notes)
# NOTE: min_height_in default = 36 in to avoid cropping; adjust via CLI.
# ------------------------------
run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_dir <- file.path(args$outdir, glue("run_{run_stamp}"))
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

n_genes <- nrow(LPP_sub)
date_str <- format(Sys.Date(), "%Y%m%d")
out_base <- file.path(run_dir, glue("{args$out_prefix}_{n_genes}_genes_{date_str}"))

height_cm_per_row <- 0.35  # visual density parameter (tweakable)
heatmap_body_cm <- n_genes * height_cm_per_row
pad_mm <- c(16, 16, 24, 16)
extra_top_anno_mm <- 3
total_vert_in <- (pad_mm[1] + pad_mm[3] + extra_top_anno_mm) / 25.4

png_width_in <- args$png_width_in
png_height_in <- max(args$min_height_in, heatmap_body_cm / 2.54 + total_vert_in) # safe default preserved

log_info(glue("Figure dimensions (in): {png_width_in} x {round(png_height_in,1)}"))

# build clade annotation object (stable color mapping)
col_map <- list(Clade = all_clade_colors)
clade_annot <- HeatmapAnnotation(
  which = "column",
  Clade = factor(clade_by_col, levels = names(all_clade_colors)),
  col = col_map,
  annotation_height = unit(2, "mm"),
  simple_anno_size_adjust = TRUE,
  gap = unit(0.3, "mm"),
  show_annotation_name = FALSE,
  show_legend = TRUE
)

bottom_pad <- HeatmapAnnotation(.pad = anno_empty(border = FALSE), annotation_height = unit(18, "mm"))

hm <- Heatmap(
  LPP_sub,
  name = "LPP",
  height = unit(heatmap_body_cm, "cm"),
  row_names_gp = gpar(fontsize = 7),
  col = colorRamp2(c(0, 1), c("white", "#0066cc")),
  cluster_rows = if (do_clustering) as.dendrogram(row_hclust) else FALSE,
  cluster_columns = FALSE,
  top_annotation = clade_annot,
  bottom_annotation = bottom_pad,
  border = TRUE,
  show_heatmap_legend = TRUE,
  show_row_names = TRUE,
  show_column_names = FALSE
)

# Export PNG (high-quality raster) and PDF (vector)
png_file <- glue("{out_base}.png")
pdf_file <- glue("{out_base}.pdf")

# Warning: 'type="cairo"' may require cairo device on some Windows setups.
png(png_file, width = png_width_in, height = png_height_in, units = "in", res = 300, type = "cairo")
draw(hm, padding = unit(pad_mm, "mm"))
dev.off()

pdf(pdf_file, width = png_width_in, height = png_height_in)
draw(hm, padding = unit(pad_mm, "mm"))
dev.off()

log_ok(glue("Saved PNG: {png_file}"))
log_ok(glue("Saved PDF: {pdf_file}"))

# ------------------------------
# Block 8: Export sorted matrix & row order (practical files for downstream use)
# ------------------------------
if (do_clustering) {
  ord <- row_hclust$order
  sorted_genes <- rownames(LPP_sub)[ord]
  LPP_sorted <- LPP_sub[sorted_genes, , drop = FALSE]
} else {
  sorted_genes <- rownames(LPP_sub)
  LPP_sorted <- LPP_sub
}

out_csv <- glue("{out_base}_sorted.csv")
fwrite(as.data.frame(cbind(gene = rownames(LPP_sorted), LPP_sorted)), out_csv)
writeLines(sorted_genes, glue("{out_base}_row_order.txt"))

log_ok(glue("Saved sorted CSV: {out_csv}"))
log_ok(glue("Saved row order: {out_base}_row_order.txt"))

# ------------------------------
# Block 9: Session info & run parameters (reproducibility)
# ------------------------------
capture.output(sessionInfo(), file = file.path(run_dir, "session_info.txt"))

params_list <- list(
  run_stamp = run_stamp,
  args = args,
  n_genes = n_genes,
  lpp_file = args$lpp,
  clades_file = args$clades,
  out_base = out_base
)
writeLines(glue::glue_collapse(paste(names(params_list), " : ", lapply(params_list, toString)), sep = "\n"), file.path(run_dir, "run_params.txt"))

log_ok("Saved session info and run parameters.")
log_ok(glue("Run complete. All outputs are in: {run_dir}"))
