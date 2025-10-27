#!/usr/bin/env Rscript
################################################################################
# Title:  Gene-list -> NPP Heatmap Generator (git-ready, documented)
# File:   gene_list_to_npp_heatmap.R
# Author: Shalev Yaacov
# Date:   2025-10-21 (revised)
#
# Purpose (one line):
#   Generate publication-ready heatmaps of z-scored Normalized Phylogenetic Profiles
#   (NPP) for a user-specified gene list; export PNG/PDF + sorted CSV + row order.
#
# Professional notes:
#   - Bioconductor packages (ComplexHeatmap, circlize) are installed via BiocManager
#     if missing. In production prefer renv/Docker for reproducibility.
#   - Default min PNG height is intentionally large (36 in) to avoid cropping on tall
#     heatmaps; configurable via --min_height_in.
#   - Script writes session_info and run_params for reproducibility.
#
#   analysis: add gene_list_to_npp_heatmap.R (NPP heatmap generator, documented)
################################################################################

# ------------------------------
# OPERATIONAL INPUT GUIDE (copy-paste examples)
# ------------------------------
# 1) Gene-list CSV (header 'gene' or single-column):
#    my_genes.csv
#    ----------------
#    gene
#    ABCA4
#    RHO
#
# 2) NPP matrix (tab/CSV) — first column MUST be gene names, subsequent columns species (taxid or sci_name)
#    example:
#    gene    9606    10090   3702
#    ABCA4   2.1     1.9     -0.1
#    RHO     1.8     1.6     -0.2
#
# 3) Clades mapping CSV (required: scientific_name,taxid,clade)
#    scientific_name,taxid,clade
#    Homo_sapiens,9606,Mammalia
#    Mus_musculus,10090,Mammalia
#
# 4) CLI usage examples:
#    Rscript gene_list_to_npp_heatmap.R --genes "ABCA4,RHO,USH2A" --npp path/to/npp.tsv --clades path/to/clades.csv
#    Rscript gene_list_to_npp_heatmap.R --csv my_genes.csv --npp npp.tsv --clades species_clades.csv --row_order input --z_clip 4
#
# 5) Outputs (under --outdir/<run_folder>/):
#    - <out_prefix>.png, .pdf, _sorted.csv, _row_order.txt, session_info.txt, run_params.txt
#
# ------------------------------
# Block 0: Dependencies & helper logging (professional)
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

# fallback clade color palette (lab can override)
all_clade_colors <- c(
  Mammalia="#E31A1C", Aves="#A6CEE3", Lepidosauria="#004529",
  Amphibia="#80CDC1", Actinopteri="#1F78B4", Platyhelminthes="#3F007D",
  Nematoda="#6A3D9A", Arthropoda="#CAB2D6", Ascomycota="#FFFF99",
  Basidiomycota="#FF7F00", Viridiplantae="#33A02C", Rhodophyta="#7F0000",
  Ochrophyta="#B2DF8A", Oomycota="#B15928", Apicomplexa="#FDBF6F",
  Amoebozoa="#EF6548", Unknown="#BBBBBB"
)

set.seed(42)

# ------------------------------
# Block 1: CLI options (user-facing)
# ------------------------------
option_list <- list(
  make_option("--genes", type = "character", default = NULL,
              help = "Comma-separated gene symbols, e.g. 'ABCA4,RHO,USH2A'"),
  make_option("--csv",   type = "character", default = "my_genes.csv",
              help = "CSV with gene list (header 'gene' or first column)"),
  make_option("--npp",   type = "character", default = "Homo_sapiens_NPP.tsv",
              help = "Path to NPP matrix (first column = gene names)"),
  make_option("--clades",type = "character", default = "species_clades.csv",
              help = "CSV with columns: scientific_name,taxid,clade"),
  make_option("--drop_first_rowcol", type = "logical", default = TRUE,
              help = "If TRUE, drop first row & column (common meta headers)"),
  make_option("--z_clip", type = "double", default = NA,
              help = "Symmetric clip for color scale (e.g., 3 -> [-3,3]); if NA use 1%/99% quantiles"),
  make_option("--row_order", type = "character", default = "cluster",
              help = "Row order: 'cluster' (default), 'input', or path to file with desired order"),
  make_option("--outdir", type = "character", default = "heatmap_runs",
              help = "Output directory base"),
  make_option("--out_prefix", type = "character", default = "npp_heatmap",
              help = "Output filename prefix"),
  make_option("--min_height_in", type = "double", default = 36,
              help = "Minimum PNG height in inches (default 36; change if desired)"),
  make_option("--png_width_in", type = "double", default = 12,
              help = "PNG width in inches (default 12)")
)
args <- parse_args(OptionParser(option_list = option_list))
log_info("Parsed CLI options.")
log_ok(glue("npp={args$npp} | clades={args$clades} | drop_first_rowcol={args$drop_first_rowcol} | row_order={args$row_order}"))

if (args$min_height_in > 20) {
  log_warn(glue("min_height_in is {args$min_height_in} inches — unusually large; keep if intentional."))
}

# ------------------------------
# Block 2: load gene list (support --genes or --csv)
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
  log_warn(glue("Removed {sum(duplicated(genes))} duplicate gene(s)."))
  genes <- unique(genes)
}
log_ok(glue("Loaded {length(genes)} genes (examples: {paste(head(genes,6), collapse=', ')})"))

# ------------------------------
# Block 3: load NPP matrix (format expectations & quick checks)
# ------------------------------
if (!file.exists(args$npp)) log_err(glue("NPP file not found: {args$npp}"))
log_info(glue("Reading NPP matrix: {args$npp}"))
npp_dt <- suppressWarnings(fread(args$npp, header = TRUE, data.table = FALSE))
if (ncol(npp_dt) < 2) log_err("NPP must have >=2 columns (gene + >=1 species).")

rn <- as.character(npp_dt[[1]])
npp_dt[[1]] <- NULL
NPP <- as.matrix(npp_dt)
rownames(NPP) <- rn

# optional drop of first row&col (common meta header); safe-guarded
if (isTRUE(args$drop_first_rowcol)) {
  if (nrow(NPP) >= 2 && ncol(NPP) >= 2) {
    NPP <- NPP[-1, -1, drop = FALSE]
    log_ok("Dropped first row & first column from NPP (per --drop_first_rowcol).")
  } else {
    log_warn("Requested drop_first_rowcol but matrix too small; skipping.")
  }
}

if (!is.numeric(NPP[1,1])) log_err("NPP content is not numeric — check separators/format.")
na_count <- sum(is.na(NPP))
if (na_count > 0) log_warn(glue("NPP contains {na_count} NA cells; correlations will use pairwise.complete.obs."))
rng <- range(NPP, na.rm = TRUE)
log_ok(glue("NPP loaded: {nrow(NPP)} genes x {ncol(NPP)} species. Range: {round(rng[1],3)}..{round(rng[2],3)}"))

# ------------------------------
# Block 4: clade mapping & fixed phylogenetic order
# ------------------------------
if (!file.exists(args$clades)) log_err(glue("Clades file not found: {args$clades}"))
stc <- suppressWarnings(fread(args$clades, header = TRUE, data.table = FALSE))
req_cols <- c("scientific_name","taxid","clade")
if (!all(req_cols %in% names(stc))) log_err("Clades file must include: scientific_name,taxid,clade")
stc$taxid <- as.character(stc$taxid)

# choose matching key
n_taxid <- sum(colnames(NPP) %in% stc$taxid)
n_scien <- sum(colnames(NPP) %in% stc$scientific_name)
match_mode <- if (n_taxid >= n_scien) "taxid" else "scientific_name"
log_info(glue("Column-ID match mode: {match_mode} (taxid_hits={n_taxid}, scien_hits={n_scien})"))

if (match_mode == "scientific_name" && n_scien > 0) {
  mapping <- setNames(stc$taxid, stc$scientific_name)
  new_cols <- as.character(mapping[colnames(NPP)])
  repl_idx <- which(!is.na(new_cols))
  colnames(NPP)[repl_idx] <- new_cols[repl_idx]
}

# reorder columns according to clades file (fixed phylogenetic order)
ordered_taxa <- as.character(stc$taxid)
keep_taxa <- intersect(ordered_taxa, colnames(NPP))
if (length(keep_taxa) == 0) log_warn("No column names matched clades taxid list; proceeding with original column order.")
NPP <- NPP[, keep_taxa, drop = FALSE]

# clade annotation vector aligned to columns
clade_by_col <- stc$clade[match(colnames(NPP), stc$taxid)]
clade_by_col[is.na(clade_by_col) | clade_by_col == ""] <- "Unknown"
clades_factor <- factor(clade_by_col, levels = names(all_clade_colors))
clades_annotation <- HeatmapAnnotation(
  which = "column",
  Clade = clades_factor,
  col = list(Clade = all_clade_colors),
  annotation_height = unit(2, "mm"),
  simple_anno_size_adjust = TRUE,
  gap = unit(0.3, "mm"),
  show_annotation_name = FALSE,
  show_legend = TRUE
)
unknown_n <- sum(clade_by_col == "Unknown")
if (unknown_n > 0) log_warn(glue("There are {unknown_n} 'Unknown' clade annotations; consider updating the clades CSV."))

log_ok(glue("Columns aligned to clade mapping — matched species: {length(keep_taxa)}"))

# ------------------------------
# Block 5: subset by gene list (case-insensitive rescue & maintain order)
# ------------------------------
present_exact <- intersect(genes, rownames(NPP))
missing_now <- setdiff(genes, rownames(NPP))

if (length(missing_now)) {
  rn_lower <- tolower(rownames(NPP))
  gn_lower <- tolower(missing_now)
  hit_idx <- match(gn_lower, rn_lower, nomatch = NA_integer_)
  rescued <- !is.na(hit_idx)
  if (any(rescued)) {
    rescued_in <- missing_now[rescued]
    rescued_map <- rownames(NPP)[hit_idx[rescued]]
    log_warn(glue("Rescued {sum(rescued)} gene(s) by case-insensitive match (examples: {paste(head(paste0(rescued_in,'→',rescued_map),6), collapse=', ')} )"))
    present_exact <- c(present_exact, rescued_map)
    missing_now <- setdiff(missing_now, rescued_in)
  }
}

present <- unique(present_exact)
missing <- missing_now

if (length(present) == 0) log_err("None of the provided genes were found in NPP rownames.")
if (length(missing) > 0) log_warn(glue("Missing {length(missing)} gene(s) dropped (examples: {paste(head(missing,6), collapse=', ')})"))

NPP_sub <- NPP[present, , drop = FALSE]
ord_input <- match(tolower(rownames(NPP_sub)), tolower(genes))
NPP_sub <- NPP_sub[order(ord_input, na.last = TRUE), , drop = FALSE]
log_ok(glue("Submatrix prepared: {nrow(NPP_sub)} genes x {ncol(NPP_sub)} species."))

# ------------------------------
# Block 5b: row order modes (input, cluster, file)
# ------------------------------
read_order_file <- function(path, present_genes) {
  if (!file.exists(path)) log_err(glue("Row-order file not found: {path}"))
  log_info(glue("Reading row order file: {path}"))
  odf <- suppressWarnings(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
  if (nrow(odf) == 0) log_err("Row-order file is empty.")
  order_vec <- odf[[1]] %>% as.character() %>% trimws()
  order_vec <- order_vec[order_vec != "" & !is.na(order_vec)]
  order_vec <- order_vec[order_vec %in% present_genes]
  unique(order_vec)
}

row_order_mode <- args$row_order
desired_row_order <- NULL

if (tolower(row_order_mode) == "input") {
  desired_row_order <- rownames(NPP_sub)
  log_ok("Row order mode: INPUT (using provided gene list order).")
} else if (tolower(row_order_mode) == "cluster") {
  log_ok("Row order mode: CLUSTER (hierarchical clustering).")
} else {
  desired_row_order <- read_order_file(row_order_mode, rownames(NPP_sub))
  if (length(desired_row_order) == 0) log_err("Row-order file yielded empty intersection with gene list.")
  NPP_sub <- NPP_sub[desired_row_order, , drop = FALSE]
  log_ok(glue("Row order mode: FILE (matched {nrow(NPP_sub)} genes from file)."))
}

# ------------------------------
# Block 6: clustering (compute only if mode==cluster and n>=2)
# ------------------------------
row_hclust <- NULL
do_clustering <- FALSE

if (tolower(row_order_mode) == "cluster") {
  if (nrow(NPP_sub) < 2) {
    log_warn("Less than 2 genes -> skipping clustering; using input order.")
  } else {
    row_cor <- suppressWarnings(cor(t(NPP_sub), method = "pearson", use = "pairwise.complete.obs"))
    if (any(is.na(row_cor))) log_warn("Correlation matrix contains NA values; check zero-variance rows.")
    row_dist <- as.dist(1 - row_cor)
    row_hclust <- hclust(row_dist, method = "average")
    do_clustering <- TRUE
    log_ok("Row clustering computed (1 - Pearson; average linkage).")
  }
} else {
  log_ok("Clustering skipped per row_order mode.")
}

# ------------------------------
# Block 7: color mapping for z-scores (z_clip or 1%/99% quantiles)
# ------------------------------
calc_breaks <- function(mat, zclip) {
  v <- as.numeric(mat); v <- v[is.finite(v)]
  if (length(v) == 0) return(c(-1,0,1))
  if (!is.na(zclip)) {
    m <- max(abs(zclip), 1)
    return(c(-m,0,m))
  } else {
    q <- quantile(v, probs = c(0.01,0.99), na.rm = TRUE, names = FALSE)
    m <- max(abs(q)); m <- max(m, 1)
    return(c(-m,0,m))
  }
}
brks <- calc_breaks(NPP_sub, args$z_clip)
col_fun <- colorRamp2(brks, c("#2166AC","#FFFFFF","#B2182B"))

# ------------------------------
# Block 8: assemble heatmap & export (geometry and safe defaults)
# ------------------------------
n_genes <- nrow(NPP_sub)
run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_dir_base <- args$outdir
run_dir <- file.path(run_dir_base, glue("npp_runs_{n_genes}g_{run_stamp}"))
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

out_base <- file.path(run_dir, args$out_prefix)

# geometry
height_cm_per_row <- 0.35
heatmap_body_cm <- n_genes * height_cm_per_row
pad_mm <- c(16,16,24,16)
extra_top_anno_mm <- 3
total_vert_in <- (pad_mm[1] + pad_mm[3] + extra_top_anno_mm) / 25.4
png_width_in <- args$png_width_in
png_height_in <- max(args$min_height_in, heatmap_body_cm / 2.54 + total_vert_in)
log_info(glue("Figure size (in): {png_width_in} x {round(png_height_in,1)}"))

# bottom spacer
bottom_pad <- HeatmapAnnotation(.pad = anno_empty(border = FALSE), annotation_height = unit(18,"mm"))

# prepare heatmap args
hm_args <- list(
  NPP_sub,
  name = "NPP (z)",
  height = unit(heatmap_body_cm, "cm"),
  row_names_gp = gpar(fontsize = 7),
  col = col_fun,
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  top_annotation = clades_annotation,
  bottom_annotation = bottom_pad,
  border = TRUE,
  show_heatmap_legend = TRUE,
  show_row_names = TRUE,
  show_column_names = FALSE
)

if (do_clustering) {
  hm_args$cluster_rows <- as.dendrogram(row_hclust)
} else {
  hm_args$cluster_rows <- FALSE
  hm_args$row_order <- rownames(NPP_sub)
}

hm <- do.call(Heatmap, hm_args)

# Export PNG & PDF
png_file <- glue("{out_base}.png")
pdf_file <- glue("{out_base}.pdf")

# Note: 'type="cairo"' may require Cairo on some Windows setups; warn in README if needed.
png(png_file, width = png_width_in, height = png_height_in, units = "in", res = 300, type = "cairo")
draw(hm, padding = unit(pad_mm, "mm"))
dev.off()

pdf(pdf_file, width = png_width_in, height = png_height_in, onefile = FALSE)
draw(hm, padding = unit(pad_mm, "mm"))
dev.off()

log_ok(glue("Saved PNG: {png_file}"))
log_ok(glue("Saved PDF: {pdf_file}"))

# ------------------------------
# Block 9: export sorted matrix & row order (handle non-cluster case)
# ------------------------------
if (do_clustering && !is.null(row_hclust)) {
  ord <- row_hclust$order
  sorted_genes <- rownames(NPP_sub)[ord]
  NPP_sorted <- NPP_sub[sorted_genes, , drop = FALSE]
} else {
  sorted_genes <- rownames(NPP_sub)
  NPP_sorted <- NPP_sub
}

out_csv <- glue("{out_base}_sorted.csv")
data.table::fwrite(as.data.frame(cbind(gene = rownames(NPP_sorted), NPP_sorted)), out_csv)
writeLines(sorted_genes, glue("{out_base}_row_order.txt"))

log_ok(glue("Saved sorted CSV: {out_csv}"))
log_ok(glue("Saved row order: {out_base}_row_order.txt"))

# ------------------------------
# Block 10: session info & run params (reproducibility)
# ------------------------------
capture.output(sessionInfo(), file = file.path(run_dir, "session_info.txt"))

params_list <- list(
  run_stamp = run_stamp,
  args = args,
  n_genes = n_genes,
  npp_file = args$npp,
  clades_file = args$clades,
  out_base = out_base
)
writeLines(glue::glue_collapse(paste(names(params_list), " : ", lapply(params_list, toString)), sep = "\n"),
           file.path(run_dir, "run_params.txt"))

log_ok("Saved session_info and run_params.")
log_ok(glue("Run complete. Outputs are in: {run_dir}"))

################################################################################
# END
################################################################################
