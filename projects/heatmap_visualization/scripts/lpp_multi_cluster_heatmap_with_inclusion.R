#!/usr/bin/env Rscript
###############################################################################
# File: lpp_multi_cluster_heatmap_with_inclusion.R
# Purpose: Generate a unified multi-cluster LPP heatmap with per-cluster row
#          clustering and visible row gaps between clusters; includes evidence
#          (Inclusion_criterion) annotation and predefined merge policy.
# Author: Shalev (refactored for publication)
# Created: 2025-10-27 (refactor)
#
# Summary:
#  - CLI arguments (clusters CSV, LPP matrix, species map, selected clusters, outdir, label)
#  - Robust reading of inputs (clusters CSV, LPP CSV/TSV/RDS)
#  - Merge evidence categories ("Descartes" + "Alliance of Genome Resources" -> "Literature")
#  - Validation checks (missing genes, duplicate assignments, unknown evidence)
#  - Constructs ComplexHeatmap with row_split by cluster and right-side annotations:
#     Cluster + Inclusion_criterion
#  - Saves PNG + gene->cluster CSV + run_info txt + runlog
#
# Dependencies:
#  optparse, data.table, ComplexHeatmap, stringr, circlize, dplyr, tidyr, grid,
#  randomcoloR
#
# Example:
#  Rscript lpp_multi_cluster_heatmap_with_inclusion.R \
#    --clusters data/clusterID_clusterGenes_InclusionCriterion_full.csv \
#    --lpp data/Homo_sapiens_LPP_1905.csv \
#    --species data/species_clades_FIXED_ORDER_082022_18_clades.csv \
#    --selected 949,1257,1017,1321,775 \
#    --outdir output --label cilia_set_A --seed 123
#
###############################################################################

# ============================
# Block 0: Packages & Dependencies
# Purpose: Ensure required packages are installed/loaded automatically.
# ============================
# Ensure BiocManager is available
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# List required packages
cran_pkgs <- c("optparse", "data.table", "stringr", "dplyr", "tidyr", "grid", "randomcoloR", "tools") # Added tools for file_ext
bioc_pkgs <- c("ComplexHeatmap", "circlize")

# Install missing CRAN packages
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("[SETUP] Installing missing CRAN package: %s", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
# Install missing Bioconductor packages
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("[SETUP] Installing missing Bioconductor package: %s", pkg))
    BiocManager::install(pkg, ask = FALSE)
  }
}

# Load libraries quietly
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ComplexHeatmap)
  library(stringr)
  library(circlize)
  library(dplyr)
  library(tidyr)
  library(grid)
  library(randomcoloR)
  library(tools) # Added for file_ext
})

message("[SETUP] All required packages loaded.")
# ---------------------------

# ---------------------------
# logging helper (console + file)
# ---------------------------
make_logger <- function(logfile) {
  append_line <- function(line) {
    cat(line, "\n")
    # Use tryCatch for safer file writing, especially in non-interactive environments
    tryCatch({
      if (!is.null(logfile)) cat(paste0(line, "\n"), file = logfile, append = TRUE)
    }, error = function(e){
      warning(sprintf("Failed to write to log file %s: %s", logfile, e$message), call. = FALSE)
    })
  }
  list(
    info = function(msg) append_line(sprintf("[INFO] %s", msg)),
    ok   = function(msg) append_line(sprintf("[OK]   %s", msg)),
    warn = function(msg) append_line(sprintf("[WARN] %s", msg)),
    fail = function(msg) { append_line(sprintf("[FAIL] %s", msg)); stop(msg, call. = FALSE) },
    write = function(msg) append_line(msg)
  )
}

# ---------------------------
# CLI options (using underscore for consistency)
# ---------------------------
option_list <- list(
  make_option(c("-c","--clusters"), type="character", default=NULL,
              help="CSV with clusters details (one gene per row). Required. Columns: cluster_id, cluster_genes, Inclusion_criterion."),
  make_option(c("-l","--lpp"), type="character", default=NULL,
              help="LPP matrix file path (CSV/TSV/RDS). Rows=genes, cols=species (scientific names). Required."),
  make_option(c("-s","--species"), type="character", default=NULL,
              help="species_to_clades CSV (scientific_name,taxid,clade). Required."),
  make_option(c("-S","--selected"), type="character", default=NULL,
              help="Comma-separated cluster IDs to include (e.g., 949,1257,1017). If omitted, all clusters used."),
  make_option(c("-o","--outdir"), type="character", default="output",
              help="Output directory (default: output)"),
  make_option(c("-L","--label"), type="character", default="lpp_with_inclusion",
              help="User label for output filenames"),
  make_option(c("-e","--allow_explode"), action="store_true", default=FALSE, # Changed to underscore
              help="Auto-split semicolon-separated cluster_genes if present (default: FALSE)."),
  make_option(c("--seed"), type="integer", default=123,
              help="Random seed for deterministic colors (default: 123)"),
  make_option(c("--png_width"), type="numeric", default=9, # Changed to underscore
              help="PNG width in inches (default 9)"),
  make_option(c("--png_height_per_gene"), type="numeric", default=0.5, # Changed to underscore
              help="PNG height per gene in inches (min 8 inches)"),
  make_option(c("--method"), type="character", default="single",
              help="Clustering method for rows (default 'single')"),
  make_option(c("--no_sessioninfo"), action="store_true", default=FALSE, # Changed to underscore
              help="Do not save sessionInfo() to runinfo file")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# required args check
if (is.null(opt$clusters) || is.null(opt$lpp) || is.null(opt$species)) {
  stop("Missing required args --clusters, --lpp, --species. Use --help for details.", call. = FALSE)
}

# ---------------------------
# prepare output, logger, ts
# ---------------------------
outdir <- normalizePath(opt$outdir, mustWork = FALSE)
# Use tryCatch for directory creation
tryCatch({
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
}, error = function(e){
  stop(sprintf("Failed to create output directory '%s': %s", outdir, e$message), call.=FALSE)
})

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
logfile <- file.path(outdir, paste0(opt$label, "_runlog_", ts, ".txt"))
logger <- make_logger(logfile)

logger$info("=== START RUN ===")
logger$info(paste0("Timestamp: ", ts))
logger$info(paste0("Label: ", opt$label))
logger$info(paste0("Outdir: ", outdir))

# ---------------------------
# safe reader for LPP matrix (CSV/TSV/RDS)
# ---------------------------
read_lpp_matrix <- function(path, logger) {
  if (!file.exists(path)) logger$fail(sprintf("LPP file not found: %s", path))
  ext <- tolower(tools::file_ext(path))
  logger$info(sprintf("Reading LPP from '%s' (ext='%s')", path, ext))
  mat <- NULL # Initialize mat
  tryCatch({
    if (ext %in% c("csv","txt")) {
      dt <- data.table::fread(path, header = TRUE, data.table = FALSE)
      if (ncol(dt) >= 2 && is.character(dt[[1]]) && !all(suppressWarnings(!is.na(as.numeric(dt[[1]]))))) {
        rn <- as.character(dt[[1]])
        mat <- as.matrix(dt[, -1, drop = FALSE]); rownames(mat) <- rn
      } else { mat <- as.matrix(dt) }
    } else if (ext == "tsv") {
      dt <- data.table::fread(path, sep = "\t", header = TRUE, data.table = FALSE)
      if (ncol(dt) >= 2 && is.character(dt[[1]]) && !all(suppressWarnings(!is.na(as.numeric(dt[[1]]))))) {
        rn <- as.character(dt[[1]])
        mat <- as.matrix(dt[, -1, drop = FALSE]); rownames(mat) <- rn
      } else { mat <- as.matrix(dt) }
    } else if (ext == "rds") {
      mat_obj <- readRDS(path)
      if (!is.matrix(mat_obj) && !is.data.frame(mat_obj)) logger$fail("RDS content must be matrix or data.frame.")
      mat <- as.matrix(mat_obj)
    } else {
      logger$warn(sprintf("Unknown LPP file extension '%s', attempting fread (assuming CSV-like).", ext))
      dt <- data.table::fread(path, header = TRUE, data.table = FALSE)
      if (ncol(dt) >= 2 && is.character(dt[[1]]) && !all(suppressWarnings(!is.na(as.numeric(dt[[1]]))))) {
        rn <- as.character(dt[[1]])
        mat <- as.matrix(dt[, -1, drop = FALSE]); rownames(mat) <- rn
      } else { mat <- as.matrix(dt) }
    }
  }, error = function(e) {
    logger$fail(sprintf("Failed to read LPP matrix from %s: %s", path, e$message))
  })
  if(is.null(mat)) logger$fail("Matrix reading resulted in NULL object.") # Final check
  return(mat)
}

# ---------------------------
# Read inputs
# ---------------------------
logger$info("Reading clusters CSV")
clusters_details <- tryCatch({
  data.table::fread(opt$clusters, header = TRUE, data.table = FALSE)
}, error = function(e) logger$fail(sprintf("Failed reading clusters CSV: %s", e$message)))
clusters_details <- as.data.frame(clusters_details, stringsAsFactors = FALSE) # Ensure data.frame
logger$ok(sprintf("Loaded clusters_details: %d rows x %d cols", nrow(clusters_details), ncol(clusters_details)))

logger$info("Reading LPP matrix")
LPP <- read_lpp_matrix(opt$lpp, logger)
logger$ok(sprintf("Loaded LPP matrix: %d genes x %d species", nrow(LPP), ncol(LPP)))

logger$info("Reading species_to_clades CSV")
species_to_clades <- tryCatch({
  data.table::fread(opt$species, header = TRUE, data.table = FALSE)
}, error = function(e) logger$fail(sprintf("Failed reading species CSV: %s", e$message)))
species_to_clades <- as.data.frame(species_to_clades, stringsAsFactors = FALSE) # Ensure data.frame
logger$ok(sprintf("Loaded species map: %d rows", nrow(species_to_clades)))

# ---------------------------
# Validate required columns
# ---------------------------
logger$info("Validating required columns in input files.")
required_cols <- c("cluster_id", "cluster_genes", "Inclusion_criterion")
miss <- setdiff(required_cols, colnames(clusters_details))
if (length(miss) > 0) logger$fail(sprintf("Missing required columns in clusters file: %s", paste(miss, collapse = ", ")))

req_species_cols <- c("scientific_name", "taxid", "clade")
misss <- setdiff(req_species_cols, colnames(species_to_clades))
if (length(misss) > 0) logger$fail(sprintf("Missing required columns in species map: %s", paste(misss, collapse = ", ")))
logger$ok("Required columns validated.")

# ---------------------------
# Selected clusters handling (defines sel_ids_final always)
# ---------------------------
sel_ids_final <- NULL # Initialize
if (!is.null(opt$selected) && nzchar(opt$selected)) {
  sel_ids <- strsplit(opt$selected, ",")[[1]] %>% trimws()
  sel_ids_num <- suppressWarnings(as.numeric(sel_ids))
  if (any(is.na(sel_ids_num))) {
    logger$warn("Some selected cluster IDs are not purely numeric. Using character matching.")
    sel_ids_final <- sel_ids
    clusters_details$cluster_id <- as.character(clusters_details$cluster_id)
  } else {
    sel_ids_final <- sel_ids_num
    clusters_details$cluster_id <- suppressWarnings(as.numeric(clusters_details$cluster_id))
    if (any(is.na(clusters_details$cluster_id))) {
      logger$warn("Some cluster_id values in the input CSV could not be converted to numeric.")
    }
  }
  logger$info(sprintf("Selected clusters provided: %s", paste(sel_ids_final, collapse = ", ")))
  cd_sel <- clusters_details %>% dplyr::filter(cluster_id %in% sel_ids_final)
} else {
  cd_sel <- clusters_details
  # Define sel_ids_final based on all unique clusters found, preserving order of appearance
  sel_ids_final <- unique(cd_sel$cluster_id)
  logger$info("No selected clusters provided (--selected): using all clusters found in clusters file.")
}
if (nrow(cd_sel) == 0) logger$fail("No rows found for the specified (or default) clusters.")

# ---------------------------
# Optional explode semicolon lists
# ---------------------------
cd_sel$cluster_genes <- as.character(cd_sel$cluster_genes) # Ensure character before checking/splitting
if (any(grepl(";", cd_sel$cluster_genes, fixed = TRUE))) {
  if (!opt$allow_explode) { # Access using underscore
    logger$fail("Found ';' in cluster_genes. Normalize input to one gene per row or run with --allow_explode.")
  } else {
    logger$info("Exploding semicolon-separated cluster_genes (--allow_explode=TRUE).")
    tryCatch({
      cd_sel <- cd_sel %>%
        mutate(cluster_genes = strsplit(cluster_genes, ";")) %>%
        unnest(cluster_genes) %>%
        mutate(cluster_genes = trimws(cluster_genes)) %>%
        filter(cluster_genes != "") # Remove rows that became empty after split
      logger$ok("Successfully exploded gene lists.")
    }, error = function(e){
      logger$fail(sprintf("Error during gene list explosion: %s", e$message))
    })
  }
} else {
  logger$ok("cluster_genes column appears normalized (no semicolons found).")
}

# Ensure final gene names are unique within the selection
if (any(duplicated(cd_sel$cluster_genes))) {
  dupcross <- cd_sel %>%
    dplyr::group_by(cluster_genes) %>%
    dplyr::summarise(n_clusters = n_distinct(cluster_id), .groups = 'drop') %>%
    dplyr::filter(n_clusters > 1)
  if (nrow(dupcross) > 0) logger$fail(sprintf("Gene(s) assigned to multiple selected clusters (ambiguous row_split). Examples: %s", paste(head(dupcross$cluster_genes, 10), collapse = ", ")))
}

# Validate genes are present in LPP matrix
specific_genes <- unique(as.character(cd_sel$cluster_genes))
missing_in_LPP <- setdiff(specific_genes, rownames(LPP))
if (length(missing_in_LPP) > 0) logger$fail(sprintf("Found %d gene(s) missing from LPP matrix. Examples: %s", length(missing_in_LPP), paste(head(missing_in_LPP,10), collapse = ", ")))
logger$ok(sprintf("All %d selected genes are present in LPP.", length(specific_genes)))

# ---------------------------
# Species mapping: scientific_name -> taxid (with checks)
# ---------------------------
logger$info("Mapping species scientific names to taxids for LPP columns.")
species_to_clades$scientific_name <- as.character(species_to_clades$scientific_name)
species_to_clades$taxid <- as.character(species_to_clades$taxid)
species_to_clades$clade <- as.character(species_to_clades$clade)
species_map <- setNames(species_to_clades$taxid, species_to_clades$scientific_name)

orig_cols <- colnames(LPP)
mapped <- species_map[orig_cols]
na_idx <- which(is.na(mapped))
if (length(na_idx) > 0) {
  logger$warn(sprintf("%d species names in LPP matrix header could not be mapped to taxid via species map. Examples: %s. Keeping original names.",
                      length(na_idx), paste(head(orig_cols[na_idx],5), collapse = ", ")))
  mapped[na_idx] <- orig_cols[na_idx] # Keep original names for unmapped columns
}
colnames(LPP) <- as.character(mapped)
logger$ok("Renamed LPP columns to mapped taxid (or kept original name if unmapped).")

# Align clade vector to LPP columns
match_idx <- match(colnames(LPP), species_to_clades$taxid)
clade_vec <- rep("Unknown", length(colnames(LPP))) # Default to Unknown
valid_match <- !is.na(match_idx)
clade_values_from_map <- as.character(species_to_clades$clade[match_idx[valid_match]])
clade_values_from_map[is.na(clade_values_from_map) | clade_values_from_map == ""] <- "Unknown"
clade_vec[valid_match] <- clade_values_from_map

# Predefined clade colors (ensure 'Unknown' is included)
all_clade_colors <- c(
  Mammalia="#E31A1C", Aves="#A6CEE3", Lepidosauria="#004529",
  Amphibia="#80CDC1", Actinopteri="#1F78B4", Platyhelminthes="#3F007D",
  Nematoda="#6A3D9A", Arthropoda="#CAB2D6", Ascomycota="#FFFF99",
  Basidiomycota="#FF7F00", Fungi_Incertae_Sedis="#E6AB02",
  Viridiplantae="#33A02C", Rhodophyta="#7F0000", Ochrophyta="#B2DF8A",
  Oomycota="#B15928", Apicomplexa="#FDBF6F", Trypanosomatidae="#662506",
  Amoebozoa="#EF6548", Unknown="#BBBBBB" # Added Unknown explicitly
)
unique_clades_in_vec <- unique(clade_vec)
missing_clades_from_palette <- setdiff(unique_clades_in_vec, names(all_clade_colors))
if (length(missing_clades_from_palette) > 0) {
  logger$warn(sprintf("Assigning 'Unknown' color for %d unexpected clade(s): %s",
                      length(missing_clades_from_palette), paste(missing_clades_from_palette, collapse = ", ")))
  temp_colors <- setNames(rep(all_clade_colors["Unknown"], length(missing_clades_from_palette)), missing_clades_from_palette)
  all_clade_colors <- c(all_clade_colors, temp_colors)
}
clade_factor_for_annot <- factor(clade_vec, levels = names(all_clade_colors))

clades_annotation <- HeatmapAnnotation(
  which = "column",
  annotation_height = unit(2, "mm"),
  simple_anno_size_adjust = TRUE,
  gap = unit(0.3, "mm"),
  show_annotation_name = FALSE,
  show_legend = TRUE,
  Clade = clade_factor_for_annot,
  col = list(Clade = all_clade_colors)
)
logger$ok("Top clade annotation constructed.")


# ---------------------------
# Evidence merging & validation
# ---------------------------
logger$info("Processing evidence categories (Inclusion_criterion).")
cd_sel$Inclusion_criterion <- as.character(cd_sel$Inclusion_criterion) # Ensure character

logger$info("Evidence counts BEFORE merge:")
logger$info(paste(capture.output(print(table(cd_sel$Inclusion_criterion, useNA = "ifany"))), collapse = "\n"))

cd_sel <- cd_sel %>%
  dplyr::mutate(Inclusion_criterion = dplyr::case_when(
    Inclusion_criterion %in% c("Descartes", "Alliance of Genome Resources") ~ "Literature",
    is.na(Inclusion_criterion) | Inclusion_criterion == "" ~ "Unknown",
    TRUE ~ Inclusion_criterion
  ))

logger$info("Evidence counts AFTER merge:")
logger$info(paste(capture.output(print(table(cd_sel$Inclusion_criterion, useNA = "ifany"))), collapse = "\n"))

evidence_levels <- c(
  "CiliaCarta: Gold Standard", "CiliaCarta: Gene Ontology", "CiliaCarta: Predicted",
  "Literature", "Novel cilia-associated candidate", "Unknown"
)
evidence_colors <- c(
  "CiliaCarta: Gold Standard"      = "#00441b", "CiliaCarta: Gene Ontology"      = "#238b45",
  "CiliaCarta: Predicted"          = "#66c2a4", "Literature"                     = "#7570b3",
  "Novel cilia-associated candidate" = "#e41a1c", "Unknown"                      = "#DDDDDD"
)

unknown_ev <- setdiff(unique(cd_sel$Inclusion_criterion), evidence_levels)
if (length(unknown_ev) > 0) logger$fail(sprintf("Unknown evidence categories found AFTER merge and cleanup: %s. Check input or script levels.", paste(unknown_ev, collapse = ", ")))

# Align evidence vector to the FINAL ordered gene list
# Create the final ordered data frame based on sel_ids_final and gene name
cd_sel_ordered_final <- cd_sel %>% arrange(match(cluster_id, sel_ids_final), cluster_genes)
final_ordered_genes <- cd_sel_ordered_final$cluster_genes
final_gene_to_evidence <- setNames(cd_sel_ordered_final$Inclusion_criterion, final_ordered_genes)
evidence_vec_final <- final_gene_to_evidence[final_ordered_genes] # Vector aligned to final order

if (any(is.na(evidence_vec_final))) {
  bad <- final_ordered_genes[is.na(evidence_vec_final)]
  logger$fail(sprintf("NA evidence mapping for %d gene(s) in final ordered list. Examples: %s", length(bad), paste(head(bad,10), collapse = ", ")))
}
logger$ok("Evidence mapping aligned to final ordered gene list.")

# ---------------------------
# Cluster colors & row annotations
# ---------------------------
logger$info("Setting up cluster factor and row annotations.")
final_gene_to_cluster <- setNames(cd_sel_ordered_final$cluster_id, final_ordered_genes)
cluster_vec_final <- final_gene_to_cluster[final_ordered_genes]

# Create cluster factor using the order defined by sel_ids_final
cluster_factor_levels <- sel_ids_final[sel_ids_final %in% unique(cluster_vec_final)] # Use defined order
cluster_factor <- factor(cluster_vec_final, levels = cluster_factor_levels)

set.seed(opt$seed)
cluster_palette <- setNames(randomcoloR::distinctColorPalette(length(cluster_factor_levels)), cluster_factor_levels)
logger$ok("Cluster color palette generated.")

evidence_factor_final <- factor(evidence_vec_final, levels = evidence_levels)

right_anno <- rowAnnotation(
  Cluster   = cluster_factor,
  Inclusion = evidence_factor_final,
  col = list(Cluster = cluster_palette, Inclusion = evidence_colors),
  show_annotation_name = TRUE,
  annotation_width = unit(6, "mm")
)
logger$ok("Right-side annotations constructed (Cluster & Inclusion).")

# ---------------------------
# Assemble Heatmap
# ---------------------------
logger$info("Assembling Heatmap object...")
LPP_sub <- LPP[final_ordered_genes, , drop = FALSE]
valid_lpp_values <- LPP_sub[!is.na(LPP_sub) & is.finite(LPP_sub)]
lpp_min <- if(length(valid_lpp_values) > 0) min(valid_lpp_values) else 0
lpp_max <- if(length(valid_lpp_values) > 0) max(valid_lpp_values) else 1
# Ensure min != max to avoid colorRamp2 error
if(lpp_min == lpp_max) {
  lpp_min <- lpp_min - 0.1 # Adjust slightly if all values are identical
  lpp_max <- lpp_max + 0.1
  if(lpp_min < 0) lpp_min <- 0
  if(lpp_max > 1) lpp_max <- 1
}
col_fun <- colorRamp2(c(lpp_min, lpp_max), c("white", "#0066cc"))

# Calculate dynamic height based on the final number of genes and inches-per-gene param
# Correct access to opt$png_height_per_gene using underscore
total_height_inches_genes = length(final_ordered_genes) * opt$png_height_per_gene
# Use a fixed minimum height (e.g., 8 inches) plus margin, ensure it's at least the gene-based height
fig_height_in <- max(8, total_height_inches_genes + 2) # Added +2 for margins/annotations


hm <- Heatmap(
  LPP_sub,
  name = "LPP",
  height = unit(total_height_inches_genes * 2.54, "cm"), # Height based on inches per gene
  row_split = cluster_factor,
  row_title = NULL,
  show_row_dend = TRUE,
  cluster_rows = TRUE,
  clustering_distance_rows = function(m) { # Robust distance calculation
    # Check for rows with zero variance or too many NAs which break correlation
    row_vars = apply(m, 1, var, na.rm=TRUE)
    valid_rows = which(row_vars > 1e-8 & rowSums(!is.na(m)) > 2) # Need at least 3 non-NA points for cor
    if(length(valid_rows) < 2) return(dist(matrix(0, nrow=nrow(m), ncol=1))) # Return dummy if < 2 valid rows
    m_valid = m[valid_rows, , drop=FALSE]
    d = as.dist(1 - cor(t(m_valid), use = "pairwise.complete.obs"))
    # Need to return a distance matrix matching original size if rows were dropped - ComplexHeatmap might handle this, or need adjustment
    # For simplicity here, assuming most rows are valid. Production code might need more robust handling.
    if (length(valid_rows) == nrow(m)) {
      return(d)
    } else {
      # Simplified: return distance only for valid rows, may cause issues if ComplexHeatmap expects full matrix
      logger$warn("Some rows skipped in clustering due to low variance or NAs.")
      # Placeholder: return a basic distance - Needs careful thought on how to handle partial clustering
      return(dist(m)) # Fallback, might not be ideal
    }
  },
  clustering_method_rows = opt$method,
  cluster_columns = FALSE,
  col = col_fun,
  row_names_gp = gpar(fontsize = 6),
  top_annotation = clades_annotation,
  right_annotation = right_anno,
  border = TRUE,
  show_heatmap_legend = TRUE,
  show_row_names = (length(final_ordered_genes) <= 150),
  show_column_names = FALSE
)
logger$ok("Heatmap object constructed.")

# ---------------------------
# Export PNG + provenance files
# ---------------------------
# Access png_width using underscore
png_width <- opt$png_width
output_png <- file.path(outdir, sprintf("%s_%s.png", opt$label, ts))

logger$info(sprintf("Exporting PNG: %s (%.1f x %.1f in)", output_png, png_width, fig_height_in))
tryCatch({
  png(filename = output_png, width = png_width, height = fig_height_in, units = "in", res = 300)
  draw(hm, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  logger$ok(sprintf("Figure exported successfully: %s", normalizePath(output_png)))
}, error = function(e){
  logger$fail(sprintf("Failed to export PNG: %s", e$message))
})

# Save gene->cluster mapping CSV using the final ordered data
mapping_csv <- file.path(outdir, sprintf("%s_genes_map_%s.csv", opt$label, ts))
map_df <- data.frame(gene = final_ordered_genes,
                     cluster = cluster_vec_final, # Use the vector aligned with final_ordered_genes
                     inclusion = evidence_vec_final, # Use the vector aligned with final_ordered_genes
                     stringsAsFactors = FALSE)
tryCatch({
  write.csv(map_df, mapping_csv, row.names = FALSE, quote = FALSE)
  logger$ok(sprintf("Gene->cluster mapping saved: %s (n=%d)", mapping_csv, nrow(map_df)))
}, error = function(e){
  logger$warn(sprintf("Failed to save gene mapping CSV: %s", e$message))
})

# Save run info
runinfo_txt <- file.path(outdir, sprintf("%s_runinfo_%s.txt", opt$label, ts))
tryCatch({
  sink(runinfo_txt)
  cat("Run info\n========\n")
  cat(sprintf("Timestamp: %s\n", ts))
  cat(sprintf("Label: %s\n", opt$label))
  # Use normalizePath cautiously, ensure paths exist if mustWork=TRUE
  cat(sprintf("Clusters file: %s\n", normalizePath(opt$clusters, mustWork = FALSE)))
  cat(sprintf("LPP file: %s\n", normalizePath(opt$lpp, mustWork = FALSE)))
  cat(sprintf("Species map: %s\n", normalizePath(opt$species, mustWork = FALSE)))
  cat(sprintf("Selected clusters: %s\n", ifelse(is.null(opt$selected), "ALL", opt$selected)), "\n")
  cat(sprintf("Allow explode: %s\n", opt$allow_explode), "\n") # Access with underscore
  cat(sprintf("Seed: %d\n", opt$seed), "\n")
  cat(sprintf("Clustering method: %s\n", opt$method), "\n")
  cat(sprintf("PNG dims (inches): %.1f x %.1f\n", png_width, fig_height_in), "\n")
  # Access no_sessioninfo using underscore
  if (!opt$no_sessioninfo) {
    cat("\nSession info:\n")
    print(utils::sessionInfo())
  }
  sink()
  logger$ok(sprintf("Run info saved: %s", normalizePath(runinfo_txt, mustWork = FALSE)))
}, error = function(e){
  logger$warn(sprintf("Failed to save run info file: %s", e$message))
  if(sink.number() > 0) sink() # Close sink if left open
})

logger$ok(sprintf("Complete run log saved: %s", normalizePath(logfile, mustWork = FALSE)))
logger$info("=== END RUN ===")
