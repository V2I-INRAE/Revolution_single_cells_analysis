# QC metrics, outlier labeling and cell filtering.

suppressPackageStartupMessages({
  library(Seurat)
})

inspect_seurat_qc <- function(seurat_obj, qc_config) {
  metrics <- qc_config$metrics
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
    seurat_obj,
    pattern = metrics$mitochondrial_pattern
  )
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(
    seurat_obj,
    pattern = metrics$ribosomal_pattern
  )
  # percent.hb from the explicit hemoglobin gene list, never a ^HB pattern
  # (HBEGF, HBP1 and HBS1L match but are not globins). Unlike the
  # mitochondrial list, which io.R stopifnots because it rewrites those
  # rownames, partially absent genes are dropped with a note: percent.hb
  # stays a valid fraction of the remaining globins, and the pilot
  # pre-rebuild matrices contain all five, so a drop would flag a changed
  # input source. Total absence still fails: an empty selection would return
  # zeros for every cell and silently disable the hemoglobin filter.
  hb_genes <- intersect(metrics$hemoglobin_genes, rownames(seurat_obj))
  if (length(hb_genes) == 0L) {
    stop("None of the configured hemoglobin genes are in the matrix: ",
      paste(metrics$hemoglobin_genes, collapse = ", "))
  }
  missing_hb <- setdiff(metrics$hemoglobin_genes, hb_genes)
  if (length(missing_hb) > 0L) {
    message("Hemoglobin genes not in the matrix, excluded from percent.hb: ",
      paste(missing_hb, collapse = ", "))
  }
  seurat_obj[["percent.hb"]] <- PercentageFeatureSet(
    seurat_obj,
    features = hb_genes
  )
  seurat_obj$log10GenesPerUMI <- (
    log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
  )

  qs <- function(x) round(quantile(x, c(0, 0.25, 0.5, 0.75, 1)), 2)
  summary_tbl <- rbind(
    nCount_RNA = qs(seurat_obj$nCount_RNA),
    nFeature_RNA = qs(seurat_obj$nFeature_RNA),
    percent.mt = qs(seurat_obj$percent.mt),
    percent.ribo = qs(seurat_obj$percent.ribo),
    percent.hb = qs(seurat_obj$percent.hb),
    log10GenesPerUMI = qs(seurat_obj$log10GenesPerUMI)
  )
  colnames(summary_tbl) <- c("min", "q25", "median", "q75", "max")
  cat("\nQC metric summary:\n")
  print(summary_tbl)

  return(seurat_obj)
}

label_qc_outliers <- function(seurat_obj, filtering_config) {
  # Adaptive limits computed per sample: median + k * MAD, capped by the
  # fixed ceilings; cells at either limit pass.
  feature_max <- min(
    filtering_config$max_features_cap,
    median(seurat_obj$nFeature_RNA) +
      filtering_config$feature_mads * mad(seurat_obj$nFeature_RNA)
  )
  mt_max <- min(
    filtering_config$max_mito_cap,
    median(seurat_obj$percent.mt) +
      filtering_config$mito_mads * mad(seurat_obj$percent.mt)
  )
  ribo_min <- filtering_config$min_ribo_percent

  seurat_obj$flag_low_features <- seurat_obj$nFeature_RNA < filtering_config$min_features
  seurat_obj$flag_high_features <- seurat_obj$nFeature_RNA > feature_max
  seurat_obj$flag_low_complexity <- (
    seurat_obj$log10GenesPerUMI <= filtering_config$min_log10_genes_per_umi
  )
  seurat_obj$flag_high_mt <- seurat_obj$percent.mt > mt_max
  seurat_obj$flag_low_ribo <- seurat_obj$percent.ribo <= ribo_min
  # Deliberate difference from every other ceiling: hemoglobin contamination is
  # excluded at or above the fixed threshold (>=), per the run decision from
  # the pilot; all other limits exclude only strictly-above cells.
  seurat_obj$flag_high_hb <- (
    seurat_obj$percent.hb >= filtering_config$max_hb_percent
  )

  flag_cols <- c(
    "flag_low_features", "flag_high_features",
    "flag_low_complexity", "flag_high_mt", "flag_high_hb", "flag_low_ribo"
  )
  seurat_obj$qc_outlier <- (
    seurat_obj$flag_low_features | seurat_obj$flag_high_features |
      seurat_obj$flag_low_complexity | seurat_obj$flag_high_mt |
      seurat_obj$flag_high_hb | seurat_obj$flag_low_ribo
  )

  cat("\nFlagged cells:\n")
  for (col in c(flag_cols, "qc_outlier")) {
    flag <- seurat_obj[[col]][[1]]
    cat(sprintf("  %-22s %6d (%5.1f%%)\n", col, sum(flag), 100 * mean(flag)))
  }
  cat(sprintf(
    "  nFeature_RNA ceiling %.2f (median + %g MAD, cap %g)\n",
    feature_max, filtering_config$feature_mads, filtering_config$max_features_cap
  ))
  cat(sprintf(
    "  mitochondrial ceiling %.2f%% (median + %g MAD, cap %g%%)\n",
    mt_max, filtering_config$mito_mads, filtering_config$max_mito_cap
  ))
  cat(sprintf(
    "  hemoglobin ceiling %.2f%% (cells at or above excluded)\n",
    filtering_config$max_hb_percent
  ))
  cat(sprintf("  ribosomal floor %.2f%% (cells above retained)\n", ribo_min))

  seurat_obj@misc$qc <- list(
    feature_max = feature_max, mt_max = mt_max, ribo_min = ribo_min
  )
  return(seurat_obj)
}

label_retained_cells <- function(seurat_obj) {
  seurat_obj$keep_qc <- !seurat_obj$qc_outlier
  seurat_obj$keep_scrublet <- seurat_obj$keep_qc & !seurat_obj$predicted_doublets
  return(seurat_obj)
}

# Subset to Scrublet-clean cells and reapply the per-feature minimum, keeping
# per-sample count layers and the original (unrecomputed) metadata values.
filter_labeled_cells <- function(seurat_obj, qc_config) {
  clean <- subset(seurat_obj, subset = keep_scrublet)

  min_cells <- qc_config$filtering$min_cells_per_feature
  samples <- seurat_obj@misc$qc$sample_order
  metadata <- seurat_obj[[]]
  metadata$qc_umap_1 <- NULL
  metadata$qc_umap_2 <- NULL
  objects <- list()

  for (sample in samples) {
    counts <- LayerData(clean, assay = "RNA", layer = paste0("counts.", sample))
    obj <- CreateSeuratObject(counts, project = sample, min.cells = min_cells)
    meta <- metadata[colnames(obj), , drop = FALSE]
    obj <- AddMetaData(obj, meta)
    obj$keep_cell <- TRUE
    obj$doublet_filter_method <- "Scrublet"

    message(sample, " / Scrublet: ", ncol(obj), " cells; ", nrow(obj), " genes")
    objects[[sample]] <- obj
  }
  clean <- merge(objects[[1]], y = objects[-1], merge.data = FALSE)
  clean@misc$qc <- seurat_obj@misc$qc

  return(clean)
}

# Drop the configured feature names (config$qc$remove_genes) from the matrix.
# Exact membership (%in%), deliberately not grepl: grepl vectorizes over x,
# the strings being tested, not over pattern, so a multi-element pattern is
# recycled element-wise and silently mismatches; the regex alternative,
# paste(patterns, collapse = "|"), would still risk substring hits (e.g.
# MT-ND4 would also catch MT-ND4L). Every entry is a full feature name, so
# %in% is simpler and exact. Listed names absent from the matrix are
# tolerated and logged: a future input source may lack one. The removed set
# is stored in misc$qc so the run record can state the removal.
remove_listed_genes <- function(seurat_obj, qc_config) {
  remove_genes <- qc_config$remove_genes
  hits <- rownames(seurat_obj) %in% remove_genes
  removed <- rownames(seurat_obj)[hits]
  seurat_obj <- seurat_obj[!hits, ]
  seurat_obj@misc$qc$removed_genes <- removed
  message("Removed ", length(removed), " listed gene(s): ",
    paste(removed, collapse = ", "))
  missing <- setdiff(remove_genes, removed)
  if (length(missing) > 0L) {
    message("Listed gene(s) not found in the matrix, nothing to remove: ",
      paste(missing, collapse = ", "))
  }
  return(seurat_obj)
}

# Minimal handoff contract between QC and normalization: unique aligned cell
# IDs, RNA counts present and nothing normalized or scaled yet.
validate_clean_object <- function(seurat_obj) {
  DefaultAssay(seurat_obj) <- "RNA"
  rna_layers <- Layers(seurat_obj[["RNA"]])
  stopifnot(
    inherits(seurat_obj, "Seurat"),
    "RNA" %in% Assays(seurat_obj),
    anyDuplicated(colnames(seurat_obj)) == 0,
    identical(rownames(seurat_obj[[]]), colnames(seurat_obj)),
    "RNA counts are required" = any(grepl("^counts($|\\.)", rna_layers)),
    "Expected counts-only RNA input" = !any(
      grepl("^(data|scale\\.data)($|\\.)", rna_layers)
    )
  )
  invisible(TRUE)
}
