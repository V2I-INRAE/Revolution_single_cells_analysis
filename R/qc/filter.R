# QC metrics, outlier labeling and cell filtering for the REVO samples

suppressPackageStartupMessages({
  library(Seurat)
})

inspect_seurat_qc <- function(seurat_obj) {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
    seurat_obj,
    pattern = qc_params$metrics$mitochondrial_pattern
  )
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(
    seurat_obj,
    pattern = qc_params$metrics$ribosomal_pattern
  )
  ribo_threshold <- qc_params$metrics$min_ribo_percent
  seurat_obj$ribo_status <- ifelse(
    seurat_obj$percent.ribo < ribo_threshold,
    sprintf("< %.1f%%", ribo_threshold),
    sprintf(">= %.1f%%", ribo_threshold)
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
    log10GenesPerUMI = qs(seurat_obj$log10GenesPerUMI)
  )
  colnames(summary_tbl) <- c("min", "q25", "median", "q75", "max")
  cat("\nQC metric summary:\n")
  print(summary_tbl)
  cat("\nRibosomal percentage labels:\n")
  print(table(seurat_obj$ribo_status))

  return(seurat_obj)
}

# --- Per-cell QC values, appended across samples in main.R and used by
# plot_qc_comparison. `stage` is "before" or "after" outlier filtering.
collect_qc_summary <- function(
  seurat_obj,
  sample_id,
  stage,
  coordinates,
  cells = colnames(seurat_obj)
) {
  stopifnot(
    all(c("cell", "UMAP1", "UMAP2") %in% colnames(coordinates)),
    anyDuplicated(coordinates$cell) == 0
  )
  metadata <- seurat_obj[[]][cells, , drop = FALSE]
  coordinate_idx <- match(rownames(metadata), coordinates$cell)
  stopifnot(!anyNA(coordinate_idx))

  data.frame(
    sample_id = sample_id,
    cell = rownames(metadata),
    stage = stage,
    UMAP1 = coordinates$UMAP1[coordinate_idx],
    UMAP2 = coordinates$UMAP2[coordinate_idx],
    nFeature_RNA = metadata$nFeature_RNA,
    nCount_RNA = metadata$nCount_RNA,
    percent.mt = metadata$percent.mt,
    percent.ribo = metadata$percent.ribo,
    ribo_status = metadata$ribo_status,
    log10GenesPerUMI = metadata$log10GenesPerUMI,
    row.names = NULL
  )
}

mad_bounds <- function(x, n_mad, lower_floor = -Inf, upper_cap = Inf) {
  med <- median(x)
  mad_x <- mad(x)
  c(
    lower = max(lower_floor, med - n_mad * mad_x),
    upper = min(upper_cap, med + n_mad * mad_x)
  )
}

label_qc_outliers <- function(seurat_obj) {
  params <- qc_params$filtering

  feature_max <- mad_bounds(
    seurat_obj$nFeature_RNA,
    upper_cap = params$max_features_cap,
    n_mad = params$feature_mad_multiplier
  )["upper"]

  seurat_obj$flag_low_features <- (
    seurat_obj$nFeature_RNA <= params$min_features
  )
  seurat_obj$flag_high_features <- seurat_obj$nFeature_RNA > feature_max
  seurat_obj$flag_low_complexity <- (
    seurat_obj$log10GenesPerUMI <= params$min_log10_genes_per_umi
  )

  initial_outlier <- (
    seurat_obj$flag_low_features |
      seurat_obj$flag_high_features |
      seurat_obj$flag_low_complexity
  )

  mt_percent <- seurat_obj$percent.mt[!initial_outlier]
  mt_max <- mad_bounds(
    mt_percent,
    upper_cap = params$max_mito_percent_cap,
    n_mad = params$mito_mad_multiplier
  )["upper"]
  seurat_obj$flag_high_mt <- seurat_obj$percent.mt > mt_max

  flag_cols <- c(
    "flag_low_features", "flag_high_features",
    "flag_low_complexity", "flag_high_mt"
  )
  seurat_obj$qc_outlier <- rowSums(seurat_obj[[]][, flag_cols]) > 0

  cat("\nFlagged cells:\n")
  for (col in c(flag_cols, "qc_outlier")) {
    flag <- seurat_obj[[col]][[1]]
    cat(sprintf(
      "  %-22s %6d (%5.1f%%)\n",
      col, sum(flag), 100 * mean(flag)
    ))
  }
  cat(sprintf("  nFeature_RNA ceiling %.2f\n", feature_max))
  cat(sprintf("  mitochondrial ceiling %.2f%%\n", mt_max))

  return(seurat_obj)
}

filter_labeled_cells <- function(seurat_obj) {
  stopifnot(all(
    c("qc_outlier", "doublet_class") %in% colnames(seurat_obj[[]])
  ))

  is_doublet <- seurat_obj$doublet_class == "Doublet"
  seurat_obj$keep_cell <- !seurat_obj$qc_outlier & !is_doublet
  seurat_obj$exclusion_reason <- "retained"
  seurat_obj$exclusion_reason[
    seurat_obj$qc_outlier & !is_doublet
  ] <- "qc_outlier"
  seurat_obj$exclusion_reason[
    !seurat_obj$qc_outlier & is_doublet
  ] <- "doublet"
  seurat_obj$exclusion_reason[seurat_obj$qc_outlier & is_doublet] <- (
    "qc_outlier+doublet"
  )

  selected_cells <- colnames(seurat_obj)[seurat_obj$keep_cell]
  counts <- GetAssayData(
    seurat_obj,
    assay = "RNA", layer = "counts"
  )[, selected_cells, drop = FALSE]
  selected_features <- rownames(counts)[
    Matrix::rowSums(counts > 0) >= qc_params$filtering$min_cells_per_feature
  ]

  clean_seurat_obj <- subset(
    seurat_obj,
    cells = selected_cells,
    features = selected_features
  )

  cat("\nCombined filtering:\n")
  for (reason in c("qc_outlier", "doublet", "qc_outlier+doublet")) {
    n_removed <- sum(seurat_obj$exclusion_reason == reason)
    cat(sprintf(
      "  %-22s %6d (%5.1f%%)\n",
      reason,
      n_removed,
      100 * n_removed / ncol(seurat_obj)
    ))
  }
  cat(sprintf(
    "  cells retained         %6d (%5.1f%%)\n",
    ncol(clean_seurat_obj), 100 * ncol(clean_seurat_obj) / ncol(seurat_obj)
  ))
  cat(sprintf(
    "  genes retained         %6d (%5.1f%%)\n",
    nrow(clean_seurat_obj), 100 * nrow(clean_seurat_obj) / nrow(seurat_obj)
  ))

  return(clean_seurat_obj)
}
