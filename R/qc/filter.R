# QC metrics, outlier labeling and cell filtering for the REVO samples

suppressPackageStartupMessages({
  library(Seurat)
})

inspect_seurat_qc <- function(seurat_obj) {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
    seurat_obj,
    pattern = "^MT-"
  )
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(
    seurat_obj,
    pattern = "^RP[LS]"
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

  return(seurat_obj)
}

# --- Per-cell QC values, appended across samples in main.R and used by
# plot_qc_comparison. `stage` is "before" or "after" outlier filtering.
collect_qc_summary <- function(
  seurat_obj,
  sample_id,
  stage,
  cells = colnames(seurat_obj)
) {
  metadata <- seurat_obj[[]][cells, , drop = FALSE]

  data.frame(
    sample_id = sample_id,
    cell = rownames(metadata),
    stage = stage,
    nFeature_RNA = metadata$nFeature_RNA,
    nCount_RNA = metadata$nCount_RNA,
    percent.mt = metadata$percent.mt,
    percent.ribo = metadata$percent.ribo,
    log10GenesPerUMI = metadata$log10GenesPerUMI,
    row.names = NULL
  )
}

mad_bounds <- function(x, lower_floor = -Inf, upper_cap = Inf, n_mad = 5) {
  med <- median(x)
  mad_x <- mad(x)
  c(
    lower = max(lower_floor, med - n_mad * mad_x),
    upper = min(upper_cap, med + n_mad * mad_x)
  )
}

label_qc_outliers <- function(seurat_obj) {
  feature_max <- mad_bounds(
    seurat_obj$nFeature_RNA,
    upper_cap = 5000,
    n_mad = 4
  )["upper"]

  seurat_obj$flag_low_features <- seurat_obj$nFeature_RNA <= 200
  seurat_obj$flag_high_features <- seurat_obj$nFeature_RNA > feature_max
  seurat_obj$flag_low_complexity <- seurat_obj$log10GenesPerUMI <= 0.8

  initial_outlier <- (
    seurat_obj$flag_low_features |
      seurat_obj$flag_high_features |
      seurat_obj$flag_low_complexity
  )
  mt_percent <- seurat_obj$percent.mt[!initial_outlier]
  mt_max <- mad_bounds(mt_percent, upper_cap = 20)["upper"]
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
  selected_features <- rownames(counts)[Matrix::rowSums(counts > 0) > 3]

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
