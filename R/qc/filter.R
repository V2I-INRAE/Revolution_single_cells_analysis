# QC metrics and outlier labeling for the REVO samples

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
collect_qc_summary <- function(seurat_obj, sample_id, stage) {
  data.frame(
    sample_id = sample_id,
    cell = colnames(seurat_obj),
    stage = stage,
    nFeature_RNA = seurat_obj$nFeature_RNA,
    nCount_RNA = seurat_obj$nCount_RNA,
    percent.mt = seurat_obj$percent.mt,
    percent.ribo = seurat_obj$percent.ribo,
    log10GenesPerUMI = seurat_obj$log10GenesPerUMI,
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
  seurat_obj$flag_low_features <- seurat_obj$nFeature_RNA <= 200
  seurat_obj$flag_high_features <- seurat_obj$nFeature_RNA > 5500
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
  cat(sprintf("  mitochondrial ceiling %.2f%%\n", mt_max))

  return(seurat_obj)
}
