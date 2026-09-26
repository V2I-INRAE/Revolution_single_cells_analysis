# QC metrics and cell filtering for the REVO samples

suppressPackageStartupMessages({
  library(Seurat)
})

inspect_seurat_qc <- function(seurat_obj) {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[LS]")
  seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)

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

filter_outliers <- function(seurat_obj) {
  cells <- WhichCells(
    seurat_obj,
    expression = nFeature_RNA > 200 & nFeature_RNA <= 5500
  )
  counts <- GetAssayData(
    seurat_obj,
    assay = "RNA", layer = "counts"
  )[, cells, drop = FALSE]
  genes <- rownames(counts)[Matrix::rowSums(counts > 0) > 3]
  clean_seurat_obj <- subset(seurat_obj, cells = cells, features = genes)

  # Filter on mitochondrial fraction of retained genes; keep original QC
  # metadata for comparable before/after plots.
  counts <- GetAssayData(clean_seurat_obj, assay = "RNA", layer = "counts")
  mt_counts <- Matrix::colSums(counts[grep("^MT-", rownames(counts)), , drop = FALSE])
  mt_percent <- 100 * mt_counts / Matrix::colSums(counts)
  mt_max <- mad_bounds(mt_percent, upper_cap = 20)["upper"]
  clean_seurat_obj <- subset(clean_seurat_obj, cells = colnames(counts)[mt_percent <= mt_max])
  cat(sprintf(
    "QC: %d/%d cells, %d/%d genes retained (mt ceiling %.2f%%)\n",
    ncol(clean_seurat_obj), ncol(seurat_obj), nrow(clean_seurat_obj), nrow(seurat_obj), mt_max
  ))

  return(clean_seurat_obj)
}
