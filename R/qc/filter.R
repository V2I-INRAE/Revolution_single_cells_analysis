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

label_qc_outliers <- function(seurat_obj) {
  params <- qc_params$filtering

  # Fixed limits shared by every sample; cells at either upper limit pass.
  feature_max <- params$max_features

  seurat_obj$flag_low_features <- (
    seurat_obj$nFeature_RNA < params$min_features
  )
  seurat_obj$flag_high_features <- seurat_obj$nFeature_RNA > feature_max
  seurat_obj$flag_low_complexity <- (
    seurat_obj$log10GenesPerUMI <= params$min_log10_genes_per_umi
  )

  mt_max <- params$max_mito_percent
  seurat_obj$flag_high_mt <- seurat_obj$percent.mt > mt_max

  flag_cols <- c(
    "flag_low_features", "flag_high_features",
    "flag_low_complexity", "flag_high_mt"
  )
  seurat_obj$qc_outlier <- (
    seurat_obj$flag_low_features | seurat_obj$flag_high_features |
      seurat_obj$flag_low_complexity | seurat_obj$flag_high_mt
  )

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

  seurat_obj@misc$qc <- list(
    feature_max = feature_max, mt_max = mt_max
  )
  return(seurat_obj)
}

label_retained_cells <- function(seurat_obj) {
  seurat_obj$keep_qc <- !seurat_obj$qc_outlier
  seurat_obj$keep_scrublet <- (
    seurat_obj$keep_qc & !seurat_obj$predicted_doublets
  )
  return(seurat_obj)
}

filter_labeled_cells <- function(seurat_obj) {
  clean <- subset(seurat_obj, subset = keep_scrublet)

  min_cells <- seurat_obj@misc$qc$params$filtering$min_cells_per_feature
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

    message(
      sample,
      " / Scrublet",
      ": ",
      ncol(obj),
      " cells; ",
      nrow(obj),
      " genes"
    )
    objects[[sample]] <- obj
  }
  clean <- merge(objects[[1]], y = objects[-1], merge.data = FALSE)
  clean@misc$qc <- seurat_obj@misc$qc

  return(clean)
}
