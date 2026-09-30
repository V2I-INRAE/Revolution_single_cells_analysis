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

  seurat_obj@misc$qc <- list(
    feature_max = unname(feature_max), mt_max = unname(mt_max)
  )
  return(seurat_obj)
}

label_retained_cells <- function(seurat_obj) {
  seurat_obj$keep_qc <- !seurat_obj$qc_outlier
  seurat_obj$keep_doubletfinder <- (
    seurat_obj$keep_qc & seurat_obj$doublet_class == "Singlet"
  )
  seurat_obj$keep_scrublet <- (
    seurat_obj$keep_qc & !seurat_obj$predicted_doublets
  )
  return(seurat_obj)
}

# Filter the checkpoint layer by layer: a pooled gene filter changes the method.
filter_labeled_cells <- function(seurat_obj, doublet_method) {
  doublet_method <- match.arg(doublet_method, c("DoubletFinder", "Scrublet"))
  keep_col <- paste0("keep_", tolower(doublet_method))
  metadata <- seurat_obj[[]]
  min_cells <- seurat_obj@misc$qc$params$filtering$min_cells_per_feature
  samples <- seurat_obj@misc$qc$sample_order
  objects <- setNames(lapply(samples, function(sample) {
    cells <- rownames(metadata)[metadata$sample == sample & metadata[[keep_col]]]
    stopifnot("No retained cells in sample" = length(cells) > 0)
    counts <- LayerData(seurat_obj, assay = "RNA", layer = paste0("counts.", sample))
    stopifnot(all(cells %in% colnames(counts)))
    counts <- counts[, cells, drop = FALSE]
    counts <- counts[Matrix::rowSums(counts > 0) >= min_cells, , drop = FALSE]
    obj <- CreateSeuratObject(counts, project = sample)
    # Preserve input QC measurements, including after the gene filter.
    meta <- metadata[cells, setdiff(colnames(metadata), c("qc_umap_1", "qc_umap_2")), drop = FALSE]
    obj <- AddMetaData(obj, meta)
    obj$keep_cell <- meta[[keep_col]]
    obj$doublet_filter_method <- doublet_method
    message(sample, " / ", doublet_method, ": ", ncol(obj), " cells; ", nrow(obj), " genes")
    obj
  }), samples)
  clean <- merge(objects[[1]], y = objects[-1], merge.data = FALSE)
  clean@misc$qc <- seurat_obj@misc$qc
  return(clean)
}
