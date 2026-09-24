# Import and clean filtered BD Rhapsody MEX for REVO26-C into a Seurat object

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

read_filtered_matrix <- function(folder_raw, sample_id) {
  # --- get path to the zipped filtered matrix
  mex_zip <-file.path(
    folder_raw,
    sample_id,
    paste0(sample_id, "_RSEC_MolsPerCell_MEX.zip")
  )
  stopifnot("MEX zip not found" = file.exists(mex_zip))

  # --- unzip the matrix to barcodes, features and matrix count
  mex_dir <- file.path(tempdir(), paste0(sample_id, "_filtered_MEX"))
  utils::unzip(mex_zip, exdir = mex_dir)
  stopifnot(all(file.exists(file.path(
    mex_dir,
    c("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz")
  ))))

  return(mex_dir)
}

build_seurat_obj <- function(unzipped_mex_dir, sample_id) {
  counts <- Seurat::Read10X(data.dir = unzipped_mex_dir) # genes x cells, RSEC counts
  colnames(counts) <- paste(sample_id, colnames(counts), sep = "_")

  # --- Pig mitochondrial gene names do not have the MT-prefix
  mito_genes <- c(
    "ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB",
    "ND1", "ND2", "ND3", "ND4", "ND5", "ND6"
  )
  idx <- match(mito_genes, rownames(counts))
  stopifnot(!anyNA(idx))
  rownames(counts)[idx] <- paste0("MT-", mito_genes)
  stopifnot(anyDuplicated(rownames(counts)) == 0)

  seurat_obj <- Seurat::CreateSeuratObject(counts = counts, project = sample_id)
  return(seurat_obj)
}

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

# One violin figure per metric: 23 panels (5 columns), before and after
# filtering side by side in each panel. Saved to results/qc/.
plot_qc_violin <- function(qc_summary, metric) {
  stage_colours <- c(before = "grey70", after = "#2C7FB8")
  p <- ggplot(qc_summary, aes(stage, .data[[metric]], fill = stage)) +
    geom_violin(linewidth = 0.3) +
    facet_wrap(~sample_id, ncol = 5) +
    scale_fill_manual(values = stage_colours) +
    theme_classic(base_size = 12) +
    theme(
      strip.text = element_text(size = 8),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    labs(
      title = paste("QC violin -", metric),
      x = NULL, y = metric, fill = "stage"
    )
  png(file.path("results", "qc", paste0("qc_violin_", metric, ".png")),
    width = 4000, height = 3400, res = 300
  )
  print(p)
  dev.off()
}

# One scatter figure: 23 panels (5 columns), "before" points drawn under
# the "after" ones so the cells removed by filtering stay visible.
plot_qc_scatter <- function(qc_summary, x, y, name) {
  stage_colours <- c(before = "grey70", after = "#2C7FB8")
  p <- ggplot(qc_summary, aes(.data[[x]], .data[[y]], colour = stage)) +
    geom_point(size = 0.15, alpha = 0.3) +
    facet_wrap(~sample_id, ncol = 5) +
    scale_colour_manual(values = stage_colours) +
    theme_classic(base_size = 12) +
    theme(strip.text = element_text(size = 8)) +
    labs(
      title = paste("QC scatter -", y, "vs", x),
      x = x, y = y, colour = "stage"
    )
  png(file.path("results", "qc", paste0("qc_scatter_", name, ".png")),
    width = 4000, height = 3400, res = 300
  )
  print(p)
  dev.off()
}

# Cross-sample QC figures, built from the per-cell summaries collected in
# main.R (before and after filtering). Called once, outside the loop.
plot_qc_comparison <- function(qc_summary) {
  qc_summary$stage <- factor(qc_summary$stage, levels = c("before", "after"))
  metrics <- c(
    "nFeature_RNA", "nCount_RNA", "percent.mt",
    "percent.ribo", "log10GenesPerUMI"
  )
  for (metric in metrics) plot_qc_violin(qc_summary, metric)
  plot_qc_scatter(qc_summary, "nCount_RNA", "nFeature_RNA", "features")
  plot_qc_scatter(qc_summary, "nFeature_RNA", "percent.mt", "mt")
  invisible(NULL)
}

save_qc_plot <- function(plot, name, width, height, before_after, sample_id) {
  plot_dir <- file.path("results", "qc", before_after, sample_id)

  png(file.path(plot_dir, paste0(sample_id, "_", name, ".png")),
    width = width, height = height, res = 300
  )
  print(plot)
  dev.off()
}

plot_qc <- function(seurat_obj, before_after, sample_id) {
  metrics <- c(
    "nFeature_RNA", "nCount_RNA", "percent.mt",
    "percent.ribo", "log10GenesPerUMI"
  )

  violin_plots <- VlnPlot(
    seurat_obj,
    features = metrics,
    ncol = 3, pt.size = 0
  ) + plot_layout(ncol = 3)

  scatter_count_feature <- FeatureScatter(
    seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"
  )

  scatter_feature_mt <- FeatureScatter(
    seurat_obj, feature1 = "nFeature_RNA", feature2 = "percent.mt"
  )

  md_long <- do.call(rbind, lapply(metrics, function(m) {
    data.frame(metric = m, values = seurat_obj[[m]][, 1])
  }))

  qc_histograms <- ggplot(md_long, aes(md_long$values)) +
    geom_histogram(
      bins = 60, fill = "grey25", colour = "white",
      linewidth = 0.1
    ) +
    facet_wrap(~metric, scales = "free") +
    theme_classic(base_size = 12) +
    labs(
      x = NULL, y = "Cells",
      title = paste(sample_id, "- QC metric distributions")
    )

  save_qc_plot(violin_plots, "qc_violin", 3600, 2400, before_after, sample_id)
  save_qc_plot(scatter_count_feature, "scatter_count_vs_feature", 2400, 2000, before_after, sample_id)
  save_qc_plot(scatter_feature_mt, "scatter_feature_vs_mt", 2400, 2000, before_after, sample_id)
  save_qc_plot(qc_histograms, "qc_histograms", 3600, 2000, before_after, sample_id)
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
  thresholds <- list(
    nCount_RNA = mad_bounds(seurat_obj$nCount_RNA),
    nFeature_RNA = mad_bounds(seurat_obj$nFeature_RNA, lower_floor = 200),
    percent.mt = mad_bounds(seurat_obj$percent.mt, upper_cap = 20),
    log10GenesPerUMI = mad_bounds(seurat_obj$log10GenesPerUMI) # only lower bound used
  )
  cat("\nMAD thresholds:\n")
  print(round(t(sapply(thresholds, identity)), 2))

  seurat_obj$flag_low_count <- seurat_obj$nCount_RNA < thresholds$nCount_RNA["lower"]
  seurat_obj$flag_high_count <- seurat_obj$nCount_RNA > thresholds$nCount_RNA["upper"]
  seurat_obj$flag_low_genes <- seurat_obj$nFeature_RNA < thresholds$nFeature_RNA["lower"]
  seurat_obj$flag_high_genes <- seurat_obj$nFeature_RNA > thresholds$nFeature_RNA["upper"]
  seurat_obj$flag_high_mt <- seurat_obj$percent.mt > thresholds$percent.mt["upper"]
  seurat_obj$flag_low_complexity <- seurat_obj$log10GenesPerUMI < thresholds$log10GenesPerUMI["lower"]

  flag_cols <- grep("^flag_", colnames(seurat_obj[[]]), value = TRUE)
  seurat_obj$qc_outlier <- rowSums(seurat_obj[[]][, flag_cols]) > 0

  cat("\nFlagged cells:\n")
  for (col in c(flag_cols, "qc_outlier")) {
    v <- seurat_obj[[col]][[1]]
    cat(sprintf(
      "  %-20s %6d (%5.1f%%)\n",
      col, sum(v), 100 * mean(v)
    ))
  }
  cat(sprintf(
    "Cells kept if outliers removed: %d (%.1f%%)\n",
    sum(!seurat_obj$qc_outlier), 100 * mean(!seurat_obj$qc_outlier)
  ))

  n_cells_imported <- ncol(seurat_obj)
  clean_seurat_obj <- subset(seurat_obj, subset = qc_outlier == FALSE)
  cat(sprintf(
    "\nCells after filtering: %d (%.1f%% of imported)\n",
    ncol(clean_seurat_obj), 100 * ncol(clean_seurat_obj) / n_cells_imported
  ))

  return(clean_seurat_obj)
}
