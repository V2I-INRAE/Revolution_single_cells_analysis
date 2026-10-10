# Normalization diagnostics: per-sample HVG scatter, HVG overlap heatmap,
# PCA elbow and PC loadings.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggprism)
})

plot_pca_elbow <- function(seurat_obj, plot_dir, stem = "pca_elbow",
  ndims = 50, width = 8, height = 6) {
  n_computed <- ncol(Embeddings(seurat_obj, reduction = "pca"))
  p <- ElbowPlot(seurat_obj, reduction = "pca", ndims = min(ndims, n_computed)) +
    theme_prism() +
    labs(title = "PCA elbow plot")
  save_figure(p, plot_dir, stem, width, height)
  invisible(p)
}

plot_pca_loadings <- function(seurat_obj, plot_dir, stem = "pca_loadings",
  dims = 1:5, width = 12, height = 14) {
  n_computed <- ncol(Embeddings(seurat_obj, reduction = "pca"))
  dims <- dims[dims <= n_computed]
  if (length(dims) == 0L) {
    stop("None of the requested loading dimensions was computed.")
  }

  figures <- list()
  for (dim in dims) {
    p <- scCustomize::PC_Plotting(seurat_object = seurat_obj, dim_number = dim)
    p[[1]] <- p[[1]] + scale_fill_gradient2(
      low = "#3B4CC0", mid = "#F7F7F7", high = "#B40426", midpoint = 0,
      limits = c(-2.5, 2.5),
      breaks = c(-2.5, 0, 2.5), oob = scales::squish,
      name = "Scaled expression\n(clipped at ±2.5)"
    )
    # Match the heatmap's balanced genes rather than the wrapper's absolute top 30.
    p[[2]] <- VizDimLoadings(
      seurat_obj, dims = dim, nfeatures = 30, reduction = "pca", balanced = TRUE
    )
    p <- p + patchwork::plot_annotation(caption = paste(
      "Top 15 positive and 15 negative loading genes; up to 500 cells from the two PC-score extremes.",
      "Heatmap shows scaled expression, not loading values. Blue = negative; near-white = zero; red = positive.",
      sep = "\n"
    ))
    save_figure(p, plot_dir, paste0(stem, "_PC", dim), width, height)
    figures[[paste0("PC", dim)]] <- p
  }
  invisible(figures)
}

plot_variable_features <- function(seurat_obj, plot_dir, stem = "variable_features",
  width = 12, height = 9) {
  message("Plotting stored sample-specific variable features")
  assay_name <- DefaultAssay(seurat_obj)
  assay <- seurat_obj[[assay_name]]
  metadata <- seurat_obj[[]]
  figures <- list()
  for (layer in Layers(assay, search = "^counts($|\\.)")) {
    cells <- Cells(assay, layer = layer)
    sample <- unique(as.character(metadata[cells, "sample"]))
    stats <- HVFInfo(assay, method = "vst", layer = layer)
    hvgs <- VariableFeatures(assay, method = "vst", layer = layer, simplify = FALSE)[[1]]
    stopifnot(length(sample) == 1L, !anyNA(sample), !is.null(stats), length(hvgs) >= 20L)
    valid <- is.finite(stats$mean) & stats$mean > 0 &
      is.finite(stats$variance.standardized) & stats$variance.standardized > 0
    genes <- rownames(stats)[valid]
    message(sample, ": ", sum(!valid), " genes omitted from log axes (non-positive/missing statistics)")
    stopifnot(all(hvgs %in% genes))

    # Keep the normalization command history and only this sample's count layer.
    # HVFInfo then uses this layer, never the merged object's first sample.
    plot_obj <- subset(seurat_obj, cells = cells)
    plot_obj[[assay_name]] <- subset(plot_obj[[assay_name]], layers = layer, features = genes)
    VariableFeatures(plot_obj) <- hvgs
    set.seed(1234)
    p <- scCustomize::VariableFeaturePlot_scCustom(
      seurat_object = plot_obj, label = FALSE,
      y_axis_log = FALSE, pt.size = 0.4, colors_use = c("grey75", "#C44E52")
    )
    # The wrapper does not expose label layout controls; use its public labeling
    # function with zero nudges and more space for long pig Ensembl identifiers.
    p <- LabelPoints(p, points = head(hvgs, 20), repel = TRUE, xnudge = 0, ynudge = 0,
      size = 3, box.padding = 0.6, force = 3, max.overlaps = Inf,
      max.time = 3, seed = 1234) +
      scale_y_log10(expand = expansion(mult = c(0.05, 0.25))) +
      theme_classic(base_size = 12) + theme(legend.position = "bottom") +
      labs(title = sample, subtitle = paste(length(hvgs), "sample-specific HVGs; top 20 labeled"),
        x = "Mean expression (log10 scale)", y = "Standardized variance (log10 scale)",
        caption = paste(sum(!valid), "genes omitted from log axes; stored VST statistics, no reselection."))
    save_figure(p, plot_dir, paste0(stem, "_", sample), width, height)
    figures[[sample]] <- p
  }
  invisible(figures)
}

# -- To check HVG profile similitude
plot_hvg_overlap <- function(seurat_obj, plot_dir, stem = "hvg_overlap",
  width = 12, height = 8) {
  assay <- seurat_obj[[DefaultAssay(seurat_obj)]]
  consensus <- VariableFeatures(assay)
  metadata <- seurat_obj[[]]

  hvgs <- VariableFeatures(
    assay,
    method = "vst", layer = "^counts($|\\.)", simplify = FALSE
  )
  cells <- lapply(names(hvgs), function(layer) Cells(assay, layer = layer))

  names(hvgs) <- vapply(cells, function(ids) {
    sample <- unique(as.character(metadata[ids, "sample"]))
    if (length(sample) != 1L || anyNA(sample)) {
      stop("Each HVG layer must correspond to one sample.")
    }
    sample
  }, character(1))
  sample_order <- unique(as.character(metadata$sample))

  # -- Just to be sure we have one per sample and the consensus which are the
  # 3000 HVGs selected from the different samples
  if (anyDuplicated(names(hvgs)) || !setequal(names(hvgs), sample_order) ||
    any(lengths(hvgs) == 0L) || length(consensus) == 0L) {
    stop("Expected one nonempty HVG set per sample and a consensus selection.")
  }

  hvgs <- hvgs[sample_order]
  hvgs$Consensus <- consensus
  genes <- unique(unlist(hvgs, use.names = FALSE))
  membership <- t(vapply(
    hvgs, function(x) as.integer(genes %in% x),
    integer(length(genes))
  ))
  colnames(membership) <- genes

  heatmap <- pheatmap::pheatmap(
    membership,
    cluster_rows = FALSE, cluster_cols = TRUE,
    color = c("grey90", "grey20"), breaks = c(-0.5, 0.5, 1.5),
    legend_breaks = c(0, 1), legend_labels = c("Not selected", "Selected"),
    show_colnames = FALSE, border_color = NA,
    labels_row = paste0(names(hvgs), " (", lengths(hvgs), ")"),
    main = "LogNormalize HVG membership by sample",
    silent = TRUE
  )
  save_figure(heatmap$gtable, plot_dir, stem, width, height)
  invisible(membership)
}
