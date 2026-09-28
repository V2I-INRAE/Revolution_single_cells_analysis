suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggprism)
})

plot_variable_features <- function(
  seurat_obj,
  output_dir = NULL,
  width = 9,
  height = 15,
  filename = "variable_features.png"
) {

  message("Plotting variable features")

  assay <- seurat_obj[[DefaultAssay(seurat_obj)]]
  consensus <- VariableFeatures(assay)
  metadata <- seurat_obj[[]]

  if (inherits(assay, "SCTAssay")) {
    tables <- lapply(levels(assay), function(model) {
      stats <- SCTResults(assay, slot = "feature.attributes", model = model)
      cells <- rownames(
        SCTResults(assay, slot = "cell.attributes", model = model)
      )
      sample <- unique(as.character(metadata[cells, "sample"]))
      stopifnot(length(sample) == 1, !anyNA(sample))
      data.frame(
        gene = rownames(stats), sample = sample,
        mean = stats$gmean, variance = stats$residual_variance
      )
    })
    axis_labels <- c("Geometric mean expression", "Residual variance")
  } else {
    tables <- lapply(Layers(assay, search = "^counts($|\\.)"), function(layer) {
      stats <- HVFInfo(assay, method = "vst", layer = layer)
      cells <- Cells(assay, layer = layer)
      sample <- unique(as.character(metadata[cells, "sample"]))
      stopifnot(length(sample) == 1, !anyNA(sample), !is.null(stats))
      data.frame(
        gene = rownames(stats), sample = sample,
        mean = stats$mean, variance = stats$variance.standardized
      )
    })
    axis_labels <- c("Mean expression", "Standardized variance")
  }

  plot_data <- do.call(rbind, tables)
  plot_data$selected <- plot_data$gene %in% consensus
  # A log-x scatter requires a positive mean and finite statistics.
  valid <- is.finite(plot_data$mean) & plot_data$mean > 0 &
    is.finite(plot_data$variance)
  if (any(!valid)) {
    message(
      "Rows omitted from plotting (non-finite statistics or non-positive mean):"
    )
    print(table(plot_data$sample[!valid]))
  }

  plot_data <- plot_data[valid, ]
  sample_order <- unique(as.character(metadata$sample))
  plot_data$sample <- factor(plot_data$sample, levels = sample_order)
  # Draw consensus features last so they remain visible over other genes.
  plot_data <- plot_data[order(plot_data$selected), ]

  p <- ggplot(
    plot_data,
    aes(mean, plot_data$variance, colour = plot_data$selected)
  ) +
    geom_point(size = 0.35, alpha = 0.6) +
    facet_wrap(~sample, ncol = 5) +
    scale_x_log10() +
    scale_y_sqrt() +
    scale_colour_manual(
      values = c("FALSE" = "grey75", "TRUE" = "#D83746"),
      breaks = c(FALSE, TRUE),
      labels = c("Other genes", "Consensus HVGs")
    ) +
    theme_prism() +
    theme(
      strip.text = element_text(size = 10),
      axis.text = element_text(size = 8),
      legend.position = "bottom"
    ) +
    labs(
      title = paste(DefaultAssay(seurat_obj), "variable features by sample"),
      subtitle = paste(
        length(consensus), "consensus HVGs highlighted across samples"
      ),
      x = axis_labels[1],
      y = paste(axis_labels[2], "(square-root scale)"), colour = NULL
    )

  # Save plot if output directory specified
  if (!is.null(output_dir)) {
    png_file <- file.path(output_dir, filename)
    ggsave(png_file, plot = p, width = width, height = height, dpi = 300)
    message("  Saved: ", png_file)
  }

  invisible(p)
}

# -- To check HVG profile similitude
plot_hvg_overlap <- function(
  seurat_obj,
  output_dir = NULL,
  width = 12,
  height = 8,
  filename = "hvg_overlap.png"
) {

  assay <- seurat_obj[[DefaultAssay(seurat_obj)]]
  consensus <- VariableFeatures(assay)
  metadata <- seurat_obj[[]]

  if (inherits(assay, "SCTAssay")) {
    # -- simplify = FALSE to return per-model HVGs and not the global hvgs.
    hvgs <- VariableFeatures(
      assay,
      layer = levels(assay), nfeatures = length(consensus),
      simplify = FALSE
    )
    cells <- lapply(names(hvgs), function(model) {
      rownames(SCTResults(assay, slot = "cell.attributes", model = model))
    })
    method <- "SCT"
  } else {
    hvgs <- VariableFeatures(
      assay,
      method = "vst", layer = "^counts($|\\.)", simplify = FALSE
    )
    cells <- lapply(names(hvgs), function(layer) Cells(assay, layer = layer))
    method <- "LogNormalize"
  }

  names(hvgs) <- vapply(cells, function(ids) {
    sample <- unique(as.character(metadata[ids, "sample"]))
    if (length(sample) != 1L || anyNA(sample)) {
      stop("Each HVG layer/model must correspond to one sample.")
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

  pheatmap::pheatmap(
    membership,
    cluster_rows = FALSE, cluster_cols = TRUE,
    color = c("grey90", "grey20"), breaks = c(-0.5, 0.5, 1.5),
    legend_breaks = c(0, 1), legend_labels = c("Not selected", "Selected"),
    show_colnames = FALSE, border_color = NA,
    labels_row = paste0(names(hvgs), " (", lengths(hvgs), ")"),
    main = paste(method, "HVG membership by sample"),
    filename = if (is.null(output_dir)) NA else file.path(output_dir, filename),
    width = width, height = height
  )
  invisible(membership)
}
