suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggprism)
})

plot_pca_elbow <- function(
  seurat_obj,
  output_dir = NULL,
  ndims = 50,
  width = 8,
  height = 6,
  filename = "pca_elbow.png"
) {
  n_computed <- ncol(Embeddings(seurat_obj, reduction = "pca"))
  p <- ElbowPlot(
    seurat_obj,
    reduction = "pca",
    ndims = min(ndims, n_computed)
  ) +
    theme_prism() +
    labs(title = "PCA elbow plot")

  if (!is.null(output_dir)) {
    png_file <- file.path(output_dir, filename)
    ggsave(png_file, plot = p, width = width, height = height, dpi = 300)
    message("  Saved: ", png_file)
  }

  invisible(p)
}

plot_pca_loadings <- function(
  seurat_obj,
  output_dir = NULL,
  dims = 1:5,
  n_features = 30,
  width = 15,
  height = 12,
  filename = "pca_loadings.png"
) {
  n_computed <- ncol(Embeddings(seurat_obj, reduction = "pca"))
  dims <- dims[dims <= n_computed]
  if (length(dims) == 0L) {
    stop("None of the requested loading dimensions was computed.")
  }

  p <- VizDimLoadings(
    seurat_obj,
    dims = dims,
    nfeatures = n_features,
    reduction = "pca",
    ncol = min(3L, length(dims)),
    balanced = TRUE
  ) &
    theme_prism()

  if (!is.null(output_dir)) {
    png_file <- file.path(output_dir, filename)
    ggsave(png_file, plot = p, width = width, height = height, dpi = 300)
    message("  Saved: ", png_file)
  }

  invisible(p)
}

plot_grouping <- function(
  metadata,
  group_by = c("sample", "pressure", "time", "pressure_time")
) {
  group_by <- match.arg(group_by)
  required <- switch(
    group_by,
    sample = "sample",
    pressure = "pressure",
    time = "time_point",
    pressure_time = c("pressure", "time_point")
  )
  missing_columns <- setdiff(required, colnames(metadata))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing plotting metadata: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  pressure <- as.character(metadata$pressure)
  time <- as.character(metadata$time_point)
  if (group_by != "sample") {
    unknown_pressure <- setdiff(
      unique(pressure),
      c("None", "C", "positive", "negative")
    )
    if (length(unknown_pressure) > 0L || anyNA(pressure)) {
      stop("Unexpected or missing pressure values in PCA metadata.")
    }
    pressure[pressure == "None"] <- "C"
  }
  if (group_by %in% c("time", "pressure_time")) {
    unknown_time <- setdiff(unique(time), c("T0H", "T4H", "T8H", "T10H"))
    if (length(unknown_time) > 0L || anyNA(time)) {
      stop("Unexpected or missing time values in PCA metadata.")
    }
    time[time == "T8H"] <- "T10H"
  }

  if (group_by == "sample") {
    values <- as.character(metadata$sample)
    if (anyNA(values) || any(values == "")) {
      stop("Missing sample values in PCA metadata.")
    }
    values <- factor(values, levels = unique(values))
    legend_title <- "Sample"
  } else if (group_by == "pressure") {
    values <- factor(pressure, levels = c("C", "positive", "negative"))
    legend_title <- "Pressure"
  } else if (group_by == "time") {
    values <- factor(time, levels = c("T0H", "T4H", "T10H"))
    legend_title <- "Time"
  } else {
    if (any(time != "T0H" & pressure == "C")) {
      stop("Control pressure is only expected at T0H.")
    }
    values <- ifelse(time == "T0H", "C", paste(pressure, time, sep = "_"))
    values <- factor(
      values,
      levels = c(
        "C", "positive_T4H", "negative_T4H",
        "positive_T10H", "negative_T10H"
      )
    )
    legend_title <- "Pressure and time"
  }
  if (anyNA(values)) {
    stop("PCA grouping produced missing labels.")
  }

  list(values = values, title = legend_title)
}

plot_pca_grouping <- function(
  seurat_obj,
  group_by = c("sample", "pressure", "time", "pressure_time"),
  output_dir = NULL,
  seed = 1234,
  width = 15,
  height = 5,
  filename = NULL
) {
  group_by <- match.arg(group_by)
  n_computed <- ncol(Embeddings(seurat_obj, reduction = "pca"))
  if (n_computed < 6L) {
    stop("The requested PC1-6 panel requires at least six computed PCs.")
  }
  grouping <- plot_grouping(seurat_obj[[]], group_by)
  plot_obj <- seurat_obj
  plot_obj$pca_group <- grouping$values
  panels <- lapply(list(1:2, 3:4, 5:6), function(dims) {
    DimPlot(
      plot_obj,
      reduction = "pca",
      group.by = "pca_group",
      dims = dims,
      shuffle = TRUE,
      seed = seed
    ) +
      labs(colour = grouping$title)
  })
  p <- patchwork::wrap_plots(panels, ncol = 3) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  if (!is.null(output_dir)) {
    if (is.null(filename)) {
      filename <- paste0("pca_", group_by, ".png")
    }
    png_file <- file.path(output_dir, filename)
    ggsave(png_file, plot = p, width = width, height = height, dpi = 300)
    message("  Saved: ", png_file)
  }

  invisible(p)
}

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
