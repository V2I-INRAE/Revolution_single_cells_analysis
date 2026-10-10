# QC-only before/after figures from the labeled, unfiltered Seurat checkpoint.
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(scCustomize)
  library(patchwork)
})

save_qc_figure <- function(figure, name, plot_dir, width, height) {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(plot_dir, paste0(name, ".png")), figure,
    width = width, height = height, units = "in", dpi = 300, bg = "white")
}

plot_qc_comparison <- function(seurat_obj, plot_dir = "results/qc") {
  metadata <- seurat_obj[[]]
  samples <- seurat_obj@misc$qc$sample_order
  params <- seurat_obj@misc$qc$params$filtering
  stopifnot(
    !anyNA(metadata$keep_qc),
    identical(metadata$keep_qc, !metadata$qc_outlier),
    setequal(samples, metadata$sample)
  )
  seurat_obj$sample <- factor(metadata$sample, levels = samples)
  retained_cells <- rownames(metadata)[metadata$keep_qc]
  retained <- subset(seurat_obj, cells = retained_cells)
  # Plot input metrics, never values recomputed by a subsetting operation.
  metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "log10GenesPerUMI", "percent.ribo")
  retained <- AddMetaData(retained, metadata[retained_cells, metrics, drop = FALSE])

  colours <- setNames(c(
    "#3A5BA0", "#D4753C", "#5A8F5A", "#C44E52", "#7B5EA7",
    "#E8A838", "#46878F", "#B07AA1", "#2E86C1", "#8C6D31",
    "#4E9A9A", "#D98880", "#6B8E23", "#9B59B6", "#1ABC9C",
    "#86714D", "#8EC9EB", "#6E2F84", "#F5A623", "#7BC657",
    "#708090", "#8B4513", "#C5A9D8"
  )[seq_along(samples)], samples)

  specs <- list(
    genes = list(
      fun = QC_Plots_Genes,
      metric = "nFeature_RNA",
      label = "Genes detected per cell",
      low = params$min_features,
      high = params$max_features
    ),
    umi = list(
      fun = QC_Plots_UMIs,
      metric = "nCount_RNA",
      label = "UMIs per cell"
    ),
    mitochondrial = list(
      fun = QC_Plots_Mito,
      metric = "percent.mt",
      label = "Mitochondrial counts (%)",
      high = params$max_mito_percent,
      extra = list(mito_name = "percent.mt")
    ),
    complexity = list(
      fun = QC_Plots_Complexity,
      metric = "log10GenesPerUMI",
      label = "Complexity: log10(genes) / log10(UMIs)",
      low = params$min_log10_genes_per_umi
    ),
    ribosomal = list(
      fun = QC_Plots_Feature,
      metric = "percent.ribo",
      label = "Ribosomal counts (%)",
      extra = list(feature = "percent.ribo")
    )
  )

  figures <- list()
  for (name in names(specs)) {
    spec <- specs[[name]]
    limits <- range(c(metadata[[spec$metric]], spec$low, spec$high))
    if (name != "complexity") limits[1] <- 0
    panels <- lapply(list(seurat_obj, retained), function(obj) {
      set.seed(1234)
      p <- do.call(spec$fun, c(list(
        seurat_object = obj, group.by = "sample",
        plot_title = NULL, x_axis_label = NULL, y_axis_label = spec$label,
        low_cutoff = spec$low, high_cutoff = spec$high,
        colors_use = colours, pt.size = 0.05, raster = TRUE,
        raster.dpi = 300, add.noise = FALSE, layer = "counts"
      ), spec$extra))
      sizes <- table(factor(obj$sample, levels = samples))
      p + scale_x_discrete(
        limits = samples,
        drop = FALSE,
        labels = paste0(samples, "\nn=", format(sizes, trim = TRUE))
      ) +
        # Replace VlnPlot's data-only limits so cap lines outside the observed
        # range are not discarded before the plot is drawn.
        scale_y_continuous(limits = limits, expand = expansion(mult = 0.04)) +
        theme_classic(base_size = 12) +
        theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          plot.title = element_text(face = "bold", size = 13)
        )
    })
    panels[[1]] <- panels[[1]] + labs(title = "Before: all vendor-called cells")
    panels[[2]] <- panels[[2]] + labs(title = "After: QC-passing cells (no doublet removal)"
    )
    caption <- ("Original per-cell QC measurements; UMI counts do not directly control filtering.")
    if (name %in% c("genes", "mitochondrial")) {
      caption <- paste0(caption,
        "\nDashed upper line: fixed threshold shared by all samples. Cells strictly above are excluded.")
    }
    if (name == "ribosomal") {
      caption <- paste0(caption,
        "\nRibosomal percentage is diagnostic only; no ribosomal filtering is applied.")
    }
    figure <- wrap_plots(panels, ncol = 1) + plot_annotation(caption = caption)
    save_qc_figure(figure, paste0("qc_before_after_", name), plot_dir, 16, 10)
    figures[[name]] <- figure
  }
  invisible(figures)
}

# Each sample has its own diagnostic embedding. Never pool these coordinates
# into a joint UMAP or recompute them when changing the highlighted label.
plot_qc_umaps <- function(seurat_obj, plot_dir = "results/qc") {
  metadata <- seurat_obj[[]]
  required <- c(
    "sample",
    "predicted_doublets",
    "qc_outlier",
    "qc_umap_1",
    "qc_umap_2"
  )

  stopifnot(all(required %in% colnames(metadata)), !anyNA(metadata[, required]))

  samples <- seurat_obj@misc$qc$sample_order

  stopifnot(
    setequal(samples, metadata$sample),
    all(is.finite(as.matrix(metadata[, c("qc_umap_1", "qc_umap_2")])))
  )

  specs <- list(
    scrublet = list(
      title = "Scrublet", negative = "Singlet",
      positive = "Doublet", flag = metadata$predicted_doublets, colour = "#C44E52"
    ),
    outliers = list(
      title = "QC outliers", negative = "QC pass",
      positive = "QC outlier", flag = metadata$qc_outlier, colour = "#D4753C"
    )
  )

  panels <- setNames(lapply(specs, function(x) list()), names(specs))
  for (sample in samples) {
    cells <- rownames(metadata)[metadata$sample == sample]
    obj <- subset(seurat_obj, cells = cells)
    coordinates <- as.matrix(metadata[colnames(obj), c("qc_umap_1", "qc_umap_2")])
    colnames(coordinates) <- c("QCUMAP_1", "QCUMAP_2")
    obj[["qc_umap"]] <- CreateDimReducObject(
      embeddings = coordinates, key = "QCUMAP_", assay = "RNA"
    )
    for (name in names(specs)) {
      spec <- specs[[name]]
      flag <- spec$flag[match(colnames(obj), rownames(metadata))]
      subtitle <- sprintf("%s: %s / %s cells (%.1f%%)", spec$positive,
        sum(flag), length(flag), 100 * mean(flag))
      obj$qc_plot_class <- factor(ifelse(flag, spec$positive, spec$negative),
        levels = c(spec$negative, spec$positive))
      colours <- setNames(c("#A0A0A0", spec$colour), levels(obj$qc_plot_class))
      panels[[name]][[sample]] <- scplotter::CellDimPlot(
        obj, reduction = "qc_umap", group_by = "qc_plot_class",
        highlight = sprintf('qc_plot_class == "%s"', spec$positive),
        pt_size = 0.25, pt_alpha = 0.5, highlight_size = 0.25,
        highlight_alpha = 1, highlight_color = "black", highlight_stroke = 0.05,
        palcolor = colours, order = "high-top", show_stat = FALSE, label = FALSE,
        raster = FALSE, theme = ggplot2::theme_classic,
        theme_args = list(base_size = 12), legend.position = "bottom",
        legend.direction = "horizontal",
        title = sample, xlab = "UMAP 1", ylab = "UMAP 2", seed = 1234
      ) + scale_colour_manual(values = colours, limits = names(colours),
          drop = FALSE, name = NULL) +
        guides(colour = guide_legend(nrow = 1, override.aes = list(alpha = 1, size = 2))) +
        labs(subtitle = subtitle) +
        theme(plot.title = element_text(face = "bold"))
    }
  }
  n_columns <- min(4L, length(samples))
  figures <- list()
  for (name in names(specs)) {
    figure <- wrap_plots(panels[[name]], ncol = n_columns, guides = "collect") +
      plot_annotation(title = specs[[name]]$title,
        caption = "All vendor-called cells. Each panel is an independent sample UMAP; coordinates are identical across these figures.") &
      theme(legend.position = "bottom")
    save_qc_figure(figure, paste0("qc_umap_", name), plot_dir,
      width = 5 * n_columns, height = 4.2 * ceiling(length(samples) / n_columns) + 1)
    figures[[name]] <- figure
  }
  invisible(figures)
}
