# Plotting functions for the QC pipeline

suppressPackageStartupMessages({
  library(ggplot2)
})

save_qc_plot <- function(plot, name, plot_dir, width, height) {
  png(file.path(plot_dir, paste0(name, ".png")),
    width = width, height = height, res = 300
  )
  print(plot)
  dev.off()
}

# --- Build a per-sample UMAP for a binary QC classification
build_qc_umap <- function(
  summary,
  class_col,
  class_levels,
  class_colours,
  point_sizes,
  point_alphas,
  ncol
) {
  summary$sample_id <- factor(
    summary$sample_id,
    levels = unique(summary$sample_id)
  )
  summary[[class_col]] <- factor(
    summary[[class_col]],
    levels = class_levels
  )
  summary <- summary[order(summary[[class_col]]), , drop = FALSE]

  ggplot(
    summary,
    aes(
      UMAP1,
      UMAP2,
      colour = .data[[class_col]],
      size = .data[[class_col]],
      alpha = .data[[class_col]]
    )
  ) +
    geom_point() +
    facet_wrap(~sample_id, ncol = ncol) +
    coord_equal() +
    scale_colour_manual(values = class_colours, drop = FALSE) +
    scale_size_manual(values = point_sizes, guide = "none") +
    scale_alpha_manual(values = point_alphas, guide = "none")
}

# --- Plot DoubletFinder calls on the per-sample diagnostic UMAPs
plot_doublet_comparison <- function(doublet_summary) {
  plot_dir <- file.path("results", "qc")
  class_colours <- c(Singlet = "#DDDDDB", Doublet = "#D83746")
  p_umap <- build_qc_umap(
    doublet_summary,
    class_col = "doublet_class",
    class_levels = c("Singlet", "Doublet"),
    class_colours = class_colours,
    point_sizes = c(Singlet = 0.1, Doublet = 0.1),
    point_alphas = c(Singlet = 0.4, Doublet = 0.4),
    ncol = 5
  ) +
    theme_classic(base_size = 12) +
    theme(strip.text = element_text(size = 8)) +
    labs(
      title = paste(
        "DoubletFinder calls per sample",
        "(UMAP coordinates are per-sample)"
      ),
      colour = "class"
    )

  save_qc_plot(
    p_umap,
    "doublet_comparison_umap",
    plot_dir,
    4000,
    3400
  )

  invisible(NULL)
}

# --- Compare nFeature_RNA between DoubletFinder classifications
plot_doublet_nfeature_violin <- function(doublet_summary) {
  plot_dir <- file.path("results", "qc")
  class_colours <- c(Singlet = "#DDDDDB", Doublet = "#D83746")
  doublet_summary$sample_id <- factor(
    doublet_summary$sample_id,
    levels = unique(doublet_summary$sample_id)
  )
  doublet_summary$doublet_class <- factor(
    doublet_summary$doublet_class,
    levels = c("Singlet", "Doublet")
  )

  # nFeature_RNA per sample, split by classification (doublets are
  # expected to have more detected genes than singlets)
  p_violin <- ggplot(
    doublet_summary,
    aes(
      sample_id,
      nFeature_RNA,
      fill = doublet_class
    )
  ) +
    geom_violin(linewidth = 0.3) +
    scale_fill_manual(values = class_colours) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "nFeature_RNA by DoubletFinder class",
      x = "sample", fill = "class"
    )

  save_qc_plot(
    p_violin,
    "doublet_comparison_nfeature_violin",
    plot_dir,
    3200,
    1800
  )

  invisible(NULL)
}

# --- Plot cells below the ribosomal percentage threshold on the same
# per-sample UMAP coordinates used for the DoubletFinder comparison.
plot_ribo_comparison <- function(qc_summary) {
  plot_dir <- file.path("results", "qc")
  threshold <- qc_params$metrics$min_ribo_percent
  low_label <- sprintf("< %.1f%%", threshold)
  retained_label <- sprintf(">= %.1f%%", threshold)
  qc_summary <- qc_summary[qc_summary$stage == "before", ]

  status_colours <- c("#DDDDDB", "#2C7FB8")
  names(status_colours) <- c(retained_label, low_label)

  point_sizes <- c(0.1, 0.15)
  names(point_sizes) <- c(retained_label, low_label)
  point_alphas <- c(0.3, 0.7)
  names(point_alphas) <- c(retained_label, low_label)

  p_umap <- build_qc_umap(
    qc_summary,
    class_col = "ribo_status",
    class_levels = c(retained_label, low_label),
    class_colours = status_colours,
    point_sizes = point_sizes,
    point_alphas = point_alphas,
    ncol = 6
  ) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    theme_void(base_size = 12) +
    theme(
      strip.text = element_text(size = 8),
      legend.position = "bottom",
      legend.direction = "horizontal"
    ) +
    labs(
      title = paste(
        "Ribosomal percentage labels per sample",
        "(UMAP coordinates are per-sample)"
      ),
      colour = "percent.ribo"
    )

  save_qc_plot(
    p_umap,
    "ribo_comparison_umap",
    plot_dir,
    4000,
    1800
  )

  invisible(NULL)
}

plot_qc_violin <- function(qc_summary, metric) {
  stage_colours <- c(before = "#3A5BA0", after = "#D35400")
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

  invisible(NULL)
}

plot_qc_scatter <- function(qc_summary, x, y, name) {
  stage_colours <- c(before = "#3A5BA0", after = "#D35400")
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

  invisible(NULL)
}


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
