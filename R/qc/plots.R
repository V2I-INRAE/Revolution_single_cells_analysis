# Plotting functions for the QC pipeline

suppressPackageStartupMessages({
  library(ggplot2)
})

save_doublet_plot <- function(plot, name, plot_dir, width, height) {
  png(file.path(plot_dir, paste0(name, ".png")),
    width = width, height = height, res = 300
  )
  print(plot)
  dev.off()
}

# --- Plot to check doublet cells on umaps and as violin plots in each sample
plot_doublet_comparison <- function(doublet_summary) {
  plot_dir <- file.path("results", "qc")
  doublet_summary$sample_id <- factor(
    doublet_summary$sample_id,
    levels = unique(doublet_summary$sample_id)
  )
  class_colours <- c(Singlet = "#DDDDDB", Doublet = "#D83746")
  # factor levels: Singlet drawn first, Doublet last (on top of the
  # singlet cloud, where the calls are actually visible)
  doublet_summary$doublet_class <- factor(
    doublet_summary$doublet_class,
    levels = c("Singlet", "Doublet")
  )

  # per-sample UMAP facets coloured by classification
  p_umap <- ggplot(
    doublet_summary,
    aes(
      UMAP1,
      UMAP2,
      colour = doublet_class
    )
  ) +
    geom_point(size = 0.1, alpha = 0.4) +
    facet_wrap(~sample_id, ncol = 5) +
    coord_equal() +
    scale_colour_manual(values = class_colours) +
    theme_classic(base_size = 12) +
    theme(strip.text = element_text(size = 8)) +
    labs(
      title = paste(
        "DoubletFinder calls per sample",
        "(UMAP coordinates are per-sample)"
      ),
      colour = "class"
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

  save_doublet_plot(
    p_umap,
    "doublet_comparison_umap",
    plot_dir,
    4000,
    3400
  )
  save_doublet_plot(
    p_violin,
    "doublet_comparison_nfeature_violin",
    plot_dir,
    3200,
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
