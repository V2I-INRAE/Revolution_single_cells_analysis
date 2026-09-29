# Doublet-specific plots. Generic saving and binary UMAP helpers come from qc.

plot_scrublet_score_histograms <- function(scores, thresholds) {
  scores$sample_id <- factor(
    scores$sample_id,
    levels = unique(scores$sample_id)
  )
  thresholds$sample_id <- factor(
    thresholds$sample_id,
    levels = levels(scores$sample_id)
  )
  threshold_data <- rbind(
    data.frame(
      sample_id = thresholds$sample_id,
      threshold = thresholds$scrublet_current_threshold,
      threshold_type = "Current: score > 0.15"
    ),
    data.frame(
      sample_id = thresholds$sample_id,
      threshold = thresholds$scrublet_automatic_threshold,
      threshold_type = "Automatic"
    )
  )
  threshold_data <- threshold_data[is.finite(threshold_data$threshold), ]

  ggplot(
    scores,
    aes(score, after_stat(density), fill = population)
  ) +
    geom_histogram(
      binwidth = 0.01,
      boundary = 0,
      position = "identity",
      alpha = 0.45
    ) +
    geom_vline(
      data = threshold_data,
      aes(xintercept = threshold, colour = threshold_type),
      linewidth = 0.45
    ) +
    facet_wrap(~sample_id, ncol = 5, scales = "free_y") +
    scale_fill_manual(values = c(
      "Observed cells" = "#2C7FB8",
      "Simulated doublets" = "#D95F0E"
    )) +
    scale_colour_manual(values = c(
      "Current: score > 0.15" = "#7A0177",
      "Automatic" = "#238B45"
    )) +
    theme_classic(base_size = 11) +
    theme(strip.text = element_text(size = 8), legend.position = "bottom") +
    labs(
      title = "Scrublet scores for observed cells and simulated doublets",
      subtitle = "Density is normalized separately for each population",
      x = "Scrublet score",
      y = "Density",
      fill = NULL,
      colour = "Threshold"
    )
}

plot_pk_sweeps <- function(bcmvn) {
  bcmvn$sample_id <- factor(bcmvn$sample_id, levels = unique(bcmvn$sample_id))
  ggplot(bcmvn, aes(pK_numeric, BCmetric)) +
    geom_line(linewidth = 0.35, colour = "#4D4D4D") +
    geom_point(size = 0.7, colour = "#4D4D4D") +
    geom_point(
      data = bcmvn[bcmvn$selected, ],
      aes(colour = selected_at_boundary),
      size = 2
    ) +
    facet_wrap(~sample_id, ncol = 5, scales = "free_y") +
    scale_colour_manual(
      values = c(`FALSE` = "#2C7FB8", `TRUE` = "#D7301F"),
      labels = c(`FALSE` = "Selected", `TRUE` = "Selected at boundary")
    ) +
    theme_classic(base_size = 11) +
    theme(strip.text = element_text(size = 8), legend.position = "bottom") +
    labs(
      title = "DoubletFinder pK selection curves",
      subtitle = "The maximum BCmvn marks the selected pK; it is not accuracy",
      x = "pK",
      y = "BCmvn",
      colour = NULL
    )
}

plot_continuous_score_umap <- function(
  cell_data,
  score_col,
  title,
  legend_title
) {
  plot_data <- data.frame(
    sample_id = cell_data$sample_id,
    UMAP1 = cell_data$UMAP1,
    UMAP2 = cell_data$UMAP2,
    score = cell_data[[score_col]]
  )
  plot_data$sample_id <- factor(
    plot_data$sample_id,
    levels = unique(plot_data$sample_id)
  )

  ggplot(plot_data, aes(UMAP1, UMAP2, colour = score)) +
    geom_point(size = 0.12, alpha = 0.65) +
    facet_wrap(~sample_id, ncol = 5) +
    coord_equal() +
    scale_colour_viridis_c(option = "magma") +
    theme_void(base_size = 10) +
    theme(
      strip.text = element_text(size = 8),
      legend.position = "bottom",
      plot.margin = margin(10, 10, 25, 10)
    ) +
    labs(
      title = title,
      subtitle = "Coordinates are shared between methods within each sample",
      colour = legend_title
    )
}

plot_monitor_binary_umap <- function(cell_data, class_col, title) {
  plot_data <- data.frame(
    sample_id = cell_data$sample_id,
    UMAP1 = cell_data$UMAP1,
    UMAP2 = cell_data$UMAP2,
    doublet_class = cell_data[[class_col]]
  )
  plot_data <- plot_data[!is.na(plot_data$doublet_class), ]
  build_qc_umap(
    plot_data,
    class_col = "doublet_class",
    class_levels = c("Singlet", "Doublet"),
    class_colours = c(Singlet = "#DDDDDB", Doublet = "#D83746"),
    point_sizes = c(Singlet = 0.1, Doublet = 0.15),
    point_alphas = c(Singlet = 0.35, Doublet = 0.75),
    ncol = 5
  ) +
    theme_void(base_size = 11) +
    theme(strip.text = element_text(size = 8), legend.position = "bottom") +
    labs(
      title = title,
      subtitle = "UMAP coordinates are computed independently per sample",
      colour = "Class"
    )
}

plot_nfeature_summary <- function(summary) {
  summary$sample_id <- factor(
    summary$sample_id,
    levels = unique(summary$sample_id)
  )
  summary$class <- factor(summary$class, levels = c("Singlet", "Doublet"))
  summary$panel <- ifelse(
    summary$method == "DoubletFinder",
    "DoubletFinder — selected pK",
    paste("Scrublet —", summary$rule)
  )
  ggplot(
    summary,
    aes(sample_id, median, colour = class, group = class)
  ) +
    geom_linerange(
      aes(ymin = q1, ymax = q3),
      position = position_dodge(width = 0.45),
      linewidth = 0.7
    ) +
    geom_point(position = position_dodge(width = 0.45), size = 1.8) +
    geom_text(
      aes(y = q3, label = paste0("n=", n)),
      position = position_dodge(width = 0.45),
      angle = 90,
      hjust = -0.1,
      size = 2.2,
      show.legend = FALSE
    ) +
    facet_wrap(~panel, ncol = 1) +
    scale_colour_manual(values = c(Singlet = "#636363", Doublet = "#D83746")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 9),
      legend.position = "bottom",
      plot.margin = margin(15, 15, 15, 30)
    ) +
    labs(
      title = "Detected genes by doublet classification",
      subtitle = "Points are medians; ranges are Q1-Q3; labels are group sizes",
      x = "Sample",
      y = "nFeature_RNA",
      colour = "Class"
    )
}
