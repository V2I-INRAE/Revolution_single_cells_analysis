# The driver sources R/utils/save_plots.R before calling this module.
plot_pseudobulk_diagnostics <- function(pb, output_dir, label, seed = 1234) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  sample_info <- pb$sample_info
  sample_info$library_before_filter <- as.numeric(colSums(pb$counts)[sample_info$sample])
  mds <- limma::plotMDS(pb$dge, top = 500, gene.selection = "pairwise",
    method = "logFC", prior.count = 2, dim.plot = c(1, 2), plot = FALSE)
  mds_data <- data.frame(sample_info, MDS1 = mds$x, MDS2 = mds$y)
  csv <- file.path(output_dir, c(sample_qc = "sample_qc.csv",
    mds_coordinates = "mds_coordinates.csv"))
  names(csv) <- c("sample_qc", "mds_coordinates")
  write.csv(sample_info, csv[["sample_qc"]], row.names = FALSE)
  write.csv(mds_data, csv[["mds_coordinates"]], row.names = FALSE)
  rds <- file.path(output_dir, "mds.rds")
  saveRDS(list(mds = mds, sample_info = sample_info), rds)

  condition_colors <- setNames(viz_palette[seq_along(levels(sample_info$condition))],
    levels(sample_info$condition))
  pig_colors <- setNames(viz_palette[seq_along(levels(sample_info$pig))], levels(sample_info$pig))
  axis_labels <- sprintf("Leading logFC dimension %d\n(%.1f%%)",
    1:2, 100 * mds$var.explained[1:2])
  mds_plots <- lapply(c("condition", "pig"), function(group) {
    colors <- if (group == "condition") condition_colors else pig_colors
    plotthis::ScatterPlot(mds_data, x = "MDS1", y = "MDS2", color_by = group,
      color_name = if (group == "condition") "Condition" else "Pig", palcolor = colors,
      size_by = 3, alpha = 1, border_size = 0.3, legend.position = "bottom",
      aspect.ratio = NULL,
      theme = "theme_classic", theme_args = list(base_size = 12), combine = FALSE,
      title = if (group == "condition") "Samples by condition" else "Samples by pig",
      xlab = axis_labels[1], ylab = axis_labels[2], seed = seed)[[1]] +
      ggrepel::geom_text_repel(data = mds_data,
        ggplot2::aes(MDS1, MDS2, label = sample), inherit.aes = FALSE,
        size = 3, seed = seed, max.overlaps = Inf, min.segment.length = 0,
        box.padding = 0.6, force = 10, force_pull = 0.1, max.iter = 100000, max.time = Inf) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.15)) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.15)) +
      ggplot2::guides(color = ggplot2::guide_legend(nrow = 1),
        fill = ggplot2::guide_legend(nrow = 1)) +
      ggplot2::theme(legend.direction = "horizontal") +
      ggplot2::coord_equal()
  })
  mds_figure <- patchwork::wrap_plots(mds_plots, ncol = 2) + patchwork::plot_annotation(
    title = paste0(label, ": pseudobulk sample MDS"),
    subtitle = sprintf("%d samples from %d pigs. Identical coordinates in both panels; animal effects retained.",
      nrow(sample_info), length(unique(sample_info$pig))),
    caption = "TMM-normalized expression; leading logFC distances using the top 500 genes per sample pair.",
    tag_levels = "A")

  qc_data <- sample_info
  qc_data$sample <- factor(qc_data$sample,
    levels = with(sample_info, sample[order(pig, condition)]))
  qc_data$library_millions <- qc_data$library_size / 1e6
  qc_data$effective_library_millions <- qc_data$effective_library_size / 1e6
  metrics <- c(n_cells = paste(label, "cells"), library_millions = "Library size (million counts)",
    tmm_factor = "TMM normalization factor", effective_library_millions = "Effective library size (millions)")
  metric_plots <- lapply(names(metrics), function(metric) {
    p <- plotthis::BarPlot(qc_data, x = "sample", y = metric,
      fill_by = "condition", fill_name = "Condition", palcolor = condition_colors,
      x_text_angle = 60, aspect.ratio = 0.55, y_min = 0,
      y_max = 1.12 * max(qc_data[[metric]]), theme = "theme_classic",
      theme_args = list(base_size = 12), combine = FALSE,
      xlab = "Sample (grouped by pig)", ylab = metrics[[metric]],
      title = metrics[[metric]], seed = seed)[[1]]
    if (metric == "tmm_factor") {
      p <- p + ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey40")
    }
    p
  })
  qc_figure <- (patchwork::wrap_plots(metric_plots, ncol = 2, guides = "collect") +
    patchwork::plot_annotation(
      title = paste0(label, ": sample-level count and normalization diagnostics"),
      subtitle = paste("All", nrow(sample_info), "samples retained. Libraries use the",
        format(nrow(pb$dge), big.mark = ","), "genes retained by filterByExpr."),
      caption = "Effective library size = filtered library size x TMM factor. Dashed line marks TMM factor 1, not an exclusion cutoff.",
      tag_levels = "A")) & ggplot2::theme(legend.position = "bottom")
  save_figure(mds_figure, output_dir, "sample_mds", width = 14, height = 6.5)
  ggplot2::ggsave(file.path(output_dir, "sample_mds.pdf"), mds_figure, width = 14, height = 6.5)
  save_figure(qc_figure, output_dir, "sample_qc", width = 14, height = 10)
  ggplot2::ggsave(file.path(output_dir, "sample_qc.pdf"), qc_figure, width = 14, height = 10)
  list(csv = csv, png = setNames(file.path(output_dir, c("sample_mds.png", "sample_qc.png")),
    c("sample_mds", "sample_qc")), rds = rds)
}
