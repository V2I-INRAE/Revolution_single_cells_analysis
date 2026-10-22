# The driver sources R/utils/save_plots.R before calling this module.
source("R/DEG/publication_style.R")
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

  render_pseudobulk_diagnostics(mds_data, mds, nrow(pb$dge), output_dir, label, seed)
  list(csv = csv, png = setNames(file.path(output_dir, c("sample_mds.png", "sample_qc.png")),
    c("sample_mds", "sample_qc")), rds = rds)
}

# Rendering is separate so saved coordinates can be restyled without rerunning MDS.
render_pseudobulk_diagnostics <- function(mds_data, mds, n_genes, output_dir, label, seed = 1234) {
  sample_info <- mds_data
  condition_colors <- setNames(viz_palette[seq_along(levels(sample_info$condition))],
    levels(sample_info$condition))
  pig_colors <- setNames(viz_palette[seq_along(levels(sample_info$pig))], levels(sample_info$pig))
  axis_labels <- sprintf("Leading logFC\ndimension %d (%.1f%%)",
    1:2, 100 * mds$var.explained[1:2])
  # Label-only starting offsets for the four residual collisions observed at
  # final print size. Points, distances, sample identities and axes are untouched.
  label_offsets <- switch(label,
    "Endothelial cluster 4" = c("REVO29-P4" = 0.4, "REVO26-N4" = -0.3),
    "Epithelial cluster 5" = c("REVO31-N4" = 0.35, "REVO31-P4" = -0.2),
    "Granulocytes cluster 1" = c("REVO30-N4" = 0.35, "REVO26-P4" = -0.25),
    "Macrophage_Monocytes cluster 3" = c("REVO27-P4" = 0.7, "REVO31-N4" = -0.3),
    numeric())
  nudge_y <- unname(label_offsets[match(mds_data$sample, names(label_offsets))])
  nudge_y[is.na(nudge_y)] <- 0
  mds_plots <- lapply(c("condition", "pig"), function(group) {
    colors <- if (group == "condition") condition_colors else pig_colors
    plotthis::ScatterPlot(mds_data, x = "MDS1", y = "MDS2", color_by = group,
      color_name = if (group == "condition") "Condition" else "Pig", palcolor = colors,
      size_by = 3.2, alpha = 1, border_size = 0.3, legend.position = "bottom",
      aspect.ratio = NULL,
      theme = "theme_classic", theme_args = list(base_size = 10), combine = FALSE,
      title = if (group == "condition") "Samples by condition" else "Samples by pig",
      xlab = axis_labels[1], ylab = axis_labels[2], seed = seed)[[1]] +
      ggrepel::geom_text_repel(data = mds_data,
        ggplot2::aes(MDS1, MDS2, label = sub("^REVO", "", sample)), inherit.aes = FALSE,
        size = 9 / ggplot2::.pt, segment.size = 0.25, seed = seed, max.overlaps = Inf, min.segment.length = 0,
        point.size = 3.2, point.padding = 0.2, nudge_y = nudge_y,
        box.padding = 0.25, force = 5, force_pull = 0.1, max.iter = 1000000, max.time = Inf) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.25)) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.25)) +
      ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
        fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE)) +
      deg_theme() +
      ggplot2::theme(legend.direction = "horizontal", plot.title.position = "panel") +
      ggplot2::coord_equal()
  })
  mds_figure <- patchwork::wrap_plots(mds_plots, ncol = 2) + patchwork::plot_annotation(
    title = paste0(label, ": pseudobulk sample MDS"),
    subtitle = sprintf("%d samples from %d pigs. Identical coordinates; animal effects retained.",
      nrow(sample_info), length(unique(sample_info$pig))),
    caption = deg_wrap(paste0("TMM-normalized expression; leading logFC distances using the top 500 genes per sample pair.\n",
      "Sample labels omit the shared REVO prefix (e.g. 26-P4 = REVO26-P4).")),
    tag_levels = "A", theme = deg_theme())

  qc_data <- sample_info
  qc_data$sample <- factor(qc_data$sample,
    levels = with(sample_info, sample[order(pig, condition)]))
  qc_data$library_millions <- qc_data$library_size / 1e6
  qc_data$effective_library_millions <- qc_data$effective_library_size / 1e6
  metrics <- c(n_cells = "Cells", library_millions = "Library size (million counts)",
    tmm_factor = "TMM factor", effective_library_millions = "Effective library size (millions)")
  metric_plots <- lapply(names(metrics), function(metric) {
    p <- plotthis::BarPlot(qc_data, x = "sample", y = metric,
      fill_by = "condition", fill_name = "Condition", palcolor = condition_colors,
      x_text_angle = 90, aspect.ratio = NULL, y_min = 0,
      y_max = 1.12 * max(qc_data[[metric]]), theme = "theme_classic",
      theme_args = list(base_size = 10), combine = FALSE,
      xlab = "Sample (grouped by pig)", ylab = metrics[[metric]],
      title = metrics[[metric]], seed = seed)[[1]] + deg_theme() +
      ggplot2::theme(plot.title.position = "panel", axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
    if (metric == "tmm_factor") {
      p <- p + ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.25)
    }
    p
  })
  qc_figure <- (patchwork::wrap_plots(metric_plots, ncol = 2, guides = "collect") +
    patchwork::plot_annotation(
      title = paste0(label, ": sample QC"),
      subtitle = paste("All", nrow(sample_info), "samples retained. Libraries use the",
        format(n_genes, big.mark = ","), "genes retained by filterByExpr."),
      caption = deg_wrap("Effective library size = filtered library size x TMM factor. Dashed line marks TMM factor 1, not an exclusion cutoff."),
      tag_levels = "A", theme = deg_theme())) & ggplot2::theme(legend.position = "bottom")
  save_deg_figure(mds_figure, output_dir, "sample_mds", height = 135)
  save_deg_figure(qc_figure, output_dir, "sample_qc", height = 180)
}
