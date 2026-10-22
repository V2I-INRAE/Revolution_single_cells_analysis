# The driver sources R/utils/save_plots.R before calling this module.
plot_pseudobulk_de <- function(pb, output_dir, label, mapping,
                              fdr_cutoff = 0.05, logfc_cutoff = 0.5, seed = 1234) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  rd_bu <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
  sample_info <- pb$sample_info
  sample_info$condition <- factor(sample_info$condition, levels = c("Control", "P4", "N4"))
  sample_info <- sample_info[order(sample_info$condition, sample_info$pig), ]
  log2_cpm <- edgeR::cpm(pb$dge, log = TRUE, prior.count = 2)
  passes <- vapply(pb$results, function(tab) {
    genes <- tab$gene[tab$FDR <= fdr_cutoff & abs(tab$logFC) >= logfc_cutoff]
    rownames(log2_cpm) %in% genes
  }, logical(nrow(log2_cpm)))
  rownames(passes) <- rownames(log2_cpm)
  selected <- rownames(log2_cpm)[rowSums(passes) > 0]
  z_scores <- log2_cpm[selected, sample_info$sample, drop = FALSE]
  if (length(selected)) z_scores <- t(scale(t(z_scores)))
  row_hc <- if (length(selected) > 1) {
    stats::hclust(stats::dist(z_scores, method = "euclidean"), method = "complete")
  } else FALSE
  if (length(selected)) {
    annotation <- ComplexHeatmap::HeatmapAnnotation(
      Condition = sample_info$condition, Pig = sample_info$pig,
      col = list(Condition = setNames(viz_palette[1:3], levels(sample_info$condition)),
        Pig = setNames(viz_palette[seq_along(levels(sample_info$pig))], levels(sample_info$pig))),
      annotation_name_gp = grid::gpar(fontsize = 10))
    heatmap <- ComplexHeatmap::Heatmap(z_scores, name = "Log2 CPM\nrow z-score",
      col = circlize::colorRamp2(seq(-2.5, 2.5, length.out = 11), rd_bu),
      cluster_rows = row_hc, row_dend_reorder = FALSE,
      cluster_columns = FALSE, cluster_column_slices = FALSE,
      column_split = sample_info$condition, column_gap = grid::unit(3, "mm"),
      top_annotation = annotation, show_row_names = FALSE,
      column_names_rot = 60, column_names_gp = grid::gpar(fontsize = 10),
      rect_gp = grid::gpar(col = NA), use_raster = TRUE,
      raster_device = "agg_png", raster_quality = 4, raster_by_magick = FALSE,
      heatmap_legend_param = list(at = c(-2.5, 0, 2.5), labels = c("<= -2.5", "0", ">= 2.5")))
  }
  for (extension in c("png", "pdf")) {
    path <- file.path(output_dir, paste0("heatmap_groups.", extension))
    if (extension == "png") {
      grDevices::png(path, width = 11, height = 10, units = "in", res = 300, type = "cairo")
    } else {
      grDevices::pdf(path, width = 11, height = 10)
    }
    title <- paste0(label, ": ", format(length(selected), big.mark = ","),
      " DE genes across ", nrow(sample_info), " samples")
    if (length(selected)) {
      ComplexHeatmap::draw(heatmap, column_title = title,
        padding = grid::unit(c(14, 4, 5, 4), "mm"))
    } else {
      grid::grid.newpage()
      grid::grid.text(title, y = 0.9, gp = grid::gpar(fontsize = 16))
      grid::grid.text("No genes pass the joint FDR and absolute log2 fold-change thresholds.",
        gp = grid::gpar(fontsize = 12))
    }
    grid::grid.text(sprintf("FDR <= %g and |log2 FC| >= %g in at least one contrast. Pig effects retained.",
      fdr_cutoff, logfc_cutoff), x = 0.5, y = grid::unit(4, "mm"), gp = grid::gpar(fontsize = 9))
    grDevices::dev.off()
  }
  gene_order <- if (length(selected) > 1) rownames(z_scores)[row_hc$order] else selected
  csv <- setNames(file.path(output_dir, paste0(c("heatmap_genes", "heatmap_zscores",
    "volcano_data", "volcano_summary"), ".csv")),
    c("heatmap_genes", "heatmap_zscores", "volcano_data", "volcano_summary"))
  write.csv(data.frame(gene = gene_order, row_order = seq_along(gene_order),
    passes[gene_order, , drop = FALSE], check.names = FALSE), csv[["heatmap_genes"]], row.names = FALSE)
  write.csv(data.frame(gene = gene_order, z_scores[gene_order, , drop = FALSE],
    check.names = FALSE), csv[["heatmap_zscores"]], row.names = FALSE)

  class_colors <- c("Not significant" = "grey75", "Significant, small effect" = "black",
    "Down" = rd_bu[1], "Up" = rd_bu[11])
  volcano_data <- lapply(names(pb$results), function(comparison) {
    tab <- pb$results[[comparison]]
    safe_symbol <- mapping$safe_name[match(tab$gene, mapping$gene)]
    tab$display_gene <- ifelse(!is.na(safe_symbol) & nzchar(safe_symbol), safe_symbol, tab$gene)
    tab$contrast <- comparison
    tab$neg_log10_fdr <- -log10(tab$FDR)
    tab$class <- "Not significant"
    tab$class[tab$FDR <= fdr_cutoff] <- "Significant, small effect"
    tab$class[tab$FDR <= fdr_cutoff & tab$logFC <= -logfc_cutoff] <- "Down"
    tab$class[tab$FDR <= fdr_cutoff & tab$logFC >= logfc_cutoff] <- "Up"
    tab$class <- factor(tab$class, levels = names(class_colors))
    tab$label <- FALSE
    ranked <- order(tab$FDR, tab$PValue, -abs(tab$logFC), tab$gene)
    for (direction in c("Down", "Up")) {
      indices <- ranked[tab$class[ranked] == direction]
      tab$label[head(indices, 5)] <- TRUE
    }
    tab
  })
  names(volcano_data) <- names(pb$results)
  plot_data <- do.call(rbind, volcano_data)
  write.csv(plot_data, csv[["volcano_data"]], row.names = FALSE)
  class_counts <- as.data.frame(table(plot_data$contrast, plot_data$class), responseName = "n_genes")
  names(class_counts)[1:2] <- c("contrast", "class")
  write.csv(class_counts, csv[["volcano_summary"]], row.names = FALSE)
  x_limit <- max(abs(plot_data$logFC)) * 1.08
  y_limit <- max(plot_data$neg_log10_fdr) * 1.12
  titles <- c(P4_vs_Control = "P4 vs Control", N4_vs_Control = "N4 vs Control", P4_vs_N4 = "P4 vs N4")
  ratios <- c(P4_vs_Control = "P4 / Control", N4_vs_Control = "N4 / Control", P4_vs_N4 = "P4 / N4")
  for (comparison in names(volcano_data)) {
    tab <- volcano_data[[comparison]]
    counts <- table(tab$class)
    p <- plotthis::ScatterPlot(tab, x = "logFC", y = "neg_log10_fdr", color_by = "class",
      color_name = "", palcolor = class_colors, shape = 16,
      size_by = 0.9, alpha = 0.8, border_size = 0,
      theme = "theme_classic", theme_args = list(base_size = 12),
      legend.position = "bottom", combine = FALSE, seed = seed,
      title = paste0(label, ": ", titles[[comparison]]),
      subtitle = sprintf("Down: %s   |   Significant, small effect: %s   |   Up: %s",
        format(counts[["Down"]], big.mark = ","),
        format(counts[["Significant, small effect"]], big.mark = ","),
        format(counts[["Up"]], big.mark = ",")),
      xlab = paste0("Log2 fold change (", ratios[[comparison]], ")"),
      ylab = "-log10(BH FDR)")[[1]] +
      ggplot2::geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", color = "grey50") +
      ggplot2::geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed", color = "grey50") +
      ggrepel::geom_text_repel(data = tab[tab$label, ],
        ggplot2::aes(logFC, neg_log10_fdr, label = display_gene), inherit.aes = FALSE,
        seed = seed, size = 3.2, max.overlaps = Inf, min.segment.length = 0,
        nudge_x = ifelse(tab$logFC[tab$label] < 0, -0.7, 0.7), nudge_y = 0.3,
        bg.color = "white", bg.r = 0.12, box.padding = 0.8,
        force = 10, force_pull = 0.1, max.iter = 100000, max.time = Inf) +
      ggplot2::coord_cartesian(xlim = c(-x_limit, x_limit), ylim = c(0, y_limit)) +
      ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow = TRUE)) +
      ggplot2::labs(color = NULL, caption = sprintf("Dashed guides: log2 FC +/-%g and FDR %g. Small effect: |log2 FC| < %g.\nFDR-only DE results are unchanged; labels show up to five top genes per outer direction.",
        logfc_cutoff, fdr_cutoff, logfc_cutoff))
    save_figure(p, output_dir, paste0("volcano_", comparison), width = 8.5, height = 7)
    ggplot2::ggsave(file.path(output_dir, paste0("volcano_", comparison, ".pdf")), p,
      width = 8.5, height = 7, bg = "white")
  }
  stems <- c("heatmap_groups", paste0("volcano_", names(volcano_data)))
  list(csv = csv, png = setNames(file.path(output_dir, paste0(stems, ".png")), stems),
    n_heatmap_genes = length(selected))
}
