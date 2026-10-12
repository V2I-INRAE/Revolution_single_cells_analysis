suppressPackageStartupMessages(library(ggplot2))

marker_colours <- c(
  "#3A5BA0", "#F5A623", "#1ABC9C", "#90141A", "#EBA5AB", "#64B024",
  "#2E5111", "#8EC9EB", "#5C6C6B", "#D35400", "#8C6D31", "#9B59B6",
  "#E46571", "#C5A9D8", "#5A8F5A", "#8B4513", "#4E9A9A", "#6B8E23",
  "#D4753C", "#97767A", "#2E86C1", "#203161", "#BEACAE", "#C44E52",
  "#7B5EA7", "#E8A838", "#46878F", "#B07AA1", "#D98880", "#86714D",
  "#6E2F84", "#7BC657", "#708090", "#B8860B", "#40569A", "#931419"
)

save_marker_plot <- function(plot, stem, directory, width, height) {
  path <- file.path(directory, paste0(stem, ".png"))
  ggsave(path, plot, width = width, height = height, dpi = 300, bg = "white")
}

plot_markers_dotplot <- function(alldata, top5, sel.clust, directory) {
  dots <- DotPlot(alldata, features = rev(as.character(unique(top5$gene))), group.by = sel.clust, assay = "RNA") + coord_flip()
  dots <- dots + labs(
    title = paste(sel.clust, "| top 5 markers per cluster"),
    subtitle = "Dot size: detection; colour: scaled average expression",
    caption = "All cells; RNA LogNormalize expression. Exploratory cluster markers, not condition DE.",
    x = NULL, y = "Cluster"
  ) + theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 9))
  save_marker_plot(dots, "dotplot_top5", directory,
    width = max(12, nlevels(Idents(alldata)) * 0.35),
    height = max(6, length(unique(top5$gene)) * 0.22 + 2)
  )
}

plot_top_markers_heatmap <- function(alldata, top20, sel.clust, directory) {
  features <- as.character(unique(top20$gene))
  colours <- setNames(marker_colours[seq_along(levels(Idents(alldata)))], levels(Idents(alldata)))
  pages <- split(features, ceiling(seq_along(features) / 50))
  for (page in seq_along(pages)) {
    genes <- pages[[page]]
    heatmap <- Seurat::DoHeatmap(alldata,
      features = genes, assay = "RNA", slot = "scale.data", group.by = sel.clust,
      group.colors = colours, label = FALSE, raster = TRUE, lines.width = 1,
      disp.min = -2.5, disp.max = 2.5
    ) + scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, limits = c(-2.5, 2.5), name = "Scaled\nexpression"
    ) + labs(
      title = paste(sel.clust, "| top 20 markers per cluster | page", page),
      subtitle = "All cells, grouped by cluster; saved regressed/scaled RNA expression",
      caption = "Saved scale.data; exploratory cluster markers, not condition DE."
    ) + theme(axis.text.y = element_text(size = 9), legend.position = "right")
    save_marker_plot(heatmap, sprintf("heatmap_top20_%02d", page),
      directory, width = 20, height = max(7, length(genes) * 0.22 + 2))
  }
}

plot_marker_features <- function(alldata, features, sel.clust, reduction, directory, label) {
  present <- features[features %in% rownames(alldata[["RNA"]])]
  missing <- setdiff(features, present)
  if (length(missing) > 0L) {
    warning(label, ": not in object, skipped: ", paste(missing, collapse = ", "))
  }
  if (length(present) == 0L) return(invisible(NULL))

  # RNA data is already log-normalized. Do not log it again or use scale.data.
  panels <- lapply(present, function(gene) {
    scplotter::FeatureStatPlot(alldata, features = gene,
      reduction = reduction, plot_type = "dim", assay = "RNA", layer = "data",
      center_zero = FALSE, theme = "theme_blank"
    ) + ggplot2::labs(title = gene)
  })
  combined <- patchwork::wrap_plots(panels, ncol = 3) + patchwork::plot_annotation(
    title = paste(label, "feature panel |", sel.clust),
    caption = paste("All cells; RNA LogNormalize data on", reduction,
      "(no z-scores). Configured panels, not cluster annotations.")
  )
  save_marker_plot(combined, paste0("features_", label), directory,
    width = 3 * 5.2, height = ceiling(length(panels) / 3) * 4.8)
}
