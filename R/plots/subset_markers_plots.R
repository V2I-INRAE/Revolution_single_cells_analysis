suppressPackageStartupMessages({
  library(ggplot2)
  library(scCustomize)
})

# Maximum cells per identity drawn in the heatmap. Without a cap the larger
# subsets would render tens of thousands of columns: slow and visually
# merged. Recorded in each output .INFO.
subset_heatmap_max_cells_per_ident <- 500L

plot_subset_dotplot <- function(alldata, top5, sel.clust, directory) {
  dots <- scCustomize::DotPlot_scCustom(
    seurat_object = alldata,
    features = rev(as.character(unique(top5$gene))),
    group.by = sel.clust,
    flip_axes = TRUE,
    scale = FALSE
  ) + ggplot2::scale_colour_gradientn(
    colours = c("#FFF5F0", "#FCBBA1", "#CB181D", "#67000D"),
    limits = c(0, NA)
  ) + ggplot2::labs(
    title = paste(sel.clust, "| top 5 markers per subcluster"),
    subtitle = "Dot size: detection; colour: normalized average expression",
    caption = "Subset cells; RNA LogNormalize data, no gene-wise scaling. Colour: 0 to plot maximum.\nExploratory subcluster markers, not condition DE.",
    x = NULL, y = "Subcluster"
  ) + ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 9),
    axis.text.x = ggplot2::element_text(size = 9)
  )
  side <- max(12, nlevels(Idents(alldata)) * 0.35,
    length(unique(top5$gene)) * 0.22 + 2)
  save_figure(dots, directory, "dotplot_top5", width = side, height = side)
  ggplot2::ggsave(file.path(directory, "dotplot_top5.pdf"), plot = dots,
    width = side, height = side, units = "in", bg = "white")
}

plot_subset_heatmap <- function(
  alldata, top25, sel.clust, directory, seed
) {
  genes <- unique(as.character(top25$gene))
  idents <- levels(Idents(alldata))
  set.seed(seed)
  cells <- unlist(lapply(idents, function(ident) {
    candidates <- colnames(alldata)[Idents(alldata) == ident]
    if (length(candidates) > subset_heatmap_max_cells_per_ident) {
      sample(candidates, subset_heatmap_max_cells_per_ident)
    } else {
      candidates
    }
  }), use.names = FALSE)
  heatmap <- Seurat::DoHeatmap(alldata,
    features = genes,
    cells = cells,
    group.by = sel.clust,
    assay = "RNA",
    slot = "scale.data",
    label = TRUE,
    size = 3
  ) + ggplot2::labs(
    caption = paste("Subset cells; RNA scaled expression (scale.data),",
      "max", subset_heatmap_max_cells_per_ident, "cells per subcluster;",
      "seed", seed, ". Exploratory subcluster markers, not condition DE.")
  )
  side <- max(10, length(idents) * 0.8 + 4, length(genes) * 0.06)
  save_figure(heatmap, directory, "heatmap_top25", width = side, height = side)
  ggplot2::ggsave(file.path(directory, "heatmap_top25.pdf"), plot = heatmap,
    width = side, height = side, units = "in", bg = "white")
}
