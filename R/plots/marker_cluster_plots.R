suppressPackageStartupMessages({
  library(ggplot2)
  library(scCustomize)
})

# Shared cluster dotplot + heatmap for parent/refit and subset marker runners.
# One implementation, one visual language; callers differ only in object
# preparation and the level_label ("Cluster" vs "Subcluster") in titles.
#
# Dotplot: normalized display (scale = FALSE) with a red colour scale on RNA
# normalized data, colour limits 0 to plot maximum. Heatmap: uniformly display-scaled
# top genes (ScaleData over exactly the plotted gene set, caller-supplied
# covariates) with an RdBu palette anchored at fixed symmetric limits.
# Both exports: PNG (300 dpi) and PDF. Heatmap cells capped per identity.

cluster_heatmap_max_cells_per_ident <- 500L
cluster_heatmap_disp_limits <- c(-2.5, 2.5)

plot_cluster_dotplot <- function(alldata, top5, sel.clust, directory,
  level_label = "Cluster") {
  dots <- scCustomize::DotPlot_scCustom(
    seurat_object = alldata,
    features = rev(as.character(unique(top5$gene))),
    group.by = sel.clust,
    assay = "RNA",
    scale = FALSE,
    flip_axes = TRUE
  ) + ggplot2::scale_colour_gradientn(
    colours = c("#FFF5F0", "#FCBBA1", "#CB181D", "#67000D"),
    limits = c(0, NA)
  ) + ggplot2::labs(
    title = paste(sel.clust, "| top 5 markers per", tolower(level_label)),
    subtitle = "Dot size: detection; colour: normalized average expression",
    caption = "RNA LogNormalize data, no gene-wise scaling. Colour: 0 to plot maximum.\nExploratory cluster markers, not condition DE.",
    x = NULL, y = level_label
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

plot_cluster_heatmap <- function(alldata, genes, sel.clust, directory, seed,
  level_label = "Cluster") {
  genes <- unique(as.character(genes))
  idents <- levels(Idents(alldata))
  set.seed(seed)
  cells <- unlist(lapply(idents, function(ident) {
    candidates <- colnames(alldata)[Idents(alldata) == ident]
    if (length(candidates) > cluster_heatmap_max_cells_per_ident) {
      sample(candidates, cluster_heatmap_max_cells_per_ident)
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
    size = 3,
    disp.min = cluster_heatmap_disp_limits[[1L]],
    disp.max = cluster_heatmap_disp_limits[[2L]]
  ) + ggplot2::scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(11L, "RdBu")),
    limits = cluster_heatmap_disp_limits,
    na.value = "white",
    name = "Expression"
  ) + ggplot2::labs(
    caption = paste("RNA display-scaled expression (scale.data),",
      "max", cluster_heatmap_max_cells_per_ident, "cells per", tolower(level_label), ";",
      "seed", seed, ". Exploratory", tolower(level_label), "markers, not condition DE.")
  )
  side <- max(10, length(idents) * 0.8 + 4, length(genes) * 0.06)
  save_figure(heatmap, directory, "heatmap_top25", width = side, height = side)
  ggplot2::ggsave(file.path(directory, "heatmap_top25.pdf"), plot = heatmap,
    width = side, height = side, units = "in", bg = "white")
}
