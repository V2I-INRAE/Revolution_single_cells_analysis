# Clustering figures: one UMAP panel per Leiden resolution plus a clustree
# transition plot. Cluster numbers and colours do not establish identity.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

plot_cluster_umap_panel <- function(obj, reduction, column, title, seed) {
  labels <- levels(obj[[]][, column])
  stopifnot("More clusters than palette colours; choose a palette policy before plotting" =
    length(labels) <= length(viz_palette))
  scplotter::CellDimPlot(obj, reduction = reduction, group_by = column,
    palcolor = setNames(viz_palette[seq_along(labels)], labels),
    label = TRUE, label_insitu = TRUE, label_repel = TRUE,
    label_size = 3, label_fg = "black", label_bg = "white",
    label_repulsion = 0.2, label_pt_size = 0, pt_size = 3, pt_alpha = 1,
    raster = TRUE, raster_dpi = c(1200, 1200),
    seed = seed, order = "random", show_stat = FALSE,
    legend.position = "none", theme = "theme_blank",
    title = title, xlab = "UMAP 1", ylab = "UMAP 2")
}

plot_cluster_resolutions <- function(obj, plot_dir) {
  config <- obj@misc$clustering
  stopifnot(
    "Saved UMAP is required for plotting" =
      !is.null(config$umap) && config$umap %in% Reductions(obj)
  )
  panels <- lapply(seq_len(nrow(config$partitions)), function(i) {
    partition <- config$partitions[i, ]
    plot_cluster_umap_panel(obj, config$umap, partition$column,
      paste(config$method, "— resolution", partition$resolution),
      config$settings$seed)
  })
  umap <- patchwork::wrap_plots(panels, ncol = 2) +
    patchwork::plot_annotation(caption =
      "Same coordinates in every panel. Cluster numbers and colours do not establish correspondence between resolutions.")
  save_figure(umap, plot_dir, "umap_resolutions",
    width = 14, height = 6 * ceiling(length(panels) / 2))

  # scplotter's tree palette encodes resolution, not cluster identity.
  tree_obj <- obj
  tree <- scplotter::ClustreePlot(tree_obj,
    prefix = paste0(config$method, "_snn_res."),
    palcolor = viz_palette[seq_len(nrow(config$partitions))],
    edge_palette = "viridis", seed = config$settings$seed,
    title = paste(config$method, "— cluster transitions")) +
    labs(caption = "Node colour: resolution; labels: cluster IDs. Transitions are not a resampling stability test.")
  max_clusters <- max(vapply(config$partitions$column,
    function(column) length(unique(obj[[]][, column])), integer(1)))
  save_figure(tree, plot_dir, "clustree",
    width = max(12, max_clusters * 0.45), height = 9)
  invisible(list(umap = umap, tree = tree))
}
