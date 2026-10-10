suppressPackageStartupMessages(library(ggplot2))

clustering_colours <- c(
  "#3A5BA0", "#F5A623", "#1ABC9C", "#90141A", "#EBA5AB", "#64B024",
  "#2E5111", "#8EC9EB", "#5C6C6B", "#D35400", "#8C6D31", "#9B59B6",
  "#E46571", "#C5A9D8", "#5A8F5A", "#8B4513", "#4E9A9A", "#6B8E23",
  "#D4753C", "#97767A", "#2E86C1", "#203161", "#BEACAE", "#C44E52",
  "#7B5EA7", "#E8A838", "#46878F", "#B07AA1", "#D98880", "#86714D",
  "#6E2F84", "#7BC657", "#708090", "#B8860B", "#40569A", "#931419",
  "#C04736", "#F09F2E", "#D83746", "#DA5E28", "#C24935", "#E2DBA4",
  "#C0D0CB", "#BDD2CB", "#DDDDDB"
)

save_clustering_plot <- function(plot, stem, output_dir, width, height) {
  ggsave(file.path(output_dir, paste0(stem, ".png")), plot,
    width = width, height = height, dpi = 300, bg = "white")
}

plot_cluster_umap_panel <- function(obj, reduction, column, title, seed) {
  labels <- levels(obj[[]][, column])
  stopifnot("More clusters than palette colours; choose a palette policy before plotting" =
    length(labels) <= length(clustering_colours))
  scplotter::CellDimPlot(obj, reduction = reduction, group_by = column,
    palcolor = setNames(clustering_colours[seq_along(labels)], labels),
    label = TRUE, label_insitu = TRUE, label_repel = TRUE,
    label_size = 3, label_fg = "black", label_bg = "white",
    label_repulsion = 0.2, label_pt_size = 0, pt_size = 3, pt_alpha = 1,
    raster = TRUE, raster_dpi = c(1200, 1200),
    seed = seed, order = "random", show_stat = FALSE,
    legend.position = "none", theme = "theme_blank",
    title = title, xlab = "UMAP 1", ylab = "UMAP 2")
}

plot_cluster_umap_comparisons <- function(obj, reductions, output_dir) {
  config <- obj@misc$clustering
  for (i in seq_len(nrow(config$partitions))) {
    partition <- config$partitions[i, ]
    panels <- lapply(seq_along(reductions), function(j) {
      plot_cluster_umap_panel(obj, reductions[[j]], partition$column,
        names(reductions)[j], config$settings$seed)
    })
    plot <- patchwork::wrap_plots(panels, ncol = 2) +
      patchwork::plot_annotation(
        title = paste("CCA — fixed Leiden clusters, resolution", partition$resolution),
        caption = paste(
          "Same cells and cluster assignments; only UMAP coordinates change. Colours match within this resolution only.",
          "Axes scale independently; apparent distances and separation are not evidence of biological improvement.",
          sep = "\n"))
    save_clustering_plot(plot, paste0("umap_clusters_res_", partition$resolution),
      output_dir, width = 14, height = 7)
    rm(panels, plot)
  }
  invisible(NULL)
}

plot_cluster_resolutions <- function(obj, output_dir, plot_resolutions) {
  config <- obj@misc$clustering
  partitions <- config$partitions[config$partitions$resolution %in% plot_resolutions, ]
  stopifnot("Saved UMAP is required for plotting" =
    !is.null(config$umap) && config$umap %in% Seurat::Reductions(obj))
  panels <- lapply(seq_len(nrow(partitions)), function(i) {
    plot_cluster_umap_panel(obj, config$umap, partitions$column[i],
      paste(config$method, "— resolution", partitions$resolution[i]), config$settings$seed)
  })
  umap <- patchwork::wrap_plots(panels, ncol = 2) +
    patchwork::plot_annotation(caption =
      "Same coordinates in every panel. Cluster numbers and colours do not establish correspondence between resolutions.")
  save_clustering_plot(umap, "umap_resolutions", output_dir,
    width = 14, height = 6 * ceiling(length(panels) / 2))

  # scplotter's tree palette encodes resolution, not cluster identity.
  tree_obj <- obj
  for (column in setdiff(config$partitions$column, partitions$column)) {
    tree_obj[[column]] <- NULL
  }
  tree <- scplotter::ClustreePlot(tree_obj, prefix = paste0(config$method, "_snn_res."),
    palcolor = clustering_colours[seq_len(nrow(partitions))],
    edge_palette = "viridis", seed = config$settings$seed,
    title = paste(config$method, "— cluster transitions")) +
    labs(caption = "Node colour: resolution; labels: cluster IDs. Transitions are not a resampling stability test.")
  max_clusters <- max(vapply(partitions$column,
    function(column) length(unique(obj[[]][, column])), integer(1)))
  save_clustering_plot(tree, "clustree", output_dir,
    width = max(12, max_clusters * 0.45), height = 9)
  invisible(list(umap = umap, tree = tree))
}
