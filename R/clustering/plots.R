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
  for (extension in c("png", "svg")) {
    ggsave(file.path(output_dir, paste0(stem, ".", extension)), plot,
      width = width, height = height, dpi = 300, bg = "white")
  }
}

plot_cluster_resolutions <- function(obj, output_dir, plot_resolutions) {
  config <- obj@misc$clustering
  partitions <- config$partitions[config$partitions$resolution %in% plot_resolutions, ]
  panels <- lapply(seq_len(nrow(partitions)), function(i) {
    column <- partitions$column[i]
    labels <- levels(obj[[]][, column])
    stopifnot("More clusters than palette colours; choose a palette policy before plotting" =
      length(labels) <= length(clustering_colours))
    scplotter::CellDimPlot(obj, reduction = config$umap, group_by = column,
      palcolor = setNames(clustering_colours[seq_along(labels)], labels),
      label = TRUE, label_insitu = TRUE, label_repel = TRUE,
      label_size = 3, label_fg = "black", label_bg = "white",
      label_repulsion = 0.2, label_pt_size = 0, pt_size = 3, pt_alpha = 1,
      raster = TRUE, raster_dpi = c(1200, 1200),
      seed = config$settings$seed, order = "random", show_stat = FALSE,
      legend.position = "none", theme = "theme_blank",
      title = paste(config$method, "— resolution", partitions$resolution[i]),
      xlab = "UMAP 1", ylab = "UMAP 2")
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

plot_clustering_comparison <- function(runs, agreement, output_dir, plot_resolutions) {
  summary <- do.call(rbind, lapply(runs, function(run) {
    data.frame(method = run$method, run$diagnostics$summary)
  }))
  scores <- do.call(rbind, lapply(runs, function(run) {
    data.frame(method = run$method, run$diagnostics$silhouettes)
  }))
  summary <- summary[summary$resolution %in% plot_resolutions, ]
  scores <- scores[scores$resolution %in% plot_resolutions, ]
  agreement <- agreement[agreement$resolution_a %in% plot_resolutions &
    agreement$resolution_b %in% plot_resolutions, ]
  colours <- c(unintegrated = "#3A5BA0", harmony = "#D4753C", scvi = "#5A8F5A")
  metrics <- c(mean_silhouette = "Mean silhouette (cell-weighted)",
    mean_cluster_silhouette = "Mean of cluster silhouette means",
    n_clusters = "Number of clusters", min_cluster_size = "Smallest cluster (cells)",
    n_small_clusters = "Clusters with fewer than 10 cells",
    absent_clusters = "Clusters absent from diagnostic subset")
  panels <- lapply(names(metrics), function(metric) {
    plotthis::LinePlot(summary, x = "resolution", y = metric, group_by = "method",
      palcolor = colours, pt_size = 2, aspect.ratio = NULL,
      xlab = "Resolution", ylab = metrics[[metric]], title = metrics[[metric]])
  })
  overview <- patchwork::wrap_plots(panels, ncol = 2, guides = "collect") +
    patchwork::plot_annotation(caption =
      "Descriptive comparisons only. Different embedding geometries preclude selecting the best method by silhouette alone.")
  save_clustering_plot(overview, "metrics_by_resolution", output_dir, 14, 14)

  distribution <- plotthis::BoxPlot(scores, x = "resolution", y = "silhouette",
    group_by = "method", palcolor = colours, add_point = FALSE,
    add_errorbar = FALSE, comparisons = NULL, outlier.shape = NA,
    xlab = "Resolution", ylab = "Silhouette (shared diagnostic subset)",
    title = "Silhouette distributions", x_text_angle = 0)
  save_clustering_plot(distribution, "silhouette_distributions", output_dir, 12, 6)

  agreement$comparison <- paste(agreement$method_a, "vs", agreement$method_b)
  agreement$kind <- ifelse(agreement$method_a == agreement$method_b,
    "Adjacent resolutions (x = lower resolution)", "Methods at matching resolution")
  ari <- plotthis::LinePlot(agreement, x = "resolution_a", y = "adjusted_rand",
    group_by = "comparison", facet_by = "kind", facet_ncol = 1,
    palcolor = clustering_colours[1:6], aspect.ratio = NULL, line_width = 0,
    xlab = "Resolution", ylab = "Adjusted Rand index", title = "Partition agreement") +
    coord_cartesian(ylim = c(-1, 1))
  save_clustering_plot(ari, "partition_agreement", output_dir, 12, 9)
  invisible(overview)
}
