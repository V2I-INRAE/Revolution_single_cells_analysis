# Publication UMAP of saved clusters; no marker-expression panels.
# Called by R/plots/publications/main.R.
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})
source("R/utils/save_plots.R")

plot_umap_for_pub <- function(obj, cluster_col, output_dir) {
  metadata <- obj[[]]
  reduction <- obj@misc$clustering$umap
  coords <- Embeddings(obj, reduction)
  ids <- as.character(metadata[[cluster_col]])
  levels <- sort(unique(ids), method = "radix")
  cluster_colours <- viz_palette[c(4, 2, 3, 1, 13, 6, 10, 31)]
  if (length(levels) > length(cluster_colours)) {
    stop("Not enough distinct cluster colours")
  }
  df <- data.frame(
    x = coords[, 1],
    y = coords[, 2],
    cluster = factor(paste("Cluster", ids), levels = paste("Cluster", levels))
  )
  centres <- aggregate(cbind(x, y) ~ cluster, df, median)
  colours <- setNames(cluster_colours[seq_along(levels)], levels(df$cluster))
  set.seed(1234)
  df <- df[sample.int(nrow(df)), , drop = FALSE]
  p <- ggplot(df, aes(x, y, colour = cluster)) +
    geom_point(size = 0.35, stroke = 0) +
    scale_colour_manual(values = colours, name = NULL) +
    ggrepel::geom_text_repel(
      data = centres,
      aes(label = cluster),
      colour = "black",
      size = 5,
      seed = 1234,
      point.size = NA,
      min.segment.length = Inf,
      show.legend = FALSE
    ) +
    coord_fixed() +
    labs(
      title = paste("Clusters -", obj@misc$subcluster$subset_label),
      subtitle = paste(format(nrow(df), big.mark = ","), "cells"),
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    theme_classic(base_size = 13) +
    theme(legend.position = "bottom") +
    guides(colour = guide_legend(override.aes = list(size = 3)))
  p <- ggrastr::rasterise(p, layers = "Point", dpi = 600)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(output_dir, "UMAP.pdf")
  ggsave(out, p, width = 7, height = 7, units = "in", bg = "white")
  invisible(out)
}
