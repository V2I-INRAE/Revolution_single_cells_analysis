# Publication UMAP: labelled clusters (left) and one marker per cluster in a
# three-row grid (right). PDF only; cells are rasterized (one image per UMAP),
# text and frames stay vector. Input RDS is read only.
# Usage (from project root):
#   Rscript R/clusters/plot_for_pub.R <subclustered.rds> <cluster_column> \
#     '1=PDPN,2=ACKR1,3=GPIHBP1,4=EDNRB,5=GJA5' <output_dir>
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork) })
source("R/utils/save_plots.R") # viz_palette

# Draw cells once at the panel's printed size; return them as one image layer.
raster_cells <- function(data, mapping, scale, size, lim, side, dpi) {
  path <- tempfile(fileext = ".png")
  on.exit(unlink(path))
  p <- ggplot(data, mapping) + geom_point(size = size, stroke = 0) + scale +
    coord_fixed(xlim = lim$x, ylim = lim$y, expand = FALSE) +
    theme_void() + theme(legend.position = "none")
  ragg::agg_png(path, width = side, height = side, units = "in", res = dpi, background = "white")
  print(p)
  grDevices::dev.off()
  annotation_raster(png::readPNG(path, native = TRUE), lim$x[1], lim$x[2], lim$y[1], lim$y[2])
}

plot_for_pub <- function(input, cluster_col, markers, output_dir,
  side = 5, dpi = 600, pt_main = 0.8, pt_marker = 0.42, feature_labels = character()) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(output_dir, "clusters_and_markers.pdf")
  info <- c("# Publication cluster/marker UMAP", "Status: running", paste("Date:", Sys.time()),
    paste("Script: R/clusters/plot_for_pub.R; input (read-only):", input),
    paste("Cluster column:", cluster_col),
    paste("Markers:", paste(names(markers), markers, sep = "=", collapse = ", ")),
    paste("Display aliases:", paste(names(feature_labels), feature_labels, sep = "=", collapse = ", ")),
    "Expression: stored RNA log-normalized data; one viridis scale over all markers; no clipping.",
    sprintf("Colours: viz_palette; cells rasterized at %d dpi (one image per UMAP); seed 1234; plain text labels.", dpi),
    paste("Output:", out))
  writeLines(info, file.path(output_dir, ".INFO"))
  obj <- readRDS(input)
  reduction <- obj@misc$clustering$umap
  obj <- JoinLayers(obj, assay = "RNA", layers = "^data($|[.])", new = "data")
  ids <- names(markers)
  clusters <- as.character(obj[[cluster_col, drop = TRUE]])
  stopifnot(setequal(ids, clusters), all(markers %in% rownames(obj[["RNA"]])),
    length(ids) <= length(viz_palette))
  expr <- as.matrix(LayerData(obj, assay = "RNA", layer = "data")[markers, , drop = FALSE])
  df <- data.frame(Embeddings(obj, reduction)[, 1:2], t(expr), check.names = FALSE,
    cluster = factor(paste("Cluster", clusters), levels = paste("Cluster", ids)))
  names(df)[1:2] <- c("x", "y")

  # Square windows; marker tiles get headroom for the gene name.
  r <- apply(df[c("x", "y")], 2, range)
  mid <- colMeans(r)
  s <- max(diff(r[, 1]), diff(r[, 2])) * 1.12
  main_lim <- list(x = mid[1] + c(-0.5, 0.5) * s, y = mid[2] + c(-0.5, 0.5) * s)
  ts <- max(diff(r[, 1]) * 1.08, diff(r[, 2]) * 1.24)
  tile_lim <- list(x = mid[1] + c(-0.5, 0.5) * ts, y = r[1, 2] - 0.03 * ts + c(0, ts))
  frame <- theme_void(base_size = 13) + theme(plot.margin = margin(0, 0, 0, 0),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.8))

  cols <- setNames(viz_palette[seq_along(ids)], levels(df$cluster))
  set.seed(1234)
  centres <- aggregate(cbind(x, y) ~ cluster, df, median)
  main <- ggplot(mapping = aes(x, y)) +
    raster_cells(df[sample(nrow(df)), ], aes(x, y, colour = cluster),
      scale_colour_manual(values = cols), pt_main, main_lim, side, dpi) +
    # Empty layer that only carries the cluster legend.
    geom_point(aes(colour = cluster), data.frame(x = NA_real_, y = NA_real_, cluster = centres$cluster), na.rm = TRUE) +
    scale_colour_manual(values = cols, limits = names(cols), name = NULL) +
    ggrepel::geom_text_repel(aes(label = cluster), centres, size = 5, seed = 1234,
      point.size = NA, min.segment.length = Inf, show.legend = FALSE) +
    annotate("text", -Inf, Inf, parse = TRUE, size = 5, hjust = -0.05, vjust = 1.3,
      label = sprintf("italic(n)~'= %s cells'", format(nrow(df), big.mark = ","))) +
    coord_fixed(xlim = main_lim$x, ylim = main_lim$y, expand = FALSE) + frame +
    guides(colour = guide_legend(nrow = ceiling(length(ids) / 5), override.aes = list(size = 3.2))) +
    theme(legend.position = "bottom", legend.justification = "left",
      legend.text = element_text(size = 13), legend.key.spacing.x = unit(1, "mm"))

  # Markers share one colour scale; high-expressing cells drawn on top.
  lim <- range(expr)
  tiles <- lapply(markers, function(g) ggplot() +
    raster_cells(df[order(df[[g]]), ], aes(x, y, colour = .data[[g]]),
      scale_colour_viridis_c(limits = lim), pt_marker, tile_lim, side / 3, dpi) +
    annotate("text", -Inf, Inf, label = if (g %in% names(feature_labels)) feature_labels[[g]] else g,
      fontface = "italic", size = 4.6, hjust = -0.06, vjust = 1.2) +
    coord_fixed(xlim = tile_lim$x, ylim = tile_lim$y, expand = FALSE) + frame)
  bar <- ggplot(data.frame(v = seq(lim[1], lim[2], length.out = 256)), aes(1, v, fill = v)) +
    geom_raster() + scale_fill_viridis_c(limits = lim, guide = "none") +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(position = "right", expand = c(0, 0)) +
    frame + theme(axis.text.y.right = element_text(margin = margin(l = 1.5, unit = "mm")),
      axis.ticks.y.right = element_line(), axis.ticks.length.y.right = unit(-0.06, "in"))

  # Main panel spans three tile rows; markers fill row-wise; colour bar at right.
  nc <- ceiling(length(markers) / 3)
  i <- seq_along(markers) - 1
  design <- c(area(1, 1, 3, 1), area(i %/% nc + 1, i %% nc + 3), area(1, nc + 4, 3, nc + 4))
  title <- obj@misc$subcluster$subset_label
  fig <- wrap_plots(c(list(main), tiles, list(bar)), design = design,
    widths = unit(c(side, 0.12, rep(side / 3, nc), 0.12, 0.17), "in"),
    heights = unit(rep(side / 3, 3), "in")) +
    plot_annotation(title = paste("scRNA-seq -", title, "cells"),
      theme = theme(plot.title = element_text(size = 17, hjust = 0.5)))

  info <- c(info, paste("Reduction:", reduction, "; cells:", nrow(df)),
    paste("Clusters plotted:", paste(ids, collapse = ", ")))
  writeLines(info, file.path(output_dir, ".INFO")) # stays "running" if rendering fails
  ggsave(out, fig, width = side * (1 + nc / 3) + 1.1,
    height = side + 0.75 + 0.33 * ceiling(length(ids) / 5), bg = "white")
  info[2] <- "Status: completed (PDF written)"
  writeLines(info, file.path(output_dir, ".INFO"))
  invisible(out)
}

if (sys.nframe() == 0L) {
  a <- commandArgs(trailingOnly = TRUE)
  pairs <- do.call(rbind, strsplit(strsplit(a[3], ",")[[1]], "="))
  plot_for_pub(a[1], a[2], setNames(pairs[, 2], pairs[, 1]), a[4])
}
