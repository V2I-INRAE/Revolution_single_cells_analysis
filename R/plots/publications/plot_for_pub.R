# Publication UMAP of saved clusters; no marker-expression panels.
# Usage: Rscript R/plots/publications/plot_for_pub.R <subclustered.rds> <cluster_column> [output_dir]
suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
source("R/utils/save_plots.R")

plot_for_pub <- function(input, cluster_col, output_dir = NULL) {
  if (is.null(output_dir)) {
    parts <- strsplit(normalizePath(input, mustWork = TRUE), .Platform$file.sep)[[1L]]
    i <- match("subsets", parts)
    if (is.na(i) || length(parts) < i + 3L) stop("Cannot derive subset output directory; supply output_dir")
    output_dir <- file.path("results", parts[i - 1L], "clusters_plots", "subset", parts[i + 1L], parts[i + 2L], "publications_figures")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  stem <- paste0("clusters_umap_", cluster_col)
  out <- file.path(output_dir, paste0(stem, ".pdf"))
  info_path <- file.path(output_dir, paste0(stem, ".INFO"))
  if (file.exists(out) || file.exists(info_path)) stop("Publication output already exists: ", stem)
  info <- c("# Publication cluster UMAP", "Status: running", paste("Started:", Sys.time()),
    paste("Input (read-only):", input), paste("Cluster column:", cluster_col),
    "Saved UMAP coordinates and cluster identities; no recomputation or marker-expression panels.",
    "Seed: 1234 for point order and label placement; points rasterized at 600 dpi.")
  writeLines(info, info_path)
  completed <- FALSE
  on.exit({
    info[info == "Status: running"] <- if (completed) "Status: completed" else "Status: failed/incomplete"
    writeLines(c(info, paste("Finished:", Sys.time())), info_path)
  }, add = TRUE)
  input_info <- file.info(input)[, c("size", "mtime")]
  obj <- readRDS(input)
  metadata <- obj[[]]
  reduction <- obj@misc$clustering$umap
  if (is.null(reduction) || !reduction %in% Reductions(obj) || !cluster_col %in% names(metadata) ||
      anyNA(metadata[[cluster_col]]) || !identical(rownames(metadata), colnames(obj))) stop("Saved UMAP or cluster metadata invalid")
  coords <- Embeddings(obj, reduction)
  if (!identical(rownames(coords), colnames(obj))) stop("UMAP and cell IDs are not aligned")
  ids <- as.character(metadata[[cluster_col]])
  levels <- sort(unique(ids), method = "radix")
  if (length(levels) > length(viz_palette)) stop("Not enough cluster colours")
  df <- data.frame(x = coords[, 1], y = coords[, 2], cluster = factor(paste("Cluster", ids), levels = paste("Cluster", levels)))
  centres <- aggregate(cbind(x, y) ~ cluster, df, median)
  colours <- setNames(viz_palette[seq_along(levels)], levels(df$cluster))
  set.seed(1234)
  df <- df[sample.int(nrow(df)), , drop = FALSE]
  p <- ggplot(df, aes(x, y, colour = cluster)) + geom_point(size = 0.8, stroke = 0) +
    scale_colour_manual(values = colours, name = NULL) +
    ggrepel::geom_text_repel(data = centres, aes(label = cluster), colour = "black", size = 5,
      seed = 1234, point.size = NA, min.segment.length = Inf, show.legend = FALSE) +
    coord_fixed() + labs(title = paste("Clusters -", obj@misc$subcluster$subset_label),
      subtitle = paste(format(nrow(df), big.mark = ","), "cells"), x = "UMAP 1", y = "UMAP 2") +
    theme_classic(base_size = 13) + theme(legend.position = "bottom")
  p <- ggrastr::rasterise(p, layers = "Point", dpi = 600)
  ggsave(out, p, width = 7, height = 7, units = "in", bg = "white")
  stopifnot(identical(rawToChar(readBin(out, "raw", 5L)), "%PDF-"),
    identical(input_info, file.info(input)[, c("size", "mtime")]))
  info <- c(info, paste("Reduction:", reduction, "; cells:", nrow(df), "; clusters:", paste(levels, collapse = ", ")),
    paste("Output:", out), "Validation: PDF header checked; input file size and modification time unchanged.")
  completed <- TRUE
  invisible(out)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2L || length(args) > 3L) stop("Expected input object, cluster column and optional output directory")
  plot_for_pub(args[1], args[2], if (length(args) == 3L) args[3] else NULL)
}
