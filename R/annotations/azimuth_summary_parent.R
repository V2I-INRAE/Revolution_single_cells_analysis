# Summary of full-parent Pan-human Azimuth predictions at resolution 0.7 only.
# One table (all cells + 24 clusters) and a two-panel UMAP: clusters | final labels.
# Coordinates come from the original parent, restricted to cleaned cell barcodes.
# No marker evidence, curated verdicts, reclustering or new embedding.
# Usage: Rscript R/annotations/azimuth_summary_parent.R
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})
source("R/annotations/config.R")
source("R/utils/save_plots.R")

cfg <- annotation_config$parent_summary
azimuth <- annotation_config$azimuth
prefix <- azimuth$column_prefix
cluster_col <- paste0("harmony_snn_res.", cfg$resolution)
annotation_dir <- file.path(dirname(cfg$input), "annotations", azimuth$output_folder)
inputs <- c(cfg$input, cfg$reference, file.path(annotation_dir, "predictions.rds"),
  file.path(annotation_dir, "gene_mapping.csv"))
input_info <- file.info(inputs)[, c("size", "mtime")]

message("Loading cleaned parent (read-only): ", cfg$input)
object <- readRDS(cfg$input)
cells <- colnames(object)
metadata <- object[[]]
clusters <- as.character(metadata[[cluster_col]])
stopifnot(length(cells) == 234183L, !anyDuplicated(cells),
  identical(rownames(metadata), cells), !anyNA(metadata$orig.ident),
  !anyNA(metadata$subset), length(clusters) == length(cells), !anyNA(clusters))
cluster_levels <- as.character(sort(as.integer(unique(clusters))))
stopifnot(length(cluster_levels) == 24L)
rm(object)
invisible(gc())

predictions <- readRDS(inputs[3])
stopifnot(identical(predictions$cell, cells), !anyDuplicated(predictions$cell))
labels <- predictions[[paste0(prefix, "label")]]
confidence <- predictions[[paste0(prefix, "confidence")]]
consistent <- predictions[[paste0(prefix, "consistent")]]
stopifnot(!anyNA(labels), all(nzchar(labels)), !anyNA(confidence),
  all(is.finite(confidence)), all(confidence >= 0 & confidence <= 1),
  !anyNA(consistent))

# Cleanup discarded the old reductions but retained the original cluster metadata.
# Release its expression matrices before loading the larger original parent.
message("Loading original UMAP (read-only): ", cfg$reference)
parent <- readRDS(cfg$reference)
at <- match(cells, colnames(parent))
stopifnot(!anyNA(at), !anyDuplicated(colnames(parent)),
  identical(clusters, as.character(parent[[cluster_col, drop = TRUE]])[at]))
reduction <- parent@misc$clustering$umap
coords <- Embeddings(parent, reduction)[cells, , drop = FALSE]
stopifnot(identical(rownames(coords), cells), ncol(coords) == 2L, all(is.finite(coords)))
message("Matched ", length(cells), " cells and 24 clusters to ", reduction,
  "; excluded ", ncol(parent) - length(cells), " cells removed by cleanup; no refitting.")
rm(parent)
invisible(gc())

panel_file <- Sys.glob(file.path(dirname(dirname(azimuth$python)), "lib", "python*",
  "site-packages", "panhumanpy", "_tools", azimuth$model_version,
  "inference_feature_panel", "inference_feature_panel.txt"))
if (length(panel_file) != 1L) stop("Model gene panel not found in the uv environment")
mapping <- read.csv(inputs[4])
mapped <- sprintf("%d/%d", sum(mapping$status == "used"),
  sum(nzchar(readLines(panel_file))))

groups <- c(list(all = seq_along(cells)),
  split(seq_along(cells), factor(clusters, levels = cluster_levels)))
summary_table <- do.call(rbind, lapply(names(groups), function(k) {
  i <- groups[[k]]
  top <- sort(table(labels[i]), decreasing = TRUE)[seq_len(3L)]
  subsets <- sort(table(metadata$subset[i]), decreasing = TRUE)
  data.frame(parent_job = cfg$parent_job, input_producer_job = cfg$input_producer_job,
    resolution = cfg$resolution, cluster = k, n_cells = length(i),
    pct_of_parent = round(100 * length(i) / length(cells), 1),
    coarse_subsets = if (k == "all") "" else paste(sprintf("%s %.1f%%",
      names(subsets), 100 * as.numeric(subsets) / length(i)), collapse = "; "),
    azimuth_label_1 = names(top)[1], azimuth_pct_1 = round(100 * as.numeric(top[1]) / length(i), 1),
    azimuth_label_2 = names(top)[2], azimuth_pct_2 = round(100 * as.numeric(top[2]) / length(i), 1),
    azimuth_label_3 = names(top)[3], azimuth_pct_3 = round(100 * as.numeric(top[3]) / length(i), 1),
    pct_unassigned = round(100 * mean(labels[i] == "Unassigned"), 1),
    median_confidence = round(median(confidence[i]), 3),
    pct_consistent_hierarchy = round(100 * mean(consistent[i]), 1),
    model_genes_mapped = mapped, stringsAsFactors = FALSE)
}))

# All final labels are retained: no rare-label threshold or Other grouping.
label_levels <- names(sort(table(labels), decreasing = TRUE))
df <- data.frame(x = coords[, 1], y = coords[, 2],
  cluster = factor(clusters, levels = cluster_levels),
  label = factor(labels, levels = label_levels))
centres <- aggregate(cbind(x, y) ~ cluster, df, median)
set.seed(1234)
df <- df[sample.int(nrow(df)), , drop = FALSE]

umap_panel <- function(column, title, colours, legend_columns, centres = NULL) {
  p <- ggplot(df, aes(x, y, colour = .data[[column]])) +
    geom_point(size = 0.15, stroke = 0) +
    scale_colour_manual(values = colours, name = NULL, drop = FALSE) +
    coord_fixed() + labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    theme_classic(base_size = 11) +
    theme(legend.position = "bottom", legend.text = element_text(size = 8),
      legend.key.size = grid::unit(3.5, "mm"), plot.title = element_text(face = "bold")) +
    guides(colour = guide_legend(override.aes = list(size = 2.5), ncol = legend_columns))
  if (!is.null(centres)) {
    p <- p + ggrepel::geom_text_repel(data = centres, aes(x, y, label = cluster),
      inherit.aes = FALSE, size = 5, seed = 1234, point.size = NA,
      min.segment.length = Inf, max.overlaps = Inf)
  }
  ggrastr::rasterise(p, layers = "Point", dpi = 300)
}

cluster_colours <- setNames(viz_palette[seq_along(cluster_levels)], cluster_levels)
# Give the frequent labels distinct house colours before extending the palette.
label_colours <- c(viz_palette[-c(7, 9, 22, 23, 33, 36, 42:45)],
  grDevices::hcl.colors(length(label_levels), "Dynamic"))
label_colours <- setNames(label_colours[seq_along(label_levels)], label_levels)
if ("Unassigned" %in% label_levels) label_colours["Unassigned"] <- "grey70"
figure <- wrap_plots(
  umap_panel("cluster", paste0("Your 24 clusters (resolution ", cfg$resolution, ")"),
    cluster_colours, 6, centres),
  umap_panel("label", "Azimuth final label", label_colours, 4), ncol = 2) +
  plot_annotation(title = paste0("Cleaned parent (job ", cfg$parent_job, ")"),
    subtitle = sprintf("%s cells | median confidence %.2f | consistent hierarchy %.0f%% | model genes mapped %s",
      format(length(cells), big.mark = ","), summary_table$median_confidence[1],
      summary_table$pct_consistent_hierarchy[1], mapped),
    caption = paste("Same original parent UMAP restricted to cleaned cells; no refitting. All final labels shown.",
      "Pan-human Azimuth v1, one-to-one pig-human orthologues; exploratory cross-species predictions.", sep = "\n"),
    theme = theme(plot.title = element_text(size = 15, face = "bold")))

dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
csv <- file.path(cfg$output_dir, "azimuth_ph_v1_summary.csv")
write.csv(summary_table, csv, row.names = FALSE)
pdf <- file.path(cfg$output_dir, "UMAP_clusters_vs_azimuth.pdf")
ggsave(pdf, figure, device = grDevices::cairo_pdf,
  width = 20, height = 9 + ceiling(length(label_levels) / 4) * 0.22,
  units = "in", limitsize = FALSE, bg = "white")

reopened <- read.csv(csv, colClasses = c(cluster = "character", resolution = "character",
  parent_job = "character", input_producer_job = "character"), check.names = FALSE)
stopifnot(nrow(reopened) == 25L, identical(reopened$cluster, names(groups)),
  identical(reopened$n_cells, summary_table$n_cells),
  sum(reopened$n_cells[reopened$cluster != "all"]) == length(cells),
  reopened$n_cells[reopened$cluster == "all"] == length(cells),
  identical(reopened$azimuth_label_1, summary_table$azimuth_label_1),
  identical(reopened$median_confidence, summary_table$median_confidence),
  identical(file.info(inputs)[, c("size", "mtime")], input_info))
message("Saved and reopened ", nrow(reopened), " rows: ", csv,
  "; saved two-panel UMAP with ", length(label_levels), " final labels: ", pdf,
  "; input files unchanged.")
