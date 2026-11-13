# Summary of saved Pan-human Azimuth predictions for each subset and resolution in
# annotation_config$summary: a 2 x 2 UMAP figure (your clusters | final label /
# medium | fine) and one combined CSV with per-cluster statistics, marker detection
# and the curated verdicts/notes in R/annotations/azimuth_notes.csv.
# Usage: Rscript R/annotations/summary.R   (inputs read-only; outputs replaced)
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})
source("R/annotations/config.R")
source("R/utils/save_plots.R")

cfg <- annotation_config$summary
azimuth <- annotation_config$azimuth
prefix <- azimuth$column_prefix
# Same cluster colours as the publication UMAPs (R/plots/publications/plot_umap_for_pub.R).
cluster_colours <- viz_palette[c(4, 2, 3, 1, 13, 6, 10, 31)]
# Label colours: the house palette without its grey tones and darkest colours, so
# no label can be confused with the grey special figure labels below.
special_labels <- c("Not mapped (model table)", "Inconsistent hierarchy")
label_colours <- viz_palette[-c(7, 9, 22, 23, 33, 36, 42:45)]

# "Proposed Known Subset" column of a cluster identity report, markdown removed.
read_identities <- function(path) {
  rows <- grep("^\\| \\*\\*[0-9]+\\*\\* \\|", readLines(path), value = TRUE)
  fields <- strsplit(rows, "|", fixed = TRUE)
  identity <- vapply(fields, `[`, "", 4L)
  identity <- gsub("\\^\\{([^}]*)\\}", "\\1", gsub("\\*\\*|`", "", identity))
  setNames(trimws(identity), gsub("\\D", "", vapply(fields, `[`, "", 2L)))
}

# Percentage of cells with counts > 0, per gene and cluster, from every count layer.
detection <- function(object, genes, clusters) {
  genes <- unique(genes)
  detected <- do.call(cbind, lapply(Layers(object[["RNA"]], search = "^counts($|[.])"),
    function(layer) {
      m <- LayerData(object, assay = "RNA", layer = layer)
      out <- Matrix::Matrix(FALSE, length(genes), ncol(m), sparse = TRUE,
        dimnames = list(genes, colnames(m)))
      present <- intersect(genes, rownames(m))
      if (length(present)) out[present, ] <- m[present, , drop = FALSE] > 0
      out
    }))[, colnames(object), drop = FALSE]
  pct <- sapply(split(seq_along(clusters), clusters), function(i) {
    100 * Matrix::rowMeans(detected[, i, drop = FALSE])
  })
  pct <- matrix(pct, nrow = length(genes), dimnames = list(genes, sort(unique(clusters))))
  pct[!genes %in% rownames(object[["RNA"]]), ] <- NA
  pct
}

top_labels <- function(labels, n = 3L) {
  t <- sort(table(labels), decreasing = TRUE)[seq_len(n)]
  list(label = names(t), pct = round(100 * as.numeric(t) / length(labels), 1))
}

cluster_rows <- function(run, predictions, clusters, identities, mapped) {
  groups <- c(list(all = seq_along(clusters)),
    split(seq_along(clusters), factor(clusters,
      levels = as.character(sort(as.integer(unique(clusters)))))))
  do.call(rbind, lapply(names(groups), function(k) {
    i <- groups[[k]]
    top <- top_labels(predictions[[paste0(prefix, "label")]][i])
    data.frame(subset = run$subset, subset_job = run$job, resolution = run$resolution,
      cluster = k, n_cells = length(i),
      pct_of_subset = round(100 * length(i) / length(clusters), 1),
      your_identity = if (k == "all") "" else identities[[k]],
      azimuth_label_1 = top$label[1], azimuth_pct_1 = top$pct[1],
      azimuth_label_2 = top$label[2], azimuth_pct_2 = top$pct[2],
      azimuth_label_3 = top$label[3], azimuth_pct_3 = top$pct[3],
      pct_unassigned = round(100 * mean(predictions[[paste0(prefix, "label")]][i] == "Unassigned"), 1),
      median_confidence = round(median(predictions[[paste0(prefix, "confidence")]][i]), 3),
      pct_consistent_hierarchy = round(100 * mean(predictions[[paste0(prefix, "consistent")]][i]), 1),
      model_genes_mapped = mapped, stringsAsFactors = FALSE)
  }))
}

# Figure labels: "False" marks an inconsistent hierarchy; NA marks a valid path that
# is missing from the model's level table (e.g. plain AT2 at medium); rare labels
# become Other.
display_labels <- function(labels, min_fraction) {
  labels[is.na(labels)] <- special_labels[1]
  labels[labels == "False"] <- special_labels[2]
  counts <- table(labels)
  rare <- names(counts)[counts < min_fraction * length(labels)]
  other <- sprintf("Other (<%g%%)", 100 * min_fraction)
  labels[labels %in% rare] <- other
  kept <- names(sort(table(labels), decreasing = TRUE))
  specials <- c(other, special_labels)
  factor(labels, levels = c(setdiff(kept, specials), intersect(specials, kept)))
}

umap_panel <- function(df, column, title, colours, point_size, centres = NULL) {
  p <- ggplot(df, aes(x, y, colour = .data[[column]])) +
    geom_point(size = point_size, stroke = 0) +
    scale_colour_manual(values = colours, name = NULL, drop = FALSE) +
    coord_fixed() +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    theme_classic(base_size = 11) +
    theme(legend.position = "right", legend.text = element_text(size = 8),
      legend.key.size = grid::unit(3.5, "mm"), plot.title = element_text(face = "bold")) +
    guides(colour = guide_legend(override.aes = list(size = 2.5), ncol = 1))
  if (!is.null(centres)) {
    p <- p + ggrepel::geom_text_repel(data = centres, aes(x, y, label = cluster),
      inherit.aes = FALSE, size = 5, seed = 1234, point.size = NA,
      min.segment.length = Inf)
  }
  ggrastr::rasterise(p, layers = "Point", dpi = 300)
}

# One colour per label across the three Azimuth panels of a figure.
label_palette <- function(levels) {
  other <- grep("^Other ", levels, value = TRUE)
  regular <- setdiff(levels, c(other, special_labels))
  if (length(regular) > length(label_colours)) stop("Not enough label colours")
  c(setNames(label_colours[seq_along(regular)], regular),
    setNames(rep("grey80", length(other)), other),
    setNames(c("grey55", "grey25"), special_labels))
}

summary_figure <- function(object, predictions, cluster_col, run, overall) {
  coords <- Embeddings(object, object@misc$clustering$umap)
  clusters <- as.character(object[[cluster_col, drop = TRUE]])
  cluster_levels <- as.character(sort(as.integer(unique(clusters))))
  if (length(cluster_levels) > length(cluster_colours)) stop("Not enough cluster colours")
  df <- data.frame(x = coords[, 1], y = coords[, 2],
    cluster = factor(clusters, levels = cluster_levels))
  for (level in c("label", "medium", "fine")) {
    df[[level]] <- display_labels(predictions[[paste0(prefix, level)]], cfg$min_label_fraction)
  }
  centres <- aggregate(cbind(x, y) ~ cluster, df, median)
  set.seed(1234)
  df <- df[sample.int(nrow(df)), , drop = FALSE]
  point_size <- if (nrow(df) > 20000) 0.15 else 0.3
  colours <- label_palette(unique(c(levels(df$label), levels(df$fine), levels(df$medium))))
  panels <- list(
    umap_panel(df, "cluster", paste0("Your clusters (resolution ", run$resolution, ")"),
      setNames(cluster_colours[seq_along(cluster_levels)], cluster_levels), point_size, centres),
    umap_panel(df, "label", "Azimuth final label", colours, point_size),
    umap_panel(df, "medium", "Azimuth medium label", colours, point_size),
    umap_panel(df, "fine", "Azimuth fine label", colours, point_size))
  wrap_plots(panels, ncol = 2) + plot_annotation(
    title = sprintf("%s (job %s), resolution %s", run$subset, run$job, run$resolution),
    subtitle = sprintf(paste("%s cells | median confidence %.2f | consistent hierarchy %.0f%% |",
      "model genes mapped %s | labels below %g%% of cells shown as Other"),
      format(overall$n_cells, big.mark = ","), overall$median_confidence,
      overall$pct_consistent_hierarchy, overall$model_genes_mapped, 100 * cfg$min_label_fraction),
    caption = paste("Pan-human Azimuth v1, run locally; pig features enter through verified",
      "one-to-one Ensembl 116 orthologues, unmapped model genes are zero-filled.",
      "Exploratory cross-species labels."),
    theme = theme(plot.title = element_text(size = 15, face = "bold")))
}

panel_file <- Sys.glob(file.path(dirname(dirname(azimuth$python)), "lib", "python*",
  "site-packages", "panhumanpy", "_tools", azimuth$model_version,
  "inference_feature_panel", "inference_feature_panel.txt"))
if (length(panel_file) != 1L) stop("Model gene panel not found in the uv environment")
n_panel <- sum(nzchar(readLines(panel_file)))
notes <- read.csv(cfg$notes, colClasses = "character", check.names = FALSE)
rows <- list()
figures <- list()
for (r in seq_len(nrow(cfg$runs))) {
  run <- cfg$runs[r, ]
  data_dir <- file.path("data", paste0("job-", cfg$parent_job), "subsets", run$subset,
    paste0("job-", run$job))
  results_dir <- file.path("results", paste0("job-", cfg$parent_job), "clusters_plots",
    "subset", run$subset, paste0("job-", run$job), paste0("res", run$resolution))
  annotation_dir <- file.path(data_dir, "annotations", azimuth$output_folder)
  message("Summarising ", run$subset, " ", run$job, " res", run$resolution)

  object <- readRDS(file.path(data_dir, "subclustered.rds"))
  predictions <- readRDS(file.path(annotation_dir, "predictions.rds"))
  if (!identical(predictions$cell, colnames(object))) stop("Predictions do not match cells")
  mapping <- read.csv(file.path(annotation_dir, "gene_mapping.csv"))
  cluster_col <- paste0("harmony_snn_res.", run$resolution)
  clusters <- as.character(object[[cluster_col, drop = TRUE]])
  identities <- read_identities(file.path(results_dir, "Markers",
    paste0("cluster_identities_res", run$resolution, ".md")))
  if (!setequal(names(identities), unique(clusters))) stop("Identity report does not match clusters")

  rows[[r]] <- cluster_rows(run, predictions, clusters, identities,
    sprintf("%d/%d", sum(mapping$status == "used"), n_panel))

  run_notes <- notes[notes$subset == run$subset & notes$resolution == run$resolution, ]
  if (!setequal(run_notes$cluster, rows[[r]]$cluster) || anyDuplicated(run_notes$cluster)) {
    stop("Notes must cover every cluster once: ", run$subset, " res", run$resolution)
  }
  run_notes <- run_notes[match(rows[[r]]$cluster, run_notes$cluster), ]
  genes <- unlist(strsplit(run_notes$evidence_genes[nzchar(run_notes$evidence_genes)], ";"))
  pct <- detection(object, genes, clusters)
  rows[[r]]$verdict <- run_notes$verdict
  rows[[r]]$marker_detection_cluster <- vapply(seq_len(nrow(run_notes)), function(i) {
    g <- strsplit(run_notes$evidence_genes[i], ";")[[1]]
    if (!length(g)) return("")
    v <- pct[g, rows[[r]]$cluster[i]]
    paste(sprintf("%s %s", g, ifelse(is.na(v), "n/a", sprintf("%.0f%%", v))), collapse = "; ")
  }, "")
  rows[[r]]$note <- run_notes$note

  figure <- summary_figure(object, predictions, cluster_col, run, rows[[r]][1, ])
  figure_dir <- file.path(results_dir, "annotations", azimuth$output_folder)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(figure_dir, "UMAP_clusters_vs_azimuth.pdf"), figure,
    width = 15, height = 13, units = "in", bg = "white")
  figures[[r]] <- figure
  rm(object)
  gc()
}

summary_table <- do.call(rbind, rows)
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
csv <- file.path(cfg$output_dir, "azimuth_ph_v1_summary.csv")
write.csv(summary_table, csv, row.names = FALSE)
pdf_all <- file.path(cfg$output_dir, "azimuth_ph_v1_umaps.pdf")
grDevices::pdf(pdf_all, width = 15, height = 13, onefile = TRUE)
for (figure in figures) print(figure)
invisible(grDevices::dev.off())

# Reopen the table: per subset, cluster rows must add up to the "all" row.
reopened <- read.csv(csv, colClasses = c(cluster = "character", resolution = "character",
  subset_job = "character"), check.names = FALSE)
stopifnot(nrow(reopened) == nrow(summary_table), identical(reopened$note, summary_table$note),
  !anyNA(reopened$median_confidence), all(nzchar(reopened$verdict)))
for (key in unique(paste(reopened$subset, reopened$resolution))) {
  part <- reopened[paste(reopened$subset, reopened$resolution) == key, ]
  if (sum(part$n_cells[part$cluster != "all"]) != part$n_cells[part$cluster == "all"]) {
    stop("Cluster cell counts do not add up: ", key)
  }
}
message("Saved ", nrow(reopened), " rows to ", csv, " and ", length(figures),
  " figures (", pdf_all, " plus one per subset resolution folder).")
