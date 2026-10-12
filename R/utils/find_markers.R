# Marker entry point: compute cluster markers, save tables, then plot.
# Usage: Rscript R/utils/find_markers.R <clustered.rds> <resolution>
suppressPackageStartupMessages({
  library(Seurat)
  library(presto)
})
source("R/config.R")
source("R/plots/markers_plots.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript R/utils/find_markers.R <data/job-XXX/clustered.rds> <resolution>")
}
input <- args[[1L]]
resolution <- as.numeric(args[[2L]])
job_id <- basename(dirname(input))
cfg <- config$markers
sel.clust <- paste0(cfg$cluster_prefix, resolution)
directory <- file.path("results", job_id, "clusters_plots", "markers", as.character(resolution))
marker_file <- file.path(directory, "markers_significant.csv")

# Protect saved tables and figures from being overwritten.
if (length(list.files(directory, pattern = "\\.(csv|png)$")) > 0L) {
  stop("Marker outputs already exist; preserve them before running again.")
}
dir.create(directory, recursive = TRUE, showWarnings = FALSE)

# Keep normalized expression, saved regressed values, and the Harmony UMAP.
message("Loading marker input: ", input, " (resolution ", resolution, ")")
alldata <- readRDS(input)
DefaultAssay(alldata) <- "RNA"
layers <- c(Layers(alldata[["RNA"]], search = "^data($|\\.)"), "scale.data")
alldata <- DietSeurat(alldata, assays = "RNA", layers = layers,
  dimreducs = cfg$reduction, graphs = NULL, misc = FALSE)
alldata <- JoinLayers(alldata, assay = "RNA", layers = "^data($|\\.)", new = "data")
Idents(alldata) <- alldata[[sel.clust, drop = TRUE]]

# Discovery uses log-normalized expression, not regression residuals.
set.seed(config$run$seed)
markers <- FindAllMarkers(alldata, assay = "RNA", slot = "data",
  test.use = cfg$test_use, only.pos = cfg$only_pos,
  min.pct = cfg$min_pct, logfc.threshold = cfg$logfc_threshold,
  min.diff.pct = -Inf, max.cells.per.ident = Inf, return.thresh = Inf,
  verbose = FALSE)
markers <- markers |>
  dplyr::mutate(cluster = factor(cluster, levels = levels(Idents(alldata))),
    pct_difference = pct.1 - pct.2) |>
  dplyr::arrange(cluster, dplyr::desc(avg_log2FC), p_val_adj, gene)
significant <- dplyr::filter(markers, p_val_adj < cfg$adjusted_p_cutoff)
top5 <- significant |> dplyr::group_by(cluster) |> dplyr::slice_head(n = 5) |> dplyr::ungroup()
top20 <- significant |> dplyr::group_by(cluster) |> dplyr::slice_head(n = 20) |> dplyr::ungroup()

# Save marker results before rendering figures.
write.csv(markers, file.path(directory, "markers_prefiltered.csv"), row.names = FALSE)
write.csv(significant, marker_file, row.names = FALSE)
for (n in c(5, 20, 25)) {
  top <- significant |> dplyr::group_by(cluster) |> dplyr::slice_head(n = n) |> dplyr::ungroup()
  write.csv(top, file.path(directory, paste0("markers_top", n, ".csv")), row.names = FALSE)
}
counts <- table(Idents(alldata))
summary <- data.frame(cluster = names(counts), n_cells = as.integer(counts),
  n_prefiltered = as.integer(table(markers$cluster)),
  n_significant = as.integer(table(significant$cluster)),
  n_adjusted_p_zero = as.integer(table(markers$cluster[markers$p_val_adj == 0])))
composition <- as.data.frame(table(sample = alldata$sample,
  cluster = alldata[[sel.clust, drop = TRUE]]), responseName = "n_cells")
write.csv(summary, file.path(directory, "cluster_summary.csv"), row.names = FALSE)
write.csv(composition, file.path(directory, "sample_cluster_counts.csv"), row.names = FALSE)

# Default DotPlot, saved regressed heatmaps, and log-normalized feature panels.
plot_markers_dotplot(alldata, top5, sel.clust, directory)
plot_top_markers_heatmap(alldata, top20, sel.clust, directory)
for (panel in names(cfg$feature_panels)) {
  plot_marker_features(alldata, cfg$feature_panels[[panel]], sel.clust,
    cfg$reduction, directory, panel)
}
message("Marker tables and plots completed: ", directory)
