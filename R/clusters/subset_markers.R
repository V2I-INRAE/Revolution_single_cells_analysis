# Subset marker entry point: cluster markers for one subclustered subset,
# then a top-5 dotplot and a top-25 heatmap.
# Usage: Rscript R/clusters/subset_markers.R <data/job-XXX/subsets/<Folder>/<subrun>/subclustered.rds> <resolution>
suppressPackageStartupMessages({
  library(Seurat)
  library(presto)
})
source("R/config.R")
source("R/cell_cycle.R") # for rename_cc_features() used below
source("R/utils/save_plots.R")
source("R/plots/subset_markers_plots.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript R/clusters/subset_markers.R <data/job-XXX/subsets/<Folder>/<subrun>/subclustered.rds> <resolution>")
}
input <- args[[1L]]
resolution <- as.numeric(args[[2L]])
if (basename(input) != "subclustered.rds") {
  stop("Input must be a subclustered.rds file, got: ", input)
}
parts <- strsplit(normalizePath(input, mustWork = FALSE), .Platform$file.sep)[[1L]]
subsets_idx <- match("subsets", parts)
if (is.na(subsets_idx) || length(parts) < subsets_idx + 3L) {
  stop("Input path must look like data/job-XXX/subsets/<Folder>/<subrun>/subclustered.rds, got: ", input)
}
parent_job <- parts[subsets_idx - 1L]
folder <- parts[subsets_idx + 1L]
subrun <- parts[subsets_idx + 2L]
cfg <- config$markers
sel.clust <- paste0(cfg$cluster_prefix, resolution)
directory <- file.path("results", parent_job, "clusters_plots", "subset",
  folder, subrun, paste0("markers_res", resolution))

# Idempotent reruns: remove previous outputs of this script before writing new
# ones. ONLY files generated here are removed: the 7 marker CSVs, the dotplot
# and the heatmap. Sibling figures (.INFO, clustree, UMAPs, PCA plots) are
# never touched.
pipeline_csvs <- c("markers_prefiltered.csv", "markers_significant.csv",
  "markers_top5.csv", "markers_top20.csv", "markers_top25.csv",
  "cluster_summary.csv", "sample_cluster_counts.csv")
for (csv in pipeline_csvs) {
  csv_path <- file.path(directory, csv)
  if (file.exists(csv_path)) {
    file.remove(csv_path)
  }
}
for (stem in c("dotplot_top5.png", "heatmap_top25.png")) {
  stem_path <- file.path(directory, stem)
  if (file.exists(stem_path)) {
    file.remove(stem_path)
  }
}
dir.create(directory, recursive = TRUE, showWarnings = FALSE)

# Keep normalized and scaled expression for discovery and plotting. Unlike the
# parent entry point, scale.data layers are retained: DoHeatmap plots them.
message("Loading subset marker input: ", input, " (resolution ", resolution, ")")
alldata <- readRDS(input)
DefaultAssay(alldata) <- "RNA"
alldata <- JoinLayers(alldata, assay = "RNA", layers = "^data($|\\.)", new = "data")
scale_layers <- Layers(alldata[["RNA"]], search = "^scale\\.data($|\\.)")
if (length(scale_layers) > 1L) {
  alldata <- JoinLayers(alldata, assay = "RNA",
    layers = "^scale\\.data($|\\.)", new = "scale.data")
}
message("RNA layers kept: ",
  paste(Layers(alldata[["RNA"]], search = "^(data|scale\\.data)($|\\.)"),
    collapse = ", "))

# Translate ID-only features to gene names BEFORE discovery so the marker
# tables and the figures share one naming, as in the parent entry point.
renames_applied <- intersect(names(cfg$feature_renames), rownames(alldata[["RNA"]]))
alldata <- rename_cc_features(alldata, cfg$feature_renames)
message("Feature renames applied: ", length(renames_applied), "/", length(cfg$feature_renames))

if (!sel.clust %in% colnames(alldata[[]])) {
  stop("Cluster column not found in object: ", sel.clust)
}
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
top25 <- significant |> dplyr::group_by(cluster) |> dplyr::slice_head(n = 25) |> dplyr::ungroup()

# Save marker results before rendering figures.
write.csv(markers, file.path(directory, "markers_prefiltered.csv"), row.names = FALSE)
write.csv(significant, file.path(directory, "markers_significant.csv"), row.names = FALSE)
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

# Top-5 dotplot and top-25 heatmap.
plot_subset_dotplot(alldata, top5, sel.clust, directory)
plot_subset_heatmap(alldata, top25, sel.clust, directory, seed = config$run$seed)
message("Subset marker tables and plots completed: ", directory)
