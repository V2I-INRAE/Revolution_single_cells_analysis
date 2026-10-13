# Coarse-subset highlight entry point: one SVG per broad subset, each showing
# the subset's res-0.7 clusters in distinct palette colors on a grey background.
# Usage: Rscript R/clusters/identify_coarse_subsets.R <annotated_clustered.rds>
#
# Read-only: the object is loaded for plotting only; nothing is written back.
# Identity handling: Cluster_Highlight_Plot matches `cluster_name` against the
# active identity, so identities are set in memory to the res-0.7 cluster column.
# The file on disk is never modified.
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(scCustomize)
})
source("R/utils/save_plots.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript R/clusters/identify_coarse_subsets.R <data/job-XXX/annotated_clustered.rds>")
}
input <- args[[1L]]
cluster_col <- "harmony_snn_res.0.7"
reduction <- "umap_harmony"
seed <- 1234L # production run seed (results/job-44493014/.INFO)

# Owner-approved broad subsets (2026-10-13); every res-0.7 cluster covered.
subset_clusters <- list(
  LYM = c("1", "2", "4", "5", "7", "8"),
  MNP = c("3", "9", "13"),
  GRA = c("6", "15"),
  DCs = c("14", "18"),
  Endothelial = c("11", "12", "20"),
  Epithelial = c("16", "17", "21", "23", "24"),
  Stromal = c("10", "19"),
  Mixed = c("22")
)

out_dir <- file.path("results", basename(dirname(input)),
  "clusters_plots", "subset")

if (!file.exists(input)) stop("Input not found: ", input)
existing <- list.files(out_dir, pattern = "\\.svg$")
if (length(existing) > 0L) {
  stop("Subset SVGs already exist in ", out_dir, ": ",
    paste(existing, collapse = ", "),
    "; preserve them before running again.")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
in_mtime <- file.mtime(input)
in_size <- file.size(input)

message("Loading (read-only): ", input)
obj <- readRDS(input)
stopifnot(
  "Cluster column missing from object" = cluster_col %in% names(obj[[]]),
  "UMAP reduction missing from object" = reduction %in% Seurat::Reductions(obj),
  "Cell count differs from validated 256,550" = ncol(obj) == 256550L
)
present <- sort(unique(as.character(obj[[cluster_col, drop = TRUE]])))
stopifnot(
  "Subset map does not cover all clusters" =
    setequal(present, unlist(subset_clusters, use.names = FALSE))
)
message("All ", length(present), " res-0.7 clusters covered by the subset map.")

# In-memory only: cluster numbers become the active identity for matching.
SeuratObject::Idents(obj) <- obj[[cluster_col, drop = TRUE]]

for (subset in names(subset_clusters)) {
  clusters <- subset_clusters[[subset]]
  colors <- viz_palette[seq_along(clusters)]
  message("Subset ", subset, ": clusters ", paste(clusters, collapse = ","))
  set.seed(seed)
  # Fixed small point size with raster off: scCustomize auto-rasterizes above
  # 200k cells and the raster path ignores pt.size (byte-identical reruns).
  p <- scCustomize::Cluster_Highlight_Plot(seurat_object = obj,
    cluster_name = clusters, highlight_color = colors,
    background_color = "lightgray", reduction = reduction,
    pt.size = 0.2,
    raster = FALSE, label = TRUE) +
    ggplot2::ggtitle(paste0(subset, " \u2014 resolution 0.7 (",
      paste(clusters, collapse = ", "), ")")) +
    ggplot2::labs(caption =
      "Same coordinates as umap_resolutions.png (umap_harmony). Highlighted clusters in subset colors; other cells grey.")
  out_svg <- file.path(out_dir, paste0(subset, ".svg"))
  ggplot2::ggsave(out_svg, plot = p, width = 10, height = 8,
    units = "in", device = "svg", bg = "white")
  message("  Saved: ", out_svg, " (", file.size(out_svg), " bytes)")
  stopifnot("Empty figure written" = file.size(out_svg) > 0)
}

stopifnot(
  "Input file changed during run" =
    file.exists(input) && file.mtime(input) == in_mtime &&
      file.size(input) == in_size
)
message("Subset highlight SVGs completed: ", out_dir)
