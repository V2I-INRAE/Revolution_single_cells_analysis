# Publication figures for one subset resolution, from saved results only (no refit).
# Run after R/DEG/main.R. Usage:
#   Rscript R/plots/publications/main.R --input=<subclustered.rds> --subset=<label>
#     --cluster-col=harmony_snn_res.X [--resolution-dir=<res0.X>]
# Reads <resolution-dir>/Markers/markers_top5_plus_canonical.csv, markers_top25.csv and
# <resolution-dir>/DEG. Writes <resolution-dir>/publication_figures/: UMAP, dotplot,
# heatmap_top25, DEG_<contrast> and gsea_<collection>/GSEA_<collection>_<contrast>.
# Every run replaces these figures.
source("R/DEG/config.R")
source("R/DEG/pseudobulk_de.R") # pseudobulk_comparisons()
source("R/plots/publications/publication_paths.R")
source("R/plots/publications/plot_umap_for_pub.R")
source("R/plots/publications/plot_top_pub.R")
source("R/plots/publications/plot_heatmap_for_pub.R")
source("R/plots/publications/plot_DEG_for_pub.R")
source("R/plots/publications/plot_gsea_overview.R")

run_publication_figures <- function(
  input,
  subset_label,
  cluster_col,
  resolution_dir = NULL
) {
  resolution <- sub("^harmony_snn_res[.]", "", cluster_col)
  if (is.null(resolution_dir)) {
    resolution_dir <- subset_resolution_dir(input, resolution)
  }
  deg_dir <- file.path(resolution_dir, "DEG")
  output_dir <- file.path(resolution_dir, "publication_figures")
  contrasts <- pseudobulk_comparisons()$contrast
  obj <- readRDS(input)

  # Idempotent: remove only the figures this script writes (cluster_XX/ is the old
  # per-cluster GSEA layout).
  unlink(Sys.glob(file.path(
    output_dir,
    c("UMAP.*", "dotplot.*", "heatmap_top25.*", "DEG_*_vs_*.*")
  )))
  unlink(
    Sys.glob(file.path(output_dir, c("gsea_*", "cluster_[0-9]*"))),
    recursive = TRUE
  )

  plot_umap_for_pub(obj, cluster_col, output_dir)
  plot_top_pub(
    obj,
    file.path(resolution_dir, "Markers", "markers_top5_plus_canonical.csv"),
    cluster_col,
    output_dir
  )
  plot_heatmap_for_pub(
    obj,
    file.path(resolution_dir, "Markers", "markers_top25.csv"),
    cluster_col,
    output_dir
  )
  for (contrast in contrasts) {
    plot_DEG_for_pub(deg_dir, output_dir, contrast, subset_label, resolution)
    for (collection in c("hallmark", "gobp", "reactome")) {
      plot_gsea_overview(
        deg_dir,
        file.path(output_dir, paste0("gsea_", collection)),
        collection,
        contrast,
        subset_label
      )
    }
  }
  invisible(output_dir)
}

if (sys.nframe() == 0L) {
  usage <- paste(
    "Usage: Rscript R/plots/publications/main.R --input=<subclustered.rds>",
    "--subset=<label> --cluster-col=harmony_snn_res.X [--resolution-dir=<res0.X>]"
  )
  args <- commandArgs(trailingOnly = TRUE)
  keys <- sub("^--([^=]+)=.*$", "\\1", args)
  allowed <- c("input", "subset", "cluster-col", "resolution-dir")
  if (
    !length(args) ||
      !all(grepl("^--[a-z-]+=.+$", args)) ||
      !all(keys %in% allowed)
  ) {
    stop(usage)
  }
  opt <- as.list(setNames(sub("^--[^=]+=", "", args), keys))
  run_publication_figures(
    opt$input,
    opt$subset,
    opt$`cluster-col`,
    opt$`resolution-dir`
  )
}
