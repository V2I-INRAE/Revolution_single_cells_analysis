# Explicit completed runs only; no normalization, clustering or embedding fitting.
args <- commandArgs(trailingOnly = TRUE)
plot_resolutions <- NULL
if (length(args) > 0L && startsWith(tail(args, 1L), "--plot-resolutions=")) {
  plot_resolutions <- as.numeric(strsplit(sub("^--plot-resolutions=", "", tail(args, 1L)), ",")[[1]])
  args <- head(args, -1L)
}
if (length(args) != 1L) {
  stop("Usage: Rscript R/clustering/plot_main.R <lognorm.rds> [--plot-resolutions=0.1,0.2,0.3,0.4]")
}

suppressPackageStartupMessages(library(Seurat))
source("R/clustering/plots.R")

obj <- readRDS(args[1])
config <- obj@misc$clustering
stopifnot("Input must be a clustered Seurat checkpoint" = !is.null(config))
if (is.null(plot_resolutions)) {
  plot_resolutions <- config$partitions$resolution
}
stopifnot("Plot resolutions must be unique saved resolutions" =
  !anyNA(plot_resolutions) && !anyDuplicated(plot_resolutions) &&
    length(plot_resolutions) > 0L &&
    all(plot_resolutions %in% config$partitions$resolution))

output_dir <- file.path("results", "clustering", config$method, config$run_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
message("Plotting ", normalizePath(args[1]), "; resolutions: ",
  paste(plot_resolutions, collapse = ", "), "; output: ", output_dir)
plot_cluster_resolutions(obj, output_dir, plot_resolutions)
