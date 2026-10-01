# Plot saved clustering runs only. The joined object is never persisted.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: Rscript R/clustering/markers_main.R <new_output_dir> <analysis.rds> [<analysis.rds> ...]")
}
source("R/clustering/params.R")
source("R/clustering/io.R")
source("R/clustering/plots.R")
source("R/clustering/marker_plots.R")
runs <- read_clustering_comparison(args[-1])
output_dir <- args[1]
dir.create(dirname(output_dir), recursive = TRUE, showWarnings = FALSE)
stopifnot("Choose a new marker output directory" = dir.create(output_dir, showWarnings = FALSE))
write.csv(data.frame(method = names(runs), source = normalizePath(args[-1]),
  resolution = marker_params$resolution),
  file.path(output_dir, "marker_sources.csv"), row.names = FALSE)

for (method in names(runs)) {
  run <- runs[[method]]
  message("Reading ", method, " checkpoint: ", run$checkpoint)
  obj <- readRDS(run$checkpoint)
  column <- run$partitions$column[run$partitions$resolution == marker_params$resolution]
  stopifnot("Checkpoint must match completed run and selected assignments" =
    identical(obj@misc$clustering$method, method) && length(column) == 1L &&
      identical(obj[[]][, column, drop = FALSE], run$assignments[, column, drop = FALSE]))
  prepared <- prepare_marker_umaps(obj, marker_params$panels, marker_params$resolution)
  message("Plotting ", ncol(obj), " cells; reference clusters: ",
    length(unique(obj[[]][, column])), "; resolution: ", marker_params$resolution)
  rm(obj)
  method_dir <- file.path(output_dir, method)
  dir.create(method_dir)
  plot_marker_umaps(prepared, marker_params$panels, method_dir)
  rm(prepared)
  gc()
}
