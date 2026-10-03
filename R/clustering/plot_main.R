# Explicit completed runs only; no normalization, clustering or UMAP fitting.
args <- commandArgs(trailingOnly = TRUE)
plot_resolutions <- NULL
if (length(args) > 0L && startsWith(tail(args, 1L), "--plot-resolutions=")) {
  plot_resolutions <- as.numeric(strsplit(sub("^--plot-resolutions=", "", tail(args, 1L)), ",")[[1]])
  args <- head(args, -1L)
}
if (length(args) < 2L) {
  stop("Usage: Rscript R/clustering/plot_main.R <new_output_dir> <analysis.rds> [<analysis.rds> ...] [--plot-resolutions=0.4,0.6,0.8]")
}

source("R/clustering/io.R")
source("R/clustering/diagnostics.R")
source("R/clustering/plots.R")

runs <- read_clustering_comparison(args[-1])
if (is.null(plot_resolutions)) {
  plot_resolutions <- runs[[1]]$partitions$resolution
}
stopifnot("Plot resolutions must be unique saved resolutions" =
  !anyNA(plot_resolutions) && !anyDuplicated(plot_resolutions) &&
    length(plot_resolutions) > 0L &&
    all(plot_resolutions %in% runs[[1]]$partitions$resolution))

output_dir <- args[1]
dir.create(dirname(output_dir), recursive = TRUE, showWarnings = FALSE)

if (!dir.create(output_dir, showWarnings = FALSE)) {
  stop("Comparison output directory already exists or is not writable: ", output_dir)
}

write.csv(data.frame(method = names(runs), source = normalizePath(args[-1])),
  file.path(output_dir, "comparison_sources.csv"), row.names = FALSE)
write.csv(data.frame(resolution = plot_resolutions),
  file.path(output_dir, "plotted_resolutions.csv"), row.names = FALSE)
comparison <- compare_clustering(runs)

for (name in names(comparison)) {
  write.csv(comparison[[name]], file.path(output_dir, paste0(name, ".csv")), row.names = FALSE)
}

summary <- do.call(rbind, lapply(runs, function(run) {
  data.frame(method = run$method, run$diagnostics$summary)
}))

write.csv(summary, file.path(output_dir, "metrics_by_resolution.csv"), row.names = FALSE)

plot_clustering_comparison(runs, comparison$agreement, output_dir, plot_resolutions)

for (i in seq_along(runs)) {
  # Only one full checkpoint in memory at a time. Figures live beside its tables.
  obj <- readRDS(runs[[i]]$checkpoint)
  plot_cluster_resolutions(obj, dirname(args[-1][i]), plot_resolutions)
  rm(obj)
  gc()
}
