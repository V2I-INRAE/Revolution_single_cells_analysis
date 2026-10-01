# Explicit completed runs only; no normalization, clustering or UMAP fitting.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: Rscript R/clustering/plot_main.R <new_output_dir> <analysis.rds> [<analysis.rds> ...]")
}

source("R/clustering/io.R")
source("R/clustering/diagnostics.R")
source("R/clustering/plots.R")

runs <- read_clustering_comparison(args[-1])

output_dir <- args[1]
dir.create(dirname(output_dir), recursive = TRUE, showWarnings = FALSE)

if (!dir.create(output_dir, showWarnings = FALSE)) {
  stop("Comparison output directory already exists or is not writable: ", output_dir)
}

write.csv(data.frame(method = names(runs), source = normalizePath(args[-1])),
  file.path(output_dir, "comparison_sources.csv"), row.names = FALSE)
comparison <- compare_clustering(runs)

for (name in names(comparison)) {
  write.csv(comparison[[name]], file.path(output_dir, paste0(name, ".csv")), row.names = FALSE)
}

summary <- do.call(rbind, lapply(runs, function(run) {
  data.frame(method = run$method, run$diagnostics$summary)
}))

write.csv(summary, file.path(output_dir, "metrics_by_resolution.csv"), row.names = FALSE)

plot_clustering_comparison(runs, comparison$agreement, output_dir)

for (i in seq_along(runs)) {
  # Only one full checkpoint in memory at a time. Figures live beside its tables.
  obj <- readRDS(runs[[i]]$checkpoint)
  plot_cluster_resolutions(obj, dirname(args[-1][i]))
  rm(obj)
  gc()
}
