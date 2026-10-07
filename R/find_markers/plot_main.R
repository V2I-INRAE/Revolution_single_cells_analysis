args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript R/find_markers/plot_main.R <clustering/analysis.rds> <marker-results-directory> <resolutions: 0.1,0.2,...>")
}
source("R/find_markers/io.R")
source("R/find_markers/plots.R")

output_dir <- normalizePath(args[2], mustWork = TRUE)
resolutions <- as.numeric(strsplit(args[3], ",", fixed = TRUE)[[1]])
if (length(list.files(output_dir, pattern = "\\.(png|svg)$", recursive = TRUE)) > 0L) {
  stop("Plot outputs already exist; preserve or remove them explicitly before replotting")
}
obj <- read_marker_input(args[1])

for (resolution in resolutions) {
  directory <- file.path(output_dir, paste0("res_", format(resolution, nsmall = 1)))
  top25 <- read.csv(file.path(directory, "markers_top25.csv"), colClasses = c(gene = "character"))
  top5 <- read.csv(file.path(directory, "markers_top5.csv"), colClasses = c(gene = "character"))
  counts <- read.csv(file.path(directory, "cluster_summary.csv"))
  plot_marker_bars(top25, counts$cluster, resolution, directory)
  plot_markers_dotplot(obj, top5, resolution, directory)
  plot_top_markers_heatmap(obj, top5, resolution, directory)
  gc()
}
