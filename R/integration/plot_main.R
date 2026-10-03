# Plot saved, compatible runs without recomputing normalization or integration.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: Rscript R/integration/plot_main.R <scvi.rds> <harmony.rds> <cca.rds> <output_dir>")
}
source("R/integration/io.R")
source("R/integration/plots.R")

obj <- read_integration_comparison(args[[1]], args[[2]], args[[3]])
output_dir <- args[[4]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(obj@misc$integration_comparison,
  file.path(output_dir, "comparison_sources.csv"), row.names = FALSE)
for (group_by in c("sample", "pressure", "time", "pressure_time")) {
  plot_integration_umaps(obj, group_by = group_by, output_dir = output_dir)
}
plot_integration_umaps(obj, group_by = "sample", facet_samples = TRUE,
  output_dir = output_dir)
