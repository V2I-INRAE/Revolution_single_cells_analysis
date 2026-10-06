# Plot saved, compatible runs without recomputing normalization or integration.
args <- commandArgs(trailingOnly = TRUE)
embedding <- "both"
if (length(args) > 0L && startsWith(tail(args, 1L), "--embedding=")) {
  embedding <- match.arg(sub("^--embedding=", "", tail(args, 1L)),
    c("umap", "tsne", "both"))
  args <- head(args, -1L)
}
if (length(args) != 4L) {
  stop("Usage: Rscript R/integration/plot_main.R <scvi.rds> <harmony.rds> <cca.rds> <new_output_dir> [--embedding=umap|tsne|both]")
}
embeddings <- if (embedding == "both") c("umap", "tsne") else embedding
output_dir <- args[[4]]
if (dir.exists(output_dir)) stop("Output directory already exists: ", output_dir)
source("R/integration/io.R")
source("R/integration/plots.R")

obj <- read_integration_comparison(args[[1]], args[[2]], args[[3]], embeddings = embeddings)
dir.create(dirname(output_dir), recursive = TRUE, showWarnings = FALSE)
if (!dir.create(output_dir, showWarnings = FALSE)) {
  stop("Cannot create new output directory: ", output_dir)
}
write.csv(obj@misc$integration_comparison,
  file.path(output_dir, "comparison_sources.csv"), row.names = FALSE)
for (embedding in embeddings) {
  for (group_by in c("sample", "pressure", "time", "pressure_time")) {
    plot_integration_umaps(obj, group_by = group_by, output_dir = output_dir,
      embedding = embedding)
  }
  plot_integration_umaps(obj, group_by = "sample", facet_samples = TRUE,
    output_dir = output_dir, embedding = embedding)
}
