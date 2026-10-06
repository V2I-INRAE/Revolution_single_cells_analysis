# Add embeddings to completed checkpoints without integration or clustering.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L || !args[1] %in% c("integration", "clustering")) {
  stop("Usage: Rscript R/integration/add_tsne.R <integration|clustering> <lognorm.rds|analysis.rds> <new_output_dir>")
}
source("R/integration/reduce.R")
source("R/clustering/cluster.R")

stage <- args[1]
input_file <- normalizePath(args[2], mustWork = TRUE)
output_dir <- args[3]
if (dir.exists(output_dir)) stop("Output directory already exists: ", output_dir)

if (stage == "integration") {
  obj <- readRDS(input_file)
  config <- obj@misc$integration
  stopifnot("Expected a completed integration checkpoint" =
    config$method %in% c("scvi", "harmony", "cca") && config$route == "lognorm",
    "Keep the original run directory basename for integration provenance" =
      basename(output_dir) == config$run_id)
  obj <- run_integration_tsnes(obj, config)
} else {
  run <- readRDS(input_file)
  obj <- readRDS(run$checkpoint)
  stopifnot("Clustering bundle and checkpoint must agree" =
    identical(run$method, obj@misc$clustering$method) &&
      identical(run$settings, obj@misc$clustering$settings) &&
      identical(run$assignments,
        obj[[]][, run$partitions$column, drop = FALSE]))
  obj <- clustering_tsne(obj)
}
obj@misc$tsne_source_file <- input_file
dir.create(dirname(output_dir), recursive = TRUE, showWarnings = FALSE)
if (!dir.create(output_dir, showWarnings = FALSE)) {
  stop("Cannot create new output directory: ", output_dir)
}
checkpoint <- file.path(output_dir, "lognorm.rds")
saveRDS(obj, checkpoint)
if (stage == "clustering") {
  run$checkpoint <- normalizePath(checkpoint)
  run$tsne_source_file <- input_file
  saveRDS(run, file.path(output_dir, "analysis.rds"))
}
message("Saved t-SNE-augmented checkpoint: ", checkpoint)
