# Add extra Leiden resolutions to a saved clustered object (no pipeline rerun).
# Usage: Rscript R/clusters/add_resolution.R <data/job-XXX/clustered.rds> <res> [<res> ...]
# Example: Rscript R/clusters/add_resolution.R data/job-44493014/clustered.rds 0.8 0.9 1.2
#
# Runs FindClusters for the requested resolutions on the stored neighbor graph
# (settings from misc$clustering), overwrites the object, and updates both
# .INFO copies. Assignment CSVs and figures are NOT regenerated and go stale.

suppressPackageStartupMessages({
  library(Seurat)
})

source("R/io.R") # for save_final_object()

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: Rscript R/clusters/add_resolution.R <data/job-XXX/clustered.rds> <res> [<res> ...]")
}

input <- args[[1L]]

new_res <- suppressWarnings(as.numeric(args[-1L]))
if (any(!is.finite(new_res) | new_res <= 0)) {
  stop("Resolutions must be finite positive numbers.")
}
new_res <- sort(new_res)

message("Loading saved object: ", input)
obj <- readRDS(input)
cfg <- obj@misc$clustering
settings <- cfg$settings

todo <- setdiff(new_res, cfg$partitions$resolution)
if (length(todo) == 0L) {
  message("All requested resolutions are already present; nothing to do.")
  quit(status = 0L)
}

message("FindClusters on '", cfg$graph, "' for: ", paste(todo, collapse = ", "))
obj <- FindClusters(
  obj,
  graph.name = cfg$graph,
  resolution = todo,
  algorithm = settings$algorithm,
  leiden_method = settings$leiden_method,
  leiden_objective_function = settings$leiden_objective_function,
  random.seed = settings$seed,
  n.iter = settings$n_iter
)

partitions <- rbind(
  cfg$partitions,
  data.frame(resolution = todo, column = paste0(cfg$graph, "_res.", todo))
)
partitions <- partitions[order(partitions$resolution), , drop = FALSE]
rownames(partitions) <- NULL
settings$resolutions <- partitions$resolution
obj@misc$clustering$settings <- settings
obj@misc$clustering$partitions <- partitions

# Active identity is the last (highest) resolution, as in the pipeline.
active_column <- tail(partitions$column, 1L)
Idents(obj) <- active_column
obj$seurat_clusters <- Idents(obj)

counts <- vapply(partitions$column,
  function(column) length(unique(obj[[]][[column]])), integer(1))
message("Clusters: ", paste(sprintf("%s at %s", counts, partitions$resolution), collapse = "; "))

save_final_object(obj, input)
message("Object overwritten: ", input)

job_id <- basename(dirname(input))
info_paths <- c(file.path("data", job_id, ".INFO"), file.path("results", job_id, ".INFO"))
lines <- readLines(info_paths[[1L]])
lines[[grep("^- resolutions: ", lines)]] <-
  paste0("- resolutions: ", paste(partitions$resolution, collapse = ", "))
lines[[grep("^- Leiden clusters: ", lines)]] <-
  paste0("- Leiden clusters: ",
    paste(sprintf("%s at resolution %s", counts, partitions$resolution), collapse = "; "))
lines <- c(lines, paste0("- Extra resolutions ", paste(todo, collapse = ", "),
  " added ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  " via R/clusters/add_resolution.R; active identity: ", active_column, ". ",
  "cluster_assignments.csv and cluster figures still cover the old resolutions only."))
for (path in info_paths) writeLines(lines, path)
message("Updated both .INFO copies.")
