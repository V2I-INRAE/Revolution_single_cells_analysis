suppressPackageStartupMessages(library(Seurat))

read_clustering_input <- function(input_file, method, dims) {
  obj <- readRDS(input_file)
  reduction <- switch(method,
    unintegrated = "pca", harmony = "harmony", cca = "integrated_cca")
  cells <- colnames(obj)
  metadata <- obj[[]]
  stopifnot(
    "Cell identifiers and metadata must align" =
      !anyDuplicated(cells) && identical(cells, rownames(metadata)),
    "Sample and condition metadata must be complete" =
      all(c("sample", "pressure", "time_point") %in% names(metadata)) &&
      !anyNA(metadata[, c("sample", "pressure", "time_point")]),
    "Requested reduction is missing" = reduction %in% Reductions(obj)
  )
  embedding <- Embeddings(obj, reduction)
  stopifnot("Selected components must exist, be finite and align with cells" =
    max(dims) <= ncol(embedding) && identical(rownames(embedding), cells) &&
      all(is.finite(embedding[, dims, drop = FALSE])))
  obj
}

clustering_provenance <- function(obj, input_file) {
  # Canonical ordering makes comparison independent of checkpoint cell order.
  cells <- sort(colnames(obj), method = "radix")
  baseline <- list(
    metadata = obj[[]][cells, c("sample", "pressure", "time_point")],
    pca = Embeddings(obj, "pca")[cells, , drop = FALSE],
    loadings = Loadings(obj, "pca"), features = VariableFeatures(obj[["RNA"]])
  )
  list(input_file = normalizePath(input_file, mustWork = TRUE),
    baseline_hash = digest::digest(baseline, algo = "sha256"),
    integration = obj@misc$integration)
}

save_clustering_results <- function(obj, data_dir, output_dir) {
  obj@misc$clustering$session_info <- sessionInfo()
  assignments <- obj[[]][, obj@misc$clustering$partitions$column, drop = FALSE]
  write.csv(data.frame(cell = rownames(assignments), assignments),
    file.path(output_dir, "cluster_assignments.csv"), row.names = FALSE)
  saveRDS(obj, file.path(data_dir, "lognorm.rds"))
}
