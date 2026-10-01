suppressPackageStartupMessages(library(Seurat))

read_clustering_input <- function(input_file, method, dims) {
  obj <- readRDS(input_file)
  reduction <- switch(method,
    unintegrated = "pca", harmony = "harmony", scvi = "integrated_scvi")
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
      all(is.finite(embedding[, dims, drop = FALSE])),
    "scVI clustering must use its complete latent representation" =
      method != "scvi" || identical(dims, seq_len(ncol(embedding))))
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

save_clustering_results <- function(obj, diagnostics, provenance, data_dir, output_dir) {
  checkpoint <- file.path(data_dir, "lognorm.rds")
  saveRDS(obj, checkpoint)
  run <- list(
    method = obj@misc$clustering$method,
    settings = obj@misc$clustering$settings,
    graph_settings = obj@misc$clustering[c("algorithm", "nn_method", "distance",
      "n_trees", "n_start", "n_iter")],
    partitions = obj@misc$clustering$partitions,
    assignments = obj[[]][, obj@misc$clustering$partitions$column, drop = FALSE],
    provenance = provenance, diagnostics = diagnostics,
    checkpoint = normalizePath(checkpoint), session_info = sessionInfo()
  )
  write.csv(data.frame(cell = rownames(run$assignments), run$assignments),
    file.path(output_dir, "cluster_assignments.csv"), row.names = FALSE)
  for (name in names(diagnostics)) {
    write.csv(diagnostics[[name]], file.path(output_dir, paste0(name, ".csv")),
      row.names = FALSE)
  }
  # Written last: plotting/comparison consumes only completed run bundles.
  saveRDS(run, file.path(output_dir, "analysis.rds"))
}

read_clustering_comparison <- function(run_files) {
  runs <- lapply(run_files, readRDS)
  methods <- vapply(runs, `[[`, character(1), "method")
  stopifnot("Select one completed run per method" = !anyDuplicated(methods))
  names(runs) <- methods
  reference <- runs[[1]]
  for (run in runs) {
    stopifnot("Runs have incompatible cells, baseline data or analysis settings" =
      identical(run$provenance$baseline_hash, reference$provenance$baseline_hash) &&
      identical(run$settings, reference$settings) &&
      identical(run$graph_settings, reference$graph_settings) &&
      setequal(rownames(run$assignments), rownames(reference$assignments)) &&
      identical(run$diagnostics$diagnostic_cells, reference$diagnostics$diagnostic_cells))
  }
  runs
}
