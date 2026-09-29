# Scrublet scores and calls; method-specific filtering happens in main.R.
detect_scrublet_doublets <- function(seurat_obj) {
  result <- scrubletR::scrublet_R(
    seurat_obj,
    python_home = Sys.getenv(
      "RETICULATE_PYTHON",
      unset = file.path(getwd(), ".venv-scrublet", "bin", "python")
    ),
    expected_doublet_rate = bd_multiplet_rate(ncol(seurat_obj)),
    min_counts = 3L,
    n_prin_comps = 30L,
    return_results_only = TRUE,
    threshold = NULL
  )
  seurat_obj$doublet_scores <- result$doublet_scores
  seurat_obj$predicted_doublets <- result$doublet_scores > 0.15
  return(seurat_obj)
}
