# Scrublet scores and calls; filtering happens after the labeled checkpoint.

# Linear interpolation of the sample-specific BD multiplet percentages.
bd_multiplet_rate <- function(
  n_cells,
  table = qc_params$doublets$bd_multiplet_table
) {
  rate <- approx(
    table$cells,
    table$rate,
    xout = n_cells,
    rule = 2
  )$y
  if (n_cells > max(table$cells) || n_cells < min(table$cells)) {
    warning(
      "Cell count outside BD table range; rate clamped to ",
      min(table$rate), "-", max(table$rate), "%"
    )
  }
  return(rate / 100)
}

detect_scrublet_doublets <- function(seurat_obj) {
  rate <- bd_multiplet_rate(ncol(seurat_obj))
  result <- scrubletR::scrublet_R(
    seurat_obj,
    python_home = Sys.getenv(
      "RETICULATE_PYTHON",
      unset = file.path(getwd(), ".venv-scrublet", "bin", "python")
    ),
    expected_doublet_rate = rate,
    min_counts = 3L,
    n_prin_comps = 30L,
    return_results_only = TRUE,
    threshold = NULL
  )
  seurat_obj$doublet_scores <- result$doublet_scores
  seurat_obj$predicted_doublets <- result$doublet_scores > 0.15
  seurat_obj@misc$qc$scrublet <- list(
    expected_doublet_rate = rate, min_counts = 3L,
    n_prin_comps = 30L, threshold = 0.15, seed = 0L
  )
  return(seurat_obj)
}
