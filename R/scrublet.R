# Scrublet scores and calls; filtering happens after the labeled QC plots.

# Interpolate BD percentages; extend the last segment above the table without a cap.
bd_multiplet_rate <- function(n_cells, table) {
  rate <- approx(table$cells, table$rate, xout = n_cells, rule = 2)$y
  if (n_cells > max(table$cells)) {
    last <- tail(table, 2L)
    rate <- last$rate[2] + (n_cells - last$cells[2]) *
      diff(last$rate) / diff(last$cells)
  }
  return(rate / 100)
}

detect_scrublet_doublets <- function(seurat_obj, doublets_config) {
  scrublet_config <- doublets_config$scrublet
  rate <- bd_multiplet_rate(ncol(seurat_obj), doublets_config$bd_multiplet_table)
  message(sprintf(
    "Scrublet expected doublet rate: %.4f%% for %d input cells (BD interpolation; upper linear extrapolation, no cap).",
    100 * rate, ncol(seurat_obj)
  ))
  result <- scrubletR::scrublet_R(
    seurat_obj,
    python_home = Sys.getenv(
      "RETICULATE_PYTHON",
      unset = file.path(getwd(), ".venv-scrublet", "bin", "python")
    ),
    expected_doublet_rate = rate,
    min_counts = scrublet_config$min_counts,
    n_prin_comps = scrublet_config$n_prin_comps,
    return_results_only = TRUE,
    threshold = NULL
  )
  seurat_obj$doublet_scores <- result$doublet_scores
  seurat_obj$predicted_doublets <- result$doublet_scores > scrublet_config$threshold
  seurat_obj@misc$qc$scrublet <- c(
    expected_doublet_rate = rate,
    scrublet_config,
    rate_method = "BD linear interpolation; upper linear extrapolation without cap"
  )
  return(seurat_obj)
}
