# Standalone monitoring workflow for DoubletFinder and Scrublet diagnostics.
# Run from the analysis repository root. Optional arguments restrict the run:
#   Rscript R/doublets/monitor.R REVO26-C

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

source("R/qc/params.R")
source("R/qc/io.R")
source("R/qc/reduce.R")
source("R/qc/doublet_finder.R")
source("R/qc/plots.R")
source("R/doublets/runners.R")
source("R/doublets/collect.R")
source("R/doublets/plots.R")

output_dir <- file.path("results", "doublets")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sample_dirs <- list.dirs("data/raw_data", recursive = FALSE, full.names = FALSE)
samples <- sort(sample_dirs[file.exists(file.path(
  "data/raw_data",
  sample_dirs,
  paste0(sample_dirs, "_RSEC_MolsPerCell_MEX.zip")
))])
requested_samples <- commandArgs(trailingOnly = TRUE)
if (length(requested_samples) > 0L) {
  stopifnot("Unknown requested sample" = all(requested_samples %in% samples))
  samples <- requested_samples
}
stopifnot("No input samples found" = length(samples) > 0L)

cell_results <- vector("list", length(samples))
scrublet_scores <- vector("list", length(samples))
bcmvn_results <- vector("list", length(samples))
call_summaries <- vector("list", length(samples))
sensitivity_results <- list()

for (i in seq_along(samples)) {
  sample_id <- samples[[i]]
  cat(sprintf(
    "\n=== Doublet monitoring: %s (%d/%d) ===\n",
    sample_id, i, length(samples)
  ))

  mex_dir <- read_filtered_matrix("data/raw_data", sample_id)
  obj <- build_seurat_obj(mex_dir, sample_id)
  reduced_obj <- preprocess_for_doublets(obj)
  coordinates <- compute_diagnostic_umap(reduced_obj)

  doubletfinder <- run_doubletfinder_monitor(reduced_obj, sample_id)
  scrublet <- run_scrublet_monitor(obj, sample_id)

  cell_results[[i]] <- collect_doublet_cells(
    obj,
    sample_id,
    coordinates,
    doubletfinder,
    scrublet
  )
  scrublet_scores[[i]] <- scrublet$scores
  bcmvn_results[[i]] <- doubletfinder$bcmvn
  call_summaries[[i]] <- collect_call_summary(
    doubletfinder$settings,
    scrublet$settings
  )
  if (!is.null(doubletfinder$sensitivity)) {
    sensitivity_results[[length(sensitivity_results) + 1L]] <- (
      doubletfinder$sensitivity
    )
  }
  rm(obj, reduced_obj, doubletfinder, scrublet)
  gc()
}

cell_data <- do.call(rbind, cell_results)
score_data <- do.call(rbind, scrublet_scores)
bcmvn_data <- do.call(rbind, bcmvn_results)
call_summary <- do.call(rbind, call_summaries)
nfeature_summary <- collect_nfeature_summary(cell_data)
sensitivity_summary <- if (length(sensitivity_results) == 0L) {
  data.frame(
    sample_id = character(),
    selected_pK = numeric(),
    adjacent_pK = numeric(),
    selected_calls = integer(),
    adjacent_calls = integer(),
    changed_calls = integer(),
    jaccard = numeric()
  )
} else {
  do.call(rbind, sensitivity_results)
}

write.csv(
  call_summary,
  file.path(output_dir, "call_summary.csv"),
  row.names = FALSE
)
write.csv(
  nfeature_summary,
  file.path(output_dir, "nfeature_summary.csv"),
  row.names = FALSE
)
write.csv(
  sensitivity_summary,
  file.path(output_dir, "pk_boundary_sensitivity.csv"),
  row.names = FALSE
)
saveRDS(
  list(
    cells = cell_data,
    scrublet_scores = score_data,
    bcmvn = bcmvn_data,
    call_summary = call_summary,
    nfeature_summary = nfeature_summary,
    pk_boundary_sensitivity = sensitivity_summary,
    samples = samples,
    generated_at = Sys.time()
  ),
  file.path(output_dir, "monitor_diagnostics.rds")
)

save_qc_plot(
  plot_scrublet_score_histograms(score_data, call_summary),
  "scrublet_score_histograms",
  output_dir,
  4000,
  3400
)
save_qc_plot(
  plot_pk_sweeps(bcmvn_data),
  "doubletfinder_pk_sweeps",
  output_dir,
  4000,
  3400
)
save_qc_plot(
  plot_continuous_score_umap(
    cell_data,
    "doubletfinder_score",
    "DoubletFinder pANN on per-sample UMAP coordinates",
    "pANN"
  ),
  "doubletfinder_score_umap",
  output_dir,
  4000,
  3400
)
save_qc_plot(
  plot_continuous_score_umap(
    cell_data,
    "scrublet_score",
    "Scrublet score on per-sample UMAP coordinates",
    "Scrublet score"
  ),
  "scrublet_score_umap",
  output_dir,
  4000,
  3400
)
save_qc_plot(
  plot_monitor_binary_umap(
    cell_data,
    "doubletfinder_class",
    "DoubletFinder calls"
  ),
  "doubletfinder_calls_umap",
  output_dir,
  4000,
  3400
)
save_qc_plot(
  plot_monitor_binary_umap(
    cell_data,
    "scrublet_current_class",
    "Scrublet calls using the current score > 0.15 rule"
  ),
  "scrublet_current_calls_umap",
  output_dir,
  4000,
  3400
)
if (any(!is.na(cell_data$scrublet_automatic_class))) {
  save_qc_plot(
    plot_monitor_binary_umap(
      cell_data,
      "scrublet_automatic_class",
      "Scrublet automatic-threshold calls"
    ),
    "scrublet_automatic_calls_umap",
    output_dir,
    4000,
    3400
  )
}
save_qc_plot(
  plot_nfeature_summary(nfeature_summary),
  "nfeature_summary",
  output_dir,
  4200,
  3600
)

cat("\nDoublet monitoring completed for", length(samples), "sample(s).\n")
print(call_summary)
