# Single-entry pipeline: raw counts -> QC -> normalization -> integration -> clustering.
# Launched by scripts/sbatch_pipeline.sh; everything is configured in R/config.R.
# The final clustered object and all figures are produced in one session with
# no intermediate checkpoints; .INFO (data/ and results/) tracks each stage.

suppressPackageStartupMessages({
  library(Seurat)
})

source("R/config.R")
source("R/io.R")
source("R/filter.R")
source("R/scrublet.R")
source("R/reduce.R")
source("R/lognorm.R")
source("R/cell_cycle.R")
source("R/integrate.R")
source("R/cluster.R")
source("R/utils/save_plots.R")
source("R/plots/plots_qc_data.R")
source("R/plots/plots_normalized_data.R")
source("R/plots/plots_integrated_data.R")
source("R/plots/plots_clustered_data.R")

future::plan("sequential")
validate_config(config)

job_id <- Sys.getenv("SLURM_JOB_ID")
run_id <- if (nzchar(job_id)) paste0("job-", job_id) else
  paste(format(Sys.time(), "%Y%m%d-%H%M%S"), Sys.getpid(), sep = "-")

input_paths <- setNames(
  vapply(config$run$samples, qc_input_path, character(1),
    input_source = config$run$input_source),
  config$run$samples
)
missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs) > 0L) {
  stop("Missing QC inputs:\n", paste(missing_inputs, collapse = "\n"))
}

layout <- create_run_layout(run_id)
record <- new_run_record(config, input_paths)
record$run_id <- run_id
record$final_object <- file.path(layout$data_dir, "clustered.rds")
final_object_path <- record$final_object
write_info(record, layout$info_paths)

message("Pipeline run: ", run_id,
  "\nInput source: ", config$run$input_source,
  "\nIntegration method: ", config$integration$method,
  "\nData: ", layout$data_dir,
  "\nResults: ", layout$results_dir)

# Stage runner: on error, stamp the failing stage into both .INFO copies and
# exit nonzero. OOM/SIGKILL cannot run this; such a record stays "running" and
# is reconciled against Slurm afterwards.
run_stage <- function(stage, code) {
  tryCatch(
    force(code),
    error = function(e) {
      record <<- fail_run_record(record, stage, e)
      write_info(record, layout$info_paths)
      message("Pipeline failed at stage: ", stage, "\n", conditionMessage(e))
      quit(status = 1L)
    }
  )
}

# --- QC -----------------------------------------------------------------------
# All input cells are kept until the before/after figures are saved; the
# labeled object is then dropped and only the clean object carries forward.
qc <- run_stage("QC", local({
  labeled_objs <- list()
  for (sample in config$run$samples) {
    counts <- read_qc_counts(sample, config$run$input_source)
    obj <- build_seurat_obj_from_counts(counts, sample, config$qc)
    rm(counts)
    gc()

    obj <- inspect_seurat_qc(obj, config$qc)
    obj <- label_qc_outliers(obj, config$qc$filtering)

    reduced_obj <- preprocess_qc_umap(obj, config$qc$doublets)
    coordinates <- compute_diagnostic_umap(reduced_obj, config$qc$doublets)
    obj <- AddMetaData(obj, coordinates[colnames(obj), , drop = FALSE])

    obj <- detect_scrublet_doublets(obj, config$qc$doublets)
    labeled_objs[[sample]] <- label_retained_cells(obj)
    rm(obj, reduced_obj, coordinates)
    gc()
  }

  labeled <- concatenate_qc_objects(labeled_objs, config$qc)
  rm(labeled_objs)
  gc()
  labeled@misc$qc$run_id <- run_id
  labeled@misc$qc$input_source <- config$run$input_source
  labeled@misc$qc$input_paths <- input_paths
  labeled@misc$qc$threshold_reference <- "per_sample_median_mad_capped"
  labeled@misc$qc$feature_min_inclusive <- TRUE

  # Before/after figures compare all input cells with QC-passing cells
  # (doublets not yet removed), on the original per-cell measurements.
  plot_qc_comparison(labeled, layout$plot_dirs[["qc_plots"]])
  plot_qc_umaps(labeled, layout$plot_dirs[["qc_plots"]])

  n_input <- ncol(labeled)
  n_qc_pass <- sum(labeled$keep_qc)
  clean <- filter_labeled_cells(labeled, config$qc)
  clean <- remove_listed_genes(clean, config$qc)
  validate_clean_object(clean)
  rm(labeled)
  gc()

  list(
    obj = clean,
    lines = c(
      paste0(n_input, " input cells; ", n_qc_pass, " QC-passing; ",
        ncol(clean), " Scrublet-clean (", sprintf("%.2f%%", 100 * ncol(clean) / n_input), ")"),
      paste0(nrow(clean), " genes after the per-feature minimum; ",
        length(clean@misc$qc$removed_genes), " listed genes removed; ",
        length(config$run$samples), " samples"),
      "8 QC figures saved (6 before/after comparisons, 2 per-sample UMAP sets)"
    )
  )
}))
record <- record_stage(record, "QC", qc$lines)
write_info(record, layout$info_paths)

# --- Normalization --------------------------------------------------------------
norm <- run_stage("normalization", local({
  set.seed(config$run$seed)
  obj <- run_lognormalize(qc$obj, config$norm)
  obj <- score_cell_cycle(obj, config$cell_cycle)
  plot_cc_tail_removal(obj, layout$plot_dirs[["norm_plots"]])
  obj <- remove_cc_tail_cells(obj)
  obj <- find_hvgs(obj, config$norm)
  obj <- scale_data(obj, config$norm)
  obj <- run_pca_analysis(obj, config$norm$n_pcs)
  obj@misc$norm_feat <- list(
    run_id = run_id, params = config$norm, pipeline = "single-entry"
  )

  plot_variable_features(obj, layout$plot_dirs[["norm_plots"]])
  plot_hvg_overlap(obj, layout$plot_dirs[["norm_plots"]])
  plot_pca_elbow(obj, layout$plot_dirs[["norm_plots"]])
  plot_pca_loadings(obj, layout$plot_dirs[["norm_plots"]])

  n_figures <- length(config$run$samples) + 8L
  list(
    obj = obj,
    lines = c(
      paste0(ncol(obj), " cells; ", nrow(obj), " genes; ",
        length(VariableFeatures(obj)), " consensus HVGs; ",
        config$norm$n_pcs, " PCs"),
      paste0("regressed: ", paste(config$norm$vars_to_regress, collapse = ", ")),
      paste0("cell-cycle: ", obj@misc$cell_cycle$removal$removed,
        " tail cells removed (",
        sprintf("%.2f%%", obj@misc$cell_cycle$removal$pct_removed), "); ",
        ncol(obj), " cells retained; cutoffs in misc$cell_cycle"),
      paste0(n_figures, " normalization figures saved (per-sample HVGs, cell-cycle removal, overlap, elbow, PC1-5 loadings)")
    )
  )
}))
record <- record_stage(record, "normalization", norm$lines)
write_info(record, layout$info_paths)
qc <- NULL
gc()

# --- Integration ------------------------------------------------------------------
int <- run_stage("integration", local({
  obj <- integrate_object(norm$obj, config$integration, config$dims, config$run$seed)
  obj <- run_integration_umaps(
    obj, config$dims, config$run$seed,
    integration_reductions(config$integration$method), config$integration$umap
  )
  figures <- plot_integration_figures(
    obj, layout$plot_dirs[["integration_plots"]], seed = config$run$seed
  )

  lines <- c(
    paste0("method: ", config$integration$method,
      "; dims ", min(config$dims), "-", max(config$dims),
      "; batch variable: sample"),
    paste0("PCA and ", config$integration$method, " UMAPs: min.dist ",
      config$integration$umap$min.dist, "; n.neighbors ",
      config$integration$umap$n.neighbors, "; seed ", config$run$seed),
    if (!is.null(obj@misc$integration$harmony_convergence)) {
      obj@misc$integration$harmony_convergence
    },
    paste0(length(figures), " UMAP figures saved")
  )
  list(obj = obj, lines = lines)
}))
record <- record_stage(record, "integration", int$lines)
write_info(record, layout$info_paths)
norm <- NULL
gc()

# --- Clustering ---------------------------------------------------------------------
clu <- run_stage("clustering", local({
  method <- config$integration$method
  obj <- cluster_object(int$obj, method, config$clustering, config$dims)
  obj <- validate_clustering_umap(
    obj,
    switch(method, harmony = "umap_harmony", cca = "umap_cca"),
    obj@misc$clustering$reduction,
    config$dims
  )

  counts <- cluster_counts(obj)
  obj@misc$run <- list(
    run_id = run_id,
    pipeline = "single-entry QC -> normalization -> integration -> clustering",
    started = record$started,
    config = config
  )
  obj@misc$clustering$session_info <- sessionInfo()

  # The final object is the run's deliverable, not an intermediate: save it
  # atomically BEFORE rendering figures, so a plotting failure cannot lose it.
  save_final_object(obj, final_object_path)
  write_cluster_assignments(
    obj, file.path(layout$plot_dirs[["clusters_plots"]], "cluster_assignments.csv")
  )
  plot_cluster_resolutions(obj, layout$plot_dirs[["clusters_plots"]])

  list(
    obj = obj,
    lines = c(
      paste0("Leiden clusters: ",
        paste(sprintf("%s at resolution %s", counts, names(counts)), collapse = "; ")),
      paste0("final object saved: ", final_object_path),
      "cluster figure + assignments CSV saved (umap_resolutions, clustree)"
    )
  )
}))
record <- record_stage(record, "clustering", clu$lines)
int <- NULL
gc()
record <- finish_run_record(record)
write_info(record, layout$info_paths)

message("Pipeline completed: ", run_id,
  "\nFinal object: ", final_object_path,
  "\nFigures: ", layout$results_dir)
