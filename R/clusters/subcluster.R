# Refit one saved coarse identity using the existing pipeline functions.
# Usage: Rscript R/clusters/subcluster.R <input.rds> "<subset label>" <folder>
#   [--only-pca | --hvg=<number> --dims=<number> --k-param=<number> --resolutions=0.1,0.2,...]
# Example: data/job-44542861/clean_annotated_object.rds DC DC --hvg=1500 --dims=10 --k-param=50
# Full-run overrides are optional; --only-pca accepts no overrides.
# No normalization, QC filtering or cell-cycle scoring is repeated.

suppressPackageStartupMessages(library(Seurat))

source("R/utils/cli.R", local = TRUE)
source("R/config.R", local = TRUE)
source("R/io.R", local = TRUE) # save_final_object(), format_value()
source("R/lognorm.R", local = TRUE)
source("R/reduce.R", local = TRUE)
source("R/integrate.R", local = TRUE)
source("R/cluster.R", local = TRUE)
source("R/utils/save_plots.R", local = TRUE)
source("R/plots/plots_normalized_data.R", local = TRUE)
source("R/plots/plots_clustered_data.R", local = TRUE)

run_subclustering <- function(input, subset_label, folder, only_pca = FALSE,
  n_hvg = NULL, dims = NULL, k_param = NULL, resolutions = NULL) {
  validate_subcluster_options(only_pca, n_hvg, dims, k_param, resolutions)
  # Print each warning where it occurs instead of an end-of-run count.
  warning_options <- options(warn = 1)
  on.exit(options(warning_options), add = TRUE)

  # Keep overrides local; each call starts from the loaded configuration.
  config <- config
  if (!is.null(n_hvg)) config$norm$n_variable_features <- as.integer(n_hvg)
  if (!is.null(dims)) config$dims <- seq_len(dims)
  if (!is.null(k_param)) config$clustering$k_param <- as.integer(k_param)
  if (!is.null(resolutions)) config$clustering$resolutions <- resolutions
  if (max(config$dims) > config$norm$n_pcs) {
    stop("--dims cannot exceed the number of computed PCs (", config$norm$n_pcs, ").")
  }
  validate_config(config)
  job_id <- Sys.getenv("SLURM_JOB_ID")
  if (!nzchar(job_id)) {
    job_id <- paste0("manual-", format(Sys.time(), "%Y%m%d-%H%M%S"), "-", Sys.getpid())
  }
  run_id <- paste0("job-", job_id)
  # The parent job is unchanged when input is an already saved subset.
  parent_job <- regmatches(input, regexpr("job-[0-9]+", input))
  data_dir <- file.path("data", parent_job, "subsets", folder, run_id)
  plot_dir <- file.path("results", parent_job, "clusters_plots", "subset", folder, run_id)
  if (dir.exists(plot_dir) || dir.exists(data_dir)) {
    stop("This job's output directory already exists; refusing to overwrite it.")
  }
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  info_path <- file.path(plot_dir, ".INFO")
  info <- c(
    "# Subclustering run",
    "Status: running",
    paste0("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("Job: ", run_id),
    paste0("Input (read-only): ", input),
    paste0("Input record: ", tools::file_path_sans_ext(input), ".INFO"),
    paste0("Pipeline parent: ", parent_job, "; subset: ", subset_label),
    paste0("Mode: ", if (only_pca) "PCA only" else "full subclustering"),
    "Methods: reuse normalized RNA; refit HVGs, scaling/regression and PCA; no QC or cell-cycle rescoring.",
    paste0("HVGs requested: ", config$norm$n_variable_features,
      " (", if (is.null(n_hvg)) "configured default" else "--hvg override", ")"),
    sprintf("PCA: compute %d components; downstream dimensions %s%s",
      config$norm$n_pcs, paste(range(config$dims), collapse = ":"),
      if (only_pca) " (not run)" else if (is.null(dims)) " (configured default)" else " (--dims override)"),
    paste0("Seed: ", config$run$seed),
    paste0("Output: ", file.path(data_dir, "subclustered.rds")),
    paste0("Figures: ", plot_dir)
  )
  completed <- FALSE
  run_results <- "Validation: pending"
  writeLines(info, info_path)
  on.exit({
    info[info == "Status: running"] <- paste0("Status: ",
      if (completed) "completed" else "failed/incomplete; see execution log")
    writeLines(c(info, paste0("Finished: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
      run_results), info_path)
  }, add = TRUE)

  info <- c(info, paste0("Normalization settings: ", format_value(config$norm)),
    if (!only_pca) c(
      paste0("Integration settings: ", format_value(config$integration)),
      paste0("Clustering settings: ", format_value(config$clustering))
    ))
  writeLines(info, info_path)
  future::plan("sequential")

  input_info <- file.info(input)[, c("size", "mtime")]
  message("Loading (read-only): ", input)

  obj <- readRDS(input)

  selected_cells <- colnames(obj)[as.character(obj$subset) == subset_label]
  obj <- subset(obj, cells = selected_cells)
  DefaultAssay(obj) <- "RNA"
  gc()

  # Keep parent cluster labels outside the prefix used by ClustreePlot.
  parent_columns <- c(
    grep(
      "_snn_res\\.",
      names(obj[[]]),
      value = TRUE
    ),
    "seurat_clusters"
  )

  for (column in parent_columns) {
    obj[[paste0("parent_", column)]] <- obj[[column, drop = TRUE]]
    obj[[column]] <- NULL
  }

  obj@misc$subcluster <- list(
    input = input, subset_label = subset_label, run_id = run_id,
    only_pca = only_pca,
    overrides = list(hvg = n_hvg, dims = dims, k_param = k_param, resolutions = resolutions),
    settings = list(norm = config$norm, dims = config$dims, seed = config$run$seed,
      integration = if (!only_pca) config$integration else NULL,
      clustering = if (!only_pca) config$clustering else NULL),
    parent_analysis = obj@misc[c("integration", "clustering")]
  )

  for (name in c(Reductions(obj), Graphs(obj), Neighbors(obj))) obj[[name]] <- NULL
  rna <- obj[["RNA"]]

  for (layer in Layers(rna, search = "^scale\\.data($|\\.)")) {
    LayerData(rna, layer = layer) <- NULL
  }

  obj[["RNA"]] <- rna
  rm(rna)

  set.seed(config$run$seed)
  obj <- find_hvgs(obj, config$norm)
  obj <- scale_data(obj, config$norm)
  obj <- run_pca_analysis(obj, config$norm$n_pcs)
  run_results <- c(
    run_results,
    paste0(
      "Actual: ",
      ncol(obj),
      " cells; ",
      length(VariableFeatures(obj)), 
      " HVGs; ", 
      ncol(Embeddings(obj, "pca")), 
      " PCA components.")
  )

  if (!only_pca && max(config$dims) > ncol(Embeddings(obj, "pca"))) {
    stop("Requested downstream dimensions exceed the actual number of PCA components.")
  }

  plot_variable_features(obj, plot_dir)
  plot_pca_elbow(obj, plot_dir)
  plot_pca_loadings(obj, plot_dir, dims = 1:10)

  if (!only_pca) {
    obj <- integrate_object(obj, config$integration, config$dims, config$run$seed)
    obj <- run_integration_umaps(
      obj,
      config$dims,
      config$run$seed,
      integration_reductions(config$integration$method),
      config$integration$umap
    )
    obj <- cluster_object(obj, config$integration$method, config$clustering, config$dims)
    obj <- validate_clustering_umap(obj, "umap_harmony", "harmony", config$dims)
  }

  # Both modes save and verify the subset before returning or cluster plotting.
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  final_path <- file.path(data_dir, "subclustered.rds")
  save_final_object(obj, final_path)

  reopened <- readRDS(final_path)

  stopifnot(
    identical(colnames(reopened), selected_cells),
    identical(reopened[[]], obj[[]]), identical(Idents(reopened), Idents(obj)),
    identical(reopened@assays, obj@assays),
    identical(reopened@reductions, obj@reductions),
    identical(reopened@graphs, obj@graphs),
    identical(reopened@misc$subcluster, obj@misc$subcluster))
  rm(reopened)
  gc()

  stopifnot(identical(file.info(input)[, c("size", "mtime")], input_info))
  run_results[run_results == "Validation: pending"] <- paste(
    "Validation: saved RDS reopened; cell IDs, metadata, identities, assays/layers/counts,",
    "reductions, graphs and subclustering settings match; input size/mtime unchanged."
  )
  if (only_pca) {
    completed <- TRUE
    message("PCA-only subset saved: ", final_path, "\nFigures: ", plot_dir,
      "\nSaved RDS verified; annotated input unchanged. Stopped before integration and clustering.")
    return(invisible(final_path))
  }

  plot_cluster_resolutions(obj, plot_dir)
  completed <- TRUE

  message("Subclustering completed: ", final_path, "\nFigures: ", plot_dir,
    "\nSaved RDS verified; annotated input unchanged. Figure review left to owner.")
  invisible(final_path)
}

# Sourcing this file defines the runner without starting a run.
if (sys.nframe() == 0L) {
  args <- parse_subcluster_args(commandArgs(trailingOnly = TRUE))
  # Seurat derives command names from the call stack; keep a direct entry call.
  run_subclustering(args$input, args$subset_label, args$folder,
    only_pca = args$only_pca, n_hvg = args$n_hvg, dims = args$dims,
    k_param = args$k_param, resolutions = args$resolutions)
}
