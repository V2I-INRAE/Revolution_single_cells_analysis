# Self-contained refit of the cleaned global object.
# Input (read-only): data/job-44542861/clean_annotated_object.rds (234183 cells,
# 11 compartments; reductions/graphs/neighbors/scale layers absent by design).
# Output: data/job-44542861/refit_clean_annotated.rds + refit_clean_annotated.INFO.
# Method: reuse existing LogNormalize data; recompute HVGs, scaling/regression,
# PCA, Harmony, neighbours, Leiden and UMAPs. No QC, doublet detection or
# cell-cycle rescoring. Settings start from R/config.R with one approved
# departure: harmony max_iter=30 with early stopping (parent run used 10).
# PCA keeps the parent wrapper behaviour (Seurat internal seed 42);
# Harmony/UMAP/clustering use config seeds (1234). Sourcing defines the
# runner; the CLI takes no arguments.

suppressPackageStartupMessages(library(Seurat))

source("R/config.R", local = TRUE)
source("R/io.R", local = TRUE) # save_final_object(), format_value()
source("R/lognorm.R", local = TRUE)
source("R/reduce.R", local = TRUE)
source("R/integrate.R", local = TRUE)
source("R/cluster.R", local = TRUE)
source("R/utils/save_plots.R", local = TRUE)
source("R/plots/plots_normalized_data.R", local = TRUE)
source("R/plots/plots_integrated_data.R", local = TRUE)
source("R/plots/plots_clustered_data.R", local = TRUE)

run_refit_clean_clusters <- function(
  input = "data/job-44542861/clean_annotated_object.rds",
  output = "data/job-44542861/refit_clean_annotated.rds"
) {
  info_path <- paste0(tools::file_path_sans_ext(output), ".INFO")
  stopifnot(!file.exists(output), !file.exists(info_path),
    !file.exists(paste0(output, ".tmp")))
  plot_dir <- file.path("results", "job-44542861", "clusters_plots",
    "refit_clean_annotated")
  stopifnot(!dir.exists(plot_dir))
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  config <- config
  config$integration$harmony$max_iter <- 30L
  validate_config(config)
  future::plan("sequential")
  job_id <- Sys.getenv("SLURM_JOB_ID")
  if (!nzchar(job_id)) {
    job_id <- paste0("manual-", format(Sys.time(), "%Y%m%d-%H%M%S"), "-",
      Sys.getpid())
  }
  run_id <- paste0("job-", job_id)

  info <- c(
    "# Cleaned-object refit",
    "Status: running",
    paste0("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("Job: ", run_id),
    paste0("Input (read-only): ", input),
    paste0("Input record: ", tools::file_path_sans_ext(input), ".INFO"),
    "Pipeline ancestor: job-44542861; cleanup producer 44827588.",
    "Method: reuse existing LogNormalize data; recompute HVGs, scaling/regression, PCA, Harmony, neighbours, Leiden and UMAPs.",
    "No QC, doublet detection or cell-cycle rescoring.",
    paste0("Harmony max_iter: 30 with early stopping (parent job-44542861 used 10; approved departure)."),
    "Seeds: PCA Seurat internal 42 via run_pca_analysis (parent behaviour); Harmony/UMAP seed 1234; clustering seed 1234.",
    paste0("Normalization settings: ", format_value(config$norm)),
    paste0("Integration settings: ", format_value(config$integration)),
    paste0("Clustering settings: ", format_value(config$clustering)),
    paste0("Output: ", output),
    paste0("Figures: ", plot_dir)
  )
  completed <- FALSE
  run_results <- "Validation: pending"
  writeLines(info, info_path)
  on.exit({
    info[info == "Status: running"] <- paste0("Status: ",
      if (completed) "completed" else "failed/incomplete; see execution log")
    writeLines(c(info, paste0("Finished: ",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), run_results), info_path)
  }, add = TRUE)

  input_info <- file.info(input)[, c("size", "mtime")]
  message("Loading (read-only): ", input)
  obj <- readRDS(input)

  # Preflight the actual saved object; never silently repair it.
  expected_counts <- c(
    "B cells" = 26790L, "CD4_T_cells" = 20774L, "CD8_T_cells" = 33565L,
    "DC" = 7436L, "Endothelial" = 13836L, "Epithelial" = 9426L,
    "Granulocytes" = 23037L, "Macrophage_Monocytes" = 41524L,
    "NK_cells" = 19389L, "Stromal" = 11079L,
    "Unconventional_T_and_ILCs" = 27327L
  )
  stopifnot(
    ncol(obj) == 234183L, !anyDuplicated(colnames(obj)),
    identical(as.character(Idents(obj)), as.character(obj$subset)),
    identical(as.integer(table(obj$subset)[names(expected_counts)]),
      unname(expected_counts)),
    !anyNA(obj$sample), length(unique(obj$sample)) > 1L,
    all(config$norm$vars_to_regress %in% colnames(obj[[]])),
    all(c("sample", "cc_tail", "cc_tail_class") %in% names(obj[[]])),
    length(Layers(obj[["RNA"]], search = "^counts($|\\.)")) > 0L,
    length(Reductions(obj)) == 0L, length(Graphs(obj)) == 0L,
    length(Neighbors(obj)) == 0L,
    length(Layers(obj[["RNA"]], search = "^scale\\.data($|\\.)")) == 0L
  )
  cells <- colnames(obj)
  metadata <- obj[[]]
  count_layers <- Layers(obj[["RNA"]], search = "^counts($|\\.)")
  data_layers <- Layers(obj[["RNA"]], search = "^data($|\\.)")
  cleanup_provenance <- obj@misc$cleanup_44542861

  # Preserve any existing partition columns before clustering overwrites them.
  parent_columns <- c(
    grep("_snn_res\\.", names(obj[[]]), value = TRUE),
    "seurat_clusters"
  )
  parent_columns <- intersect(parent_columns, names(obj[[]]))
  for (column in parent_columns) {
    obj[[paste0("parent_", column)]] <- obj[[column, drop = TRUE]]
    obj[[column]] <- NULL
  }

  DefaultAssay(obj) <- "RNA"
  obj@misc$refit_clean <- list(
    input = input, output = output, run_id = run_id,
    cleanup_producer = "44827588", pipeline_ancestor = "job-44542861",
    settings = list(norm = config$norm, dims = config$dims,
      seed_run = config$run$seed, seed_pca_internal = 42L,
      integration = config$integration, clustering = config$clustering),
    parent_analysis = obj@misc[c("integration", "clustering")]
  )

  set.seed(config$run$seed)
  obj <- find_hvgs(obj, config$norm)
  obj <- scale_data(obj, config$norm)
  obj <- run_pca_analysis(obj, config$norm$n_pcs)

  plot_variable_features(obj, plot_dir)
  plot_hvg_overlap(obj, plot_dir)
  plot_pca_elbow(obj, plot_dir)
  plot_pca_loadings(obj, plot_dir, dims = 1:10)

  obj <- integrate_object(obj, config$integration, config$dims,
    config$run$seed)
  stopifnot(grepl("^Harmony converged after",
    obj@misc$integration$harmony_convergence))
  obj <- run_integration_umaps(obj, config$dims, config$run$seed,
    integration_reductions(config$integration$method), config$integration$umap)
  # Same figures as plot_integration_figures(), called per group: that wrapper
  # asserts the plot directory contains exactly its 6 files, which cannot hold
  # here since PCA/variable-feature PNGs already share this directory.
  for (group_by in c("sample", "pressure", "time", "pressure_time")) {
    plot_integration_umaps(obj, group_by = group_by, plot_dir = plot_dir,
      seed = config$run$seed)
  }
  plot_integration_umaps(obj, group_by = "sample", facet_samples = TRUE,
    plot_dir = plot_dir, seed = config$run$seed)

  obj <- cluster_object(obj, config$integration$method, config$clustering,
    config$dims)
  obj <- validate_clustering_umap(obj, "umap_harmony", "harmony", config$dims)

  # Clustering changes active identities; restore the cleaned annotation.
  Idents(obj) <- obj$subset
  obj@misc$clustering$note <- paste(
    "All resolutions retained; active identity restored to the cleaned subset",
    "annotation, not the last partition; seurat_clusters is the last tested",
    "resolution, not a chosen optimum.")

  counts <- cluster_counts(obj)
  obj@misc$refit_clean$session_info <- sessionInfo()

  # Fit genuinely happened: HVGs, PCA, Harmony, UMAPs, graphs, assignments.
  stopifnot(
    length(VariableFeatures(obj)) == config$norm$n_variable_features,
    ncol(Embeddings(obj, "pca")) == config$norm$n_pcs,
    all(is.finite(Embeddings(obj, "pca"))),
    identical(rownames(Embeddings(obj, "pca")), colnames(obj)),
    ncol(Embeddings(obj, "harmony")) == max(config$dims),
    all(is.finite(Embeddings(obj, "harmony"))),
    identical(rownames(Embeddings(obj, "harmony")), colnames(obj)),
    all(c("umap", "umap_harmony") %in% Reductions(obj)),
    identical(colnames(obj), cells),
    identical(as.character(obj$subset), as.character(metadata$subset)),
    identical(obj@misc$cleanup_44542861, cleanup_provenance),
    identical(Layers(obj[["RNA"]], search = "^counts($|\\.)"), count_layers),
    identical(Layers(obj[["RNA"]], search = "^data($|\\.)"), data_layers)
  )
  for (column in obj@misc$clustering$partitions$column) {
    stopifnot(column %in% names(obj[[]]),
      !anyNA(obj[[]][[column]]),
      length(unique(obj[[]][[column]])) > 0L)
  }

  save_final_object(obj, output)
  write_cluster_assignments(obj, file.path(plot_dir, "cluster_assignments.csv"))
  plot_cluster_resolutions(obj, plot_dir)

  reopened <- readRDS(output)
  stopifnot(
    identical(colnames(reopened), cells),
    identical(reopened[[]], obj[[]]),
    identical(Idents(reopened), Idents(obj)),
    identical(reopened@assays, obj@assays),
    identical(reopened@reductions, obj@reductions),
    identical(reopened@graphs, obj@graphs),
    identical(reopened@misc$refit_clean[names(reopened@misc$refit_clean) !=
      "session_info"],
      obj@misc$refit_clean[names(obj@misc$refit_clean) != "session_info"]))
  rm(reopened)
  gc()

  stopifnot(identical(file.info(input)[, c("size", "mtime")], input_info))
  run_results <- c(
    paste0("Actual: ", ncol(obj), " cells; ",
      length(VariableFeatures(obj)), " HVGs; ",
      ncol(Embeddings(obj, "pca")), " PCA components; ",
      obj@misc$integration$harmony_convergence, "."),
    paste0("Leiden clusters: ",
      paste(sprintf("%s at resolution %s", counts, names(counts)),
        collapse = "; ")),
    paste0("Final object saved: ", output),
    "Validation: saved RDS reopened; cell IDs, metadata, identities, assays/layers/counts, reductions, graphs and refit settings match; input size/mtime unchanged."
  )
  completed <- TRUE
  message("Refit completed: ", output, "\nFigures: ", plot_dir,
    "\nSaved RDS verified; cleaned input unchanged.")
  invisible(output)
}

# Sourcing this file defines the runner without starting a run.
if (sys.nframe() == 0L) {
  if (length(commandArgs(trailingOnly = TRUE))) {
    stop("Run without arguments; source this file to set input/output.")
  }
  run_refit_clean_clusters()
}
