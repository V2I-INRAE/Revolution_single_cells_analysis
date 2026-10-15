# Refit one saved coarse identity using the existing pipeline functions.
# Usage: Rscript R/clusters/subcluster.R <annotated.rds> "<subset label>" <folder>
# Example: data/job-44493014/annotated_clustered.rds "T cells & NK" T_cell_nk
# No normalization, QC filtering or cell-cycle scoring is repeated.

run_subclustering <- function(input, subset_label, folder) {
  parent_job <- basename(dirname(input))
  data_dir <- file.path("data", parent_job, "subsets", folder)
  plot_dir <- file.path("results", parent_job, "clusters_plots", "subset", folder)
  if (any(dir.exists(c(data_dir, plot_dir)))) {
    stop("Subset output directory already exists; preserve it before running again.")
  }

  suppressPackageStartupMessages(library(Seurat))
  source("R/config.R", local = TRUE)
  source("R/lognorm.R", local = TRUE)
  source("R/reduce.R", local = TRUE)
  source("R/integrate.R", local = TRUE)
  source("R/cluster.R", local = TRUE)
  source("R/io.R", local = TRUE) # save_final_object()
  source("R/utils/save_plots.R", local = TRUE)
  source("R/plots/plots_clustered_data.R", local = TRUE)
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
    input = input, subset_label = subset_label,
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

  # Save before plotting, then verify the actual saved subset rather than
  # repeating the original object's input/configuration checks.
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  final_path <- file.path(data_dir, "subclustered.rds")
  save_final_object(obj, final_path)

  reopened <- readRDS(final_path)

  stopifnot(
    identical(colnames(reopened), selected_cells),
    identical(reopened[[]], obj[[]]), identical(Idents(reopened), Idents(obj)),
    identical(reopened@assays, obj@assays),
    identical(reopened@reductions, obj@reductions),
    identical(reopened@graphs, obj@graphs))
  rm(reopened)
  gc()

  plot_cluster_resolutions(obj, plot_dir)
  stopifnot(identical(file.info(input)[, c("size", "mtime")], input_info))
  message("Subclustering completed: ", final_path, "\nFigures: ", plot_dir,
    "\nSaved RDS verified; annotated input unchanged. Figure review left to owner.")
  invisible(final_path)
}

# Sourcing this file defines the runner without starting a run.
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 3L) {
    stop("Usage: Rscript R/clusters/subcluster.R <annotated.rds> \"<subset label>\" <folder>")
  }
  run_subclustering(args[[1L]], args[[2L]], args[[3L]])
}
