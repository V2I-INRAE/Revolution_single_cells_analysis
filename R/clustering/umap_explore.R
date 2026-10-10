args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 1:2) {
  stop("Usage: Rscript R/clustering/umap_explore.R <cca_clustered_lognorm.rds> [test_name]")
}
input_file <- normalizePath(args[[1]], mustWork = TRUE)
stopifnot("Only CCA clustering checkpoints are supported" =
  basename(dirname(dirname(input_file))) == "cca")
suppressPackageStartupMessages(library(Seurat))
source("R/utils/umap_explore.R")
source("R/clustering/plots.R")

tests <- names(umap_exploration_tests)
array_task <- Sys.getenv("SLURM_ARRAY_TASK_ID")
variant <- if (length(args) == 2L) {
  args[[2]]
} else if (nzchar(array_task)) {
  tests[[as.integer(array_task) + 1L]]
} else {
  "min.dist_0.5"
}
if (!variant %in% tests) stop("Unknown UMAP test: ", variant)
params <- umap_exploration_params
params$reduction <- "integrated_cca"
changes <- umap_exploration_tests[[variant]]
params[names(changes)] <- changes
reduction_name <- paste0("umap_cca_", variant)
run_id <- basename(dirname(input_file))
parent_dir <- file.path("results", "clustering", "cca", run_id, "umap_explore")
dir.create(parent_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- file.path(parent_dir, variant)
if (!dir.create(output_dir, showWarnings = FALSE)) {
  stop("Cannot create new variant directory; refusing existing outputs: ", output_dir)
}
record <- file.path(output_dir, ".INFO")
job_id <- Sys.getenv("SLURM_JOB_ID")
if (nzchar(array_task)) job_id <- paste(Sys.getenv("SLURM_ARRAY_JOB_ID"), array_task, sep = "_")
format_settings <- function(settings) {
  vapply(settings, function(value) {
    if (is.null(value)) "NULL (automatic)" else paste(value, collapse = ", ")
  }, character(1))
}
settings <- format_settings(params)
writeLines(c(
  "# CCA cluster-coloured UMAP exploration", "", "## Input and producer",
  paste0("- Input: ", input_file),
  paste0("- Producer: R/clustering/umap_explore.R; job ", job_id),
  paste0("- Log: ", Sys.getenv("UMAP_LOG_FILE")),
  paste0("- Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "", "## Tested settings", paste0("- Configuration: ", variant),
  paste0("- ", names(settings), ": ", settings),
  "- Fixed saved Leiden partitions: 0.1, 0.2, 0.3, 0.4; no graph/clustering/integration refit.",
  "- Plots only: four paired PNGs; no object, coordinates or analysis bundle saved.",
  "", "## Status and validation", "- RUNNING; visual inspection pending."
), record)

tryCatch({
  input_stat <- file.info(input_file)[, c("size", "mtime")]
  obj <- readRDS(input_file)
  config <- obj@misc$clustering
  stopifnot("Input must document CCA clustering, matching run ID and saved UMAP" =
    identical(config$method, "cca") && identical(config$run_id, run_id) &&
      identical(config$reduction, params$reduction) && identical(config$umap, "umap_cca"),
    "All four requested saved resolutions are required" =
      identical(config$partitions$resolution, c(0.1, 0.2, 0.3, 0.4)))
  columns <- config$partitions$column
  stopifnot("Saved partition columns must be unique, complete factors aligned with cells" =
    !anyDuplicated(columns) && all(columns %in% colnames(obj[[]])) &&
      identical(rownames(obj[[]]), colnames(obj)) && !anyDuplicated(colnames(obj)) &&
      all(vapply(obj[[]][, columns, drop = FALSE], function(x) {
        is.factor(x) && !anyNA(x) && nlevels(x) == length(unique(x))
      }, logical(1))))
  baseline <- config$umap_settings
  commands <- Filter(function(command) identical(command@params$reduction.name, config$umap),
    obj@commands)
  stopifnot("Saved baseline settings and command must agree" =
    !is.null(baseline) && length(commands) == 1L && identical(baseline, commands[[1]]@params))
  # Automatic baseline epochs are recorded as such, not replaced with today's defaults.
  reference <- c(umap_exploration_params, list(reduction = params$reduction))
  for (name in setdiff(names(reference), "n.epochs")) {
    if (!isTRUE(all.equal(baseline[[name]], reference[[name]], check.attributes = FALSE))) {
      stop("Baseline differs from the approved reference at ", name, "; review before fitting")
    }
  }
  stopifnot("An explicit baseline epoch count must match the exploration" =
    is.null(baseline$n.epochs) || baseline$n.epochs == params$n.epochs)
  baseline_coordinates <- Embeddings(obj, config$umap)
  stopifnot("Saved baseline UMAP must be finite, two-dimensional and aligned" =
    identical(rownames(baseline_coordinates), colnames(obj)) &&
      ncol(baseline_coordinates) == 2L && all(is.finite(baseline_coordinates)))
  rm(baseline_coordinates)
  preserved_hashes <- function(obj) {
    c(cca = digest::digest(Embeddings(obj, config$reduction), algo = "sha256"),
      baseline = digest::digest(Embeddings(obj, config$umap), algo = "sha256"),
      labels = digest::digest(obj[[]][, columns, drop = FALSE], algo = "sha256"))
  }
  before <- preserved_hashes(obj)
  counts <- vapply(obj[[]][, columns, drop = FALSE], nlevels, integer(1))
  baseline_settings <- format_settings(setNames(
    lapply(names(reference), function(name) baseline[[name]]), names(reference)))
  write(c("", "## Saved baseline and preservation",
    paste0("- ", names(baseline_settings), ": ", baseline_settings),
    "- Automatic baseline epochs require the original integration log for their effective count.",
    paste0("- Cells: ", ncol(obj), "; clusters at 0.1/0.2/0.3/0.4: ", paste(counts, collapse = "/")),
    paste0("- SHA256 ", names(before), ": ", before),
    paste0("- Packages: ", paste(vapply(c("Seurat", "SeuratObject", "uwot", "scplotter"),
      function(x) paste(x, packageVersion(x)), character(1)), collapse = "; "))
  ), record, append = TRUE)
  # Capture baseline provenance above: Seurat can replace the input-reduction command log.
  obj <- run_umap_exploration(obj, params, reduction_name)
  title <- function(label, settings) {
    paste0(label, "\nmin.dist=", settings$min.dist, "; neighbours=", settings$n.neighbors,
      "; repulsion=", settings$repulsion.strength)
  }
  reductions <- setNames(c(config$umap, reduction_name),
    c(title("Saved CCA baseline", baseline), title("Tested CCA UMAP", params)))
  plot_cluster_umap_comparisons(obj, reductions, output_dir)
  stopifnot("CCA coordinates, baseline UMAP and fixed labels must remain unchanged" =
    identical(before, preserved_hashes(obj)),
    "Source checkpoint size and modification time must remain unchanged" =
      identical(input_stat, file.info(input_file)[, c("size", "mtime")]))
  figures <- paste0("umap_clusters_res_", config$partitions$resolution, ".png")
  stopifnot("All four figures must exist and be nonempty" =
    all(file.exists(file.path(output_dir, figures))) &&
      all(file.size(file.path(output_dir, figures)) > 0))
  write(c("", "## Completion",
    paste0("- COMPLETED: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "- PASS: finite/aligned tested embedding; unchanged CCA, baseline and labels (SHA256).",
    "- PASS: input file size/mtime unchanged; four nonempty PNGs; visual inspection pending.",
    paste0("- Output: ", figures)), record, append = TRUE)
  message("COMPLETED: ", output_dir)
}, error = function(error) {
  write(paste0("- FAILED: ", conditionMessage(error), "; ", Sys.time()), record, append = TRUE)
  stop(error)
})
