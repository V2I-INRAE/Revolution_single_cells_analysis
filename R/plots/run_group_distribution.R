# Read-only condition, pressure and time-point composition at a saved resolution.
# Usage: Rscript R/plots/run_group_distribution.R <input.rds> <resolution>
# Supports parent clustered/annotated objects and job-specific subclustered objects.
# Annotated inputs use harmony_snn_res.<resolution>, not the subset identities.
suppressPackageStartupMessages(library(SeuratObject))
source("R/utils/save_plots.R")
source("R/plots/plots_group_distribution.R")

# Keep the parent producer and, for subsets, the subclustering producer in the path.
group_distribution_output_dir <- function(input, resolution) {
  data_root <- paste0(normalizePath("data", mustWork = TRUE), "/")
  input <- normalizePath(input, mustWork = TRUE)
  if (!startsWith(input, data_root)) stop("Input must be inside the project data directory.")
  parts <- strsplit(substring(input, nchar(data_root) + 1L), "/", fixed = TRUE)[[1L]]
  if (!grepl("^job-[0-9]+$", parts[[1L]])) stop("Input lacks a parent job directory.")
  if (length(parts) == 2L &&
      (parts[[2L]] %in% c("clustered.rds", "annotated_clustered.rds") ||
        parts[[2L]] == "refit_clean_annotated.rds")) {
    return(file.path("results", parts[[1L]], "clusters_plots",
      "distribution", resolution))
  }
  if (length(parts) == 5L && parts[[2L]] == "subsets" &&
      grepl("^job-[0-9]+$", parts[[4L]]) && parts[[5L]] == "subclustered.rds") {
    return(file.path("results", parts[[1L]], "clusters_plots", "subset",
      parts[[3L]], parts[[4L]], "distribution", resolution))
  }
  stop("Expected data/job-XXX/{clustered,annotated_clustered}.rds or ",
    "data/job-XXX/subsets/<folder>/job-YYY/subclustered.rds.")
}

run_group_distribution <- function(input, resolution) {
  resolution <- suppressWarnings(as.numeric(resolution))
  if (length(resolution) != 1L || !is.finite(resolution) || resolution <= 0) {
    stop("Resolution must be a positive number identifying a saved clustering column.")
  }
  resolution <- as.character(resolution)
  cluster_col <- paste0("harmony_snn_res.", resolution)
  plot_dir <- group_distribution_output_dir(input, resolution)
  if (dir.exists(plot_dir)) {
    stop("Distribution directory already exists; refusing to overwrite: ", plot_dir)
  }
  input_info <- file.info(input)[, c("size", "mtime")]
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  info_path <- file.path(plot_dir, ".INFO")
  info <- c(
    "# Group-distribution figures",
    "Status: running",
    paste0("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("Job: ", Sys.getenv("SLURM_JOB_ID", unset = "manual")),
    paste0("Input (read-only): ", normalizePath(input)),
    paste0("Input bytes: ", input_info$size, "; mtime: ", input_info$mtime),
    paste0("Producer directories: ", paste(regmatches(input,
      gregexpr("job-[0-9]+", input))[[1L]], collapse = ", ")),
    paste0("Resolution: ", resolution, "; cluster column: ", cluster_col),
    "Methods: existing condition/pressure/time-point groups; 100% within-cluster stacked bars; no individual-sample plot or differential-abundance test.",
    "Plot settings: shared viz_palette and canonical group levels; PNG 300 dpi, height 6 inches, width max(10, number of clusters * 0.5).",
    "No filtering, normalization, reclustering or RDS write; no stochastic step/seed applicable.",
    paste0("Output: ", plot_dir)
  )
  completed <- FALSE
  validation <- "Validation: pending"
  writeLines(info, info_path)
  on.exit({
    info[info == "Status: running"] <- paste0("Status: ",
      if (completed) "completed" else "failed/incomplete; see execution log")
    writeLines(c(info, paste0("Finished: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
      validation), info_path)
  }, add = TRUE)

  message("Loading (read-only): ", input)
  obj <- readRDS(input)
  meta <- obj[[]]
  if (!cluster_col %in% names(meta)) stop("Cluster column missing: ", cluster_col)
  stopifnot(
    "Cell IDs missing, duplicated or misaligned" =
      !anyNA(rownames(meta)) && !anyDuplicated(rownames(meta)) &&
        identical(rownames(meta), colnames(obj)),
    "Missing or empty cluster labels" =
      !anyNA(meta[[cluster_col]]) && all(nzchar(as.character(meta[[cluster_col]])))
  )
  meta <- derive_group_columns(meta)
  sample_groups <- unique(meta[c("sample", "pressure", "time_point")])
  stopifnot("Sample IDs have inconsistent pressure/time labels" =
    !anyDuplicated(sample_groups$sample))
  for (group_col in c("condition_group", "pressure_group", "time_group")) {
    stopifnot("Missing group assignments" = !anyNA(meta[[group_col]]))
    counts <- table(factor(meta[[cluster_col]]), meta[[group_col]])
    stopifnot(
      "Group counts do not retain all cells" = sum(counts) == nrow(meta),
      "Empty cluster in group counts" = all(rowSums(counts) > 0),
      "Within-cluster proportions do not sum to one" =
        all(abs(rowSums(prop.table(counts, margin = 1)) - 1) < 1e-12)
    )
  }
  info <- c(info, paste0("Actual: ", nrow(meta), " cells; ",
    length(unique(meta[[cluster_col]])), " clusters; ", length(unique(meta$sample)), " samples."))
  writeLines(info, info_path)
  paths <- plot_group_distributions(obj, cluster_col, plot_dir)
  for (path in unlist(paths)) {
    image <- png::readPNG(path)
    stopifnot("Empty or unreadable figure" = file.size(path) > 0 && all(dim(image)[1:2] > 0))
    rm(image)
  }
  stopifnot("Input changed during plotting" =
    identical(file.info(input)[, c("size", "mtime")], input_info))
  validation <- paste("Validation: all cells assigned; group counts retain all cells;",
    "within-cluster proportions sum to one; three saved PNGs reopened; input size/mtime unchanged.",
    "Visual interpretation remains with owner.")
  completed <- TRUE
  message("Group-distribution figures completed: ", plot_dir)
  invisible(paths)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) {
    stop("Usage: Rscript R/plots/run_group_distribution.R <input.rds> <resolution>")
  }
  run_group_distribution(args[[1L]], args[[2L]])
}
