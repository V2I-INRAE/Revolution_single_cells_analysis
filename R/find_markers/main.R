args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript R/find_markers/main.R <clustering/lognorm.rds> <resolutions: 0.1,0.2,...>")
}
source("R/find_markers/params.R")
source("R/find_markers/io.R")
source("R/find_markers/find_markers.R")

suppressPackageStartupMessages(library(presto))

resolutions <- as.numeric(strsplit(args[2], ",", fixed = TRUE)[[1]])
message("Loading marker input: ", args[1])
obj <- read_marker_input(args[1])
message("Input ready: ", ncol(obj), " cells and ", nrow(obj[["RNA"]]), " RNA genes")
job_id <- Sys.getenv("SLURM_JOB_ID")
run_id <- (
  if (nzchar(job_id)) paste0("job-", job_id) else
  paste(format(Sys.time(), "%Y%m%d-%H%M%S"), Sys.getpid(), sep = "-")
)
output_dir <- file.path("results", "find_markers", "cca", run_id)
dir.create(dirname(output_dir), recursive = TRUE, showWarnings = FALSE)

if (!dir.create(output_dir, showWarnings = FALSE)) {
  stop("Cannot create new run directory: ", output_dir)
}
message("Output directory: ", output_dir)

for (resolution in resolutions) {
  message("Processing cluster resolution ", resolution)
  directory <- file.path(output_dir, paste0("res_", format(resolution, nsmall = 1)))
  dir.create(directory)

  markers <- find_all_cluster_markers(
    obj,
    resolution,
    only_pos = only_pos,
    min_pct = min_pct,
    logfc_threshold = logfc_threshold,
    test_use = test_use
  )
  message("Finding markers has been done for resolution ", resolution)
  significant <- dplyr::filter(markers, p_val_adj < adjusted_p_cutoff)
  message(nrow(markers), " cluster–gene results returned; ", nrow(significant),
    " pass adjusted p < ", adjusted_p_cutoff, ". Saving tables.")
  write.csv(
    markers,
    file.path(directory, "markers_prefiltered.csv"),
    row.names = FALSE
  )
  write.csv(
    significant,
    file.path(directory, "markers_significant.csv"),
    row.names = FALSE
  )
  write.csv(
    top_cluster_markers(significant, 25),
    file.path(directory, "markers_top25.csv"),
    row.names = FALSE
  )
  write.csv(
    top_cluster_markers(significant, 5),
    file.path(directory, "markers_top5.csv"),
    row.names = FALSE
  )

  clusters <- obj[[paste0("cca_snn_res.", resolution), drop = TRUE]]
  counts <- table(clusters)

  summary <- data.frame(
    cluster = names(counts),
    n_cells = as.integer(counts),
    n_prefiltered = as.integer(table(markers$cluster)),
    n_significant = as.integer(table(significant$cluster)),
    n_adjusted_p_zero = as.integer(table(markers$cluster[markers$p_val_adj == 0]))
  )

  composition <- as.data.frame(
    table(sample = obj$sample, cluster = clusters),
    responseName = "n_cells"
  )

  write.csv(summary, file.path(directory, "cluster_summary.csv"), row.names = FALSE)
  write.csv(composition, file.path(directory, "sample_cluster_counts.csv"), row.names = FALSE)
  message("Completed cluster resolution ", resolution, "; tables saved to ", directory)
}
message("Marker discovery completed for all ", length(resolutions), " resolutions")
