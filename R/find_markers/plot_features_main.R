args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript R/find_markers/plot_features_main.R <clustering/lognorm.rds> <resolution> <all or lineages: T_cells,Macrophages>")
}
source("R/find_markers/params.R")
source("R/find_markers/plots.R")

suppressPackageStartupMessages(library(Seurat))

resolution <- as.numeric(args[2])
lineages <- if (args[3] == "all") names(markers) else strsplit(args[3], ",", fixed = TRUE)[[1]]
if (length(lineages) == 0L) {
  stop("No lineages requested")
}
unknown <- setdiff(lineages, names(markers))
if (length(unknown) > 0L) {
  stop("Unknown lineages not in params.R markers: ", paste(unknown, collapse = ", "))
}

message("Loading object: ", args[1])
obj <- readRDS(args[1])
DefaultAssay(obj) <- "RNA"

# --- scplotter/GetAssayData cannot handle the per-sample split layers of the v5 assay.
obj <- JoinLayers(obj, assay = "RNA", layers = "data")
Idents(obj) <- obj[[paste0("cca_snn_res.", resolution), drop = TRUE]]
message(
  "Input ready: ", ncol(obj), " cells, ", length(levels(Idents(obj))),
  " clusters at resolution ", resolution
)

job_id <- Sys.getenv("SLURM_JOB_ID")
run_id <- (
  if (nzchar(job_id)) {
    paste0("job-", job_id)
  } else {
    paste(format(Sys.time(), "%Y%m%d-%H%M%S"), Sys.getpid(), sep = "-")
  }
)
run_dir <- file.path("results", "find_markers", "cca", run_id)
output_dir <- file.path(run_dir, paste0("res_", format(resolution, nsmall = 1)))
# Guard the run directory, not just the resolution subdirectory: a second
# invocation under the same SLURM_JOB_ID would otherwise overwrite the .INFO.
dir.create(dirname(run_dir), recursive = TRUE, showWarnings = FALSE)
if (!dir.create(run_dir, showWarnings = FALSE)) {
  stop("Run directory already exists: ", run_dir)
}
dir.create(output_dir)
info_path <- file.path(run_dir, "plot_features.INFO")
writeLines(c(
  "# Lineage marker feature plots (scplotter)",
  "## Inputs",
  paste("- object:", args[1]),
  paste("- producer:", run_id),
  if (nzchar(job_id)) paste("- log:", Sys.getenv("LINEAGE_PLOTS_LOG")),
  paste("- resolution:", resolution),
  paste("- lineages:", paste(lineages, collapse = ", ")),
  unlist(lapply(lineages, function(lineage) {
    paste0("- ", lineage, " markers: ", paste(markers[[lineage]], collapse = ", "))
  })),
  "## Settings",
  "- assay RNA, LogNormalize data layer; reduction umap_cca; marker list from R/find_markers/params.R",
  "- detection CSV: percentage of all cells in each cluster with RNA data > 0; rounded to two decimals",
  "## Status",
  paste("- started:", format(Sys.time()))
), info_path)
message("Output directory: ", output_dir)

for (lineage in lineages) {
  message("Plotting ", lineage, " (", length(markers[[lineage]]), " markers)")
  plot_marker_features(obj, markers[[lineage]], resolution, output_dir, label = lineage)
  pct <- write_marker_detection(obj, markers[[lineage]], resolution,
    output_dir, label = lineage)
  if (!is.null(pct)) {
    message("Detection percentage (non-zero expression) written to detection_", lineage, ".csv")
    print(round(pct * 100, 1))
  }
  absent <- setdiff(markers[[lineage]], rownames(obj[["RNA"]]))
  if (length(absent) > 0L) {
    cat(paste0("- ", lineage, " markers not in object, skipped: ",
      paste(absent, collapse = ", "), "\n"),
      file = info_path, append = TRUE)
  }
}
expected <- unlist(lapply(lineages, function(lineage) {
  file.path(output_dir, c(
    paste0("features_", lineage, ".png"),
    paste0("detection_", lineage, ".csv")
  ))
}))
bad <- expected[!file.exists(expected) | file.size(expected) == 0]
if (length(bad) > 0L) {
  stop("Missing or empty outputs: ", paste(bad, collapse = ", "))
}
cat(
  "- completed: ", format(Sys.time()), "\n## Validation\n",
  paste0("- ", list.files(output_dir), collapse = "\n"), "\n",
  file = info_path, append = TRUE, sep = ""
)
message("Done; outputs in ", output_dir)
