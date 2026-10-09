args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript R/integration/umap_explore.R <harmony_checkpoint.rds>")
}
input_file <- normalizePath(args[[1]], mustWork = TRUE)
source("R/integration/params.R")
source("R/integration/plots.R")

params <- umap_exploration_params
variant <- paste0("min.dist_", params$min.dist)
reduction_name <- paste0("umap_", params$reduction, "_", variant)
data_dir <- file.path(dirname(input_file), variant)
output_dir <- file.path("results", "integration", "harmony",
  basename(dirname(input_file)), variant)
parent_records <- file.path(c(dirname(input_file), dirname(output_dir)), ".INFO")
job_id <- Sys.getenv("SLURM_JOB_ID")
log_file <- Sys.getenv("UMAP_LOG_FILE")

if (any(dir.exists(c(data_dir, output_dir)))) {
  stop("Exploration output directory already exists; refusing to overwrite: ", variant)
}
for (directory in c(data_dir, output_dir)) {
  if (!dir.create(directory)) stop("Cannot create output directory: ", directory)
}
records <- file.path(c(data_dir, output_dir), ".INFO")
settings <- vapply(params, function(value) paste(value, collapse = ", "), character(1))
for (record in records) {
  writeLines(c(
    "# Harmony UMAP parameter exploration", "", "## Inputs and run",
    paste0("- Input: ", input_file),
    paste0("- Producer: R/integration/umap_explore.R; execution job ", job_id),
    paste0("- Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("- Log: ", log_file), "", "## Settings",
    paste0("- Sole change from original UMAP: min.dist 0.3 → ", params$min.dist, "."),
    paste0("- ", names(settings), ": ", settings),
    paste0("- New reduction: ", reduction_name),
    "- Original checkpoint and all existing reductions preserved; no integration rerun.",
    "", "## Outputs and status",
    paste0("- Extended checkpoint: ", file.path(data_dir, "lognorm.rds")),
    paste0("- Figures: ", output_dir),
    "- Six 300 dpi PNGs; unchanged PCA baseline and original plotting settings.",
    "- RUNNING; no separate object-reading or validation job."
  ), record)
}
for (record in parent_records) {
  write(c("", paste0("## UMAP exploration — ", variant),
    paste0("- Execution job: ", job_id, "; log: ", log_file),
    paste0("- Reduction: ", reduction_name, "; original object/reductions unchanged."),
    paste0("- Checkpoint and run record: ", data_dir),
    paste0("- Figures and run record: ", output_dir),
    "- RUNNING; completion status follows below."), record, append = TRUE)
}

tryCatch({
  obj <- readRDS(input_file)
  stopifnot("Expected a Harmony checkpoint" = obj@misc$integration$method == "harmony")
  if (reduction_name %in% Reductions(obj)) {
    stop("Reduction already exists; refusing to overwrite: ", reduction_name)
  }
  obj <- do.call(Seurat::RunUMAP, c(list(
    object = obj, reduction.name = reduction_name,
    reduction.key = paste0(gsub("[^[:alnum:]]", "", reduction_name), "_"),
    verbose = TRUE
  ), params))
  obj@misc$umap_explorations[[reduction_name]] <- list(
    input_file = input_file, job_id = job_id, params = params,
    data_dir = data_dir, output_dir = output_dir, session_info = sessionInfo()
  )
  saveRDS(obj, file.path(data_dir, "lognorm.rds"))
  for (record in records) {
    write("- Extended checkpoint saved; figure rendering started.", record, append = TRUE)
  }
  figures <- plot_integration_figures(obj, output_dir = output_dir,
    seed = params$seed.use, method_reduction = reduction_name)
  status <- paste0("- COMPLETED: extended checkpoint and ", length(figures),
    " PNGs; ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  for (record in c(records, parent_records)) write(status, record, append = TRUE)
  message(status)
}, error = function(error) {
  status <- paste0("- FAILED: ", conditionMessage(error), "; ",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  for (record in c(records, parent_records)) write(status, record, append = TRUE)
  stop(error)
})
