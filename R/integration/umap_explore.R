args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 1:2) {
  stop("Usage: Rscript R/integration/umap_explore.R <harmony_or_cca_checkpoint.rds> [test_name]")
}
input_file <- normalizePath(args[[1]], mustWork = TRUE)
source("R/integration/params.R")
source("R/integration/plots.R")

method <- basename(dirname(dirname(input_file)))
if (!method %in% names(umap_exploration_reductions)) stop("Unsupported integration method: ", method)
method_label <- switch(method, harmony = "Harmony", cca = "CCA")
baseline_name <- paste0("umap_", method)
tests <- names(umap_exploration_tests)
if (method == "harmony") tests <- setdiff(tests, "min.dist_0.5")
array_task <- Sys.getenv("SLURM_ARRAY_TASK_ID")
variant <- if (length(args) == 2L) {
  args[[2]]
} else if (nzchar(array_task)) {
  tests[[as.integer(array_task) + 1L]]
} else {
  "min.dist_0.5"
}
if (!variant %in% names(umap_exploration_tests)) stop("Unknown UMAP test: ", variant)
params <- umap_exploration_params
params$reduction <- umap_exploration_reductions[[method]]
changes <- umap_exploration_tests[[variant]]
params[names(changes)] <- changes
reduction_name <- paste0(baseline_name, "_", variant)
reductions <- setNames(
  c(baseline_name, reduction_name),
  c(paste0(method_label, " baseline\nmin.dist=0.3; neighbours=30; repulsion=1"),
    paste0(method_label, " tested\nmin.dist=", params$min.dist,
      "; neighbours=", params$n.neighbors, "; repulsion=", params$repulsion.strength))
)
output_dir <- file.path("results", "integration", method,
  basename(dirname(input_file)), variant)
parent_records <- if (nzchar(array_task)) character() else
  file.path(c(dirname(input_file), dirname(output_dir)), ".INFO")
job_id <- Sys.getenv("SLURM_JOB_ID")
if (nzchar(array_task)) job_id <- paste(Sys.getenv("SLURM_ARRAY_JOB_ID"), array_task, sep = "_")
log_file <- Sys.getenv("UMAP_LOG_FILE")

if (!dir.exists(output_dir)) {
  if (!dir.create(output_dir)) stop("Cannot create output directory: ", output_dir)
}
records <- file.path(output_dir, ".INFO")
settings <- vapply(params, function(value) paste(value, collapse = ", "), character(1))
for (record in records) {
  write(c(
    "", paste0("## Plots-only ", method_label, " UMAP exploration"),
    paste0("- Input: ", input_file),
    paste0("- Producer: R/integration/umap_explore.R; execution job ", job_id),
    paste0("- Configuration: ", variant),
    paste0("- Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("- Log: ", log_file), "", "## Settings",
    paste0("- Changed parameters: ", paste(names(changes), collapse = ", "), "."),
    paste0("- Reference: saved original ", baseline_name,
      "; min.dist 0.3, n.neighbors 30, repulsion.strength 1."),
    paste0("- ", names(settings), ": ", settings),
    paste0("- New reduction: ", reduction_name),
    paste0("- Panel: ", names(reductions), "; reduction: ", reductions),
    "- Original checkpoint and all existing reductions preserved; no integration rerun.",
    "", "## Outputs and status",
    "- Plots only: no updated object or embedding saved; new reduction exists only in memory.",
    paste0("- Figure: ", file.path(output_dir, "umap_sample_lognorm.png")),
    paste0("- One 300 dpi sample-coloured PNG: saved ", method_label,
      " baseline versus tested ", method_label, " UMAP; original plotting settings."),
    "- Sample PNG replacement authorized; other existing PNGs are unchanged historical outputs. Earlier .INFO entries retained.",
    "- No visual inspection requested for this exploration.",
    "- RUNNING; no separate object-reading or validation job."
  ), record, append = TRUE)
}
for (record in parent_records) {
  write(c("", paste0("## UMAP exploration — ", variant),
    paste0("- Execution job: ", job_id, "; log: ", log_file),
    paste0("- Reduction: ", reduction_name, "; original object/reductions unchanged."),
    "- Plots only; no updated object or embedding saved.",
    paste0("- Figures and run record: ", output_dir),
    "- RUNNING; completion status follows below."), record, append = TRUE)
}

tryCatch({
  obj <- readRDS(input_file)
  stopifnot("Checkpoint method must match its input directory" = obj@misc$integration$method == method,
    "Checkpoint must contain integration coordinates and the saved baseline UMAP" =
      all(c(params$reduction, baseline_name) %in% Reductions(obj)))
  obj <- run_umap_exploration(obj, params, reduction_name)
  for (record in records) {
    write("- New UMAP ready in memory; figure rendering started.", record, append = TRUE)
  }
  plot_integration_umaps(obj, group_by = "sample", output_dir = output_dir,
    seed = params$seed.use, filename = "umap_sample_lognorm.png",
    reductions = reductions)
  status <- paste0("- COMPLETED: one sample comparison PNG; no object/embedding saved; ",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  for (record in c(records, parent_records)) write(status, record, append = TRUE)
  message(status)
}, error = function(error) {
  status <- paste0("- FAILED: ", conditionMessage(error), "; ",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  for (record in c(records, parent_records)) write(status, record, append = TRUE)
  stop(error)
})
