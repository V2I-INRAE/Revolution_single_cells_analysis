args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !tolower(args[[1]]) %in% c("harmony", "cca", "scvi")) {
  stop("Usage: Rscript R/integration/main.R <harmony|cca|scvi> <lognorm_input.rds>")
}
method <- tolower(args[[1]])
input_file <- args[[2]]
if (!file.exists(input_file)) stop("Input checkpoint not found: ", input_file)
input_file <- normalizePath(input_file)

source("R/integration/params.R")
source("R/integration/io.R")
source("R/integration/harmony.R")
source("R/integration/cca.R")
source("R/integration/scvi.R")
source("R/integration/reduce.R")
source("R/integration/plots.R")

job_id <- Sys.getenv("SLURM_JOB_ID")
log_file <- list.files("logs",
  pattern = paste0("-integration-", method, "-", job_id, "[.]log$"),
  full.names = TRUE)
run_id <- if (nzchar(job_id)) {
  paste0("job-", job_id)
} else {
  paste(format(Sys.time(), "%Y%m%d-%H%M%S"), Sys.getpid(), sep = "-")
}
data_dir <- file.path("data", "integration", method, run_id)
output_dir <- file.path("results", "integration", method, run_id)
routes <- "lognorm"
dims <- integration_params$dims
scvi_dims <- if (method == "scvi") seq_len(scvi_params$n_latent) else dims
seed <- 1234

reductions <- c(umap = "pca", switch(method,
  harmony = c(umap_harmony = "harmony"),
  cca = c(umap_cca = "integrated_cca"),
  scvi = c(umap_scvi = "integrated_scvi")
))

for (directory in c(data_dir, output_dir)) {
  dir.create(dirname(directory), showWarnings = FALSE, recursive = TRUE)
  if (!dir.create(directory, showWarnings = FALSE)) {
    stop("Cannot create new run directory (exists or is not writable): ", directory)
  }
}
message("Integration: ", method, "; route: lognorm; run: ", run_id,
  "\nInput: ", input_file, "\nData: ", data_dir, "\nResults: ", output_dir)

for (directory in c(data_dir, output_dir)) {
  writeLines(c(
    paste0("# Local provenance — `", directory, "`"),
    "", "## Run",
    paste0("- Producer: ", run_id, "; R/integration/main.R ", method, "."),
    paste0("- Input: ", input_file, "."),
    paste0("- Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "."),
    "", "## Settings",
    paste0("- Sample correction; axes ", paste(range(dims), collapse = "–"),
      "; seed ", seed, "; method ", method, "."),
    if (method == "harmony") paste0("- Harmony max_iter: ", harmony_params$max_iter, "."),
    "- PCA and selected-method UMAPs; sample/pressure/time/pressure_time and sample facets.",
    "- Six 300 dpi PNGs; diagnostic sampling and mixing metrics disabled.",
    "", "## Outputs and status",
    paste0("- Object: ", file.path(data_dir, "lognorm.rds"), "."),
    paste0("- Figures: ", output_dir, "."),
    if (length(log_file)) paste0("- Log: ", paste(log_file, collapse = ", "), ".") else
      "- Log: console output (no job-specific log found).",
    "- RUNNING; completion status follows below."
  ), file.path(directory, ".INFO"))
}

for (route in routes) {
  obj <- read_integration_input(input_file, dims = dims)

  if (method == "harmony") {
    obj <- run_harmony_integration(obj, dims = dims,
      max_iter = harmony_params$max_iter, seed = seed)
    harmony_log <- if (length(log_file) == 1L) readLines(log_file, warn = FALSE) else character()
    convergence <- grep("^Harmony converged after [0-9]+ iterations$", harmony_log, value = TRUE)
    harmony_status <- if (length(convergence)) {
      tail(convergence, 1L)
    } else if (any(harmony_log == paste0("Harmony ", harmony_params$max_iter, "/", harmony_params$max_iter))) {
      paste0("Harmony did not converge within ", harmony_params$max_iter, " iterations")
    } else {
      "Harmony convergence not verified: see execution log"
    }
    for (directory in c(data_dir, output_dir)) {
      write(paste0("- ", harmony_status, "."), file.path(directory, ".INFO"), append = TRUE)
    }
  }
  if (method == "cca") {
    obj <- run_cca_integration(obj, dims = dims, seed = seed)
  }
  if (method == "scvi") {
    obj <- run_scvi_integration(obj, params = scvi_params)
  }
  obj <- run_integration_umaps(
    obj, dims = dims, scvi_dims = scvi_dims, seed = seed,
    reductions = reductions
  )

  obj@misc$integration <- list(
    input_file = input_file, route = route, batch_col = "sample",
    method = method, run_id = run_id, reductions = reductions,
    dims = dims, scvi_dims = if (method == "scvi") scvi_dims else NULL,
    harmony_params = if (method == "harmony") harmony_params else NULL,
    harmony_convergence = if (method == "harmony") harmony_status else NULL,
    seed = seed,
    purpose = "sensitivity analysis against unintegrated PCA",
    session_info = sessionInfo()
  )

  saveRDS(obj, file.path(data_dir, paste0(route, ".rds")))

  for (directory in c(data_dir, output_dir)) {
    write("- Final object saved; UMAP figure rendering started.",
      file.path(directory, ".INFO"), append = TRUE)
  }
  figures <- plot_integration_figures(obj, output_dir = output_dir, seed = seed)
  for (directory in c(data_dir, output_dir)) {
    write(paste0("- COMPLETED: final object and ", length(figures),
      " UMAP PNGs; ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "."),
      file.path(directory, ".INFO"), append = TRUE)
  }
  message("Integration and UMAP figures complete: ", output_dir)

  rm(obj)
  gc()
}
