args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !tolower(args[[1]]) %in% c("harmony", "cca", "scvi")) {
  stop("Usage: Rscript R/integration/main.R <harmony|cca|scvi>")
}
method <- tolower(args[[1]])

source("R/integration/io.R")
source("R/integration/harmony.R")
source("R/integration/cca.R")
source("R/integration/scvi.R")
source("R/integration/reduce.R")
source("R/integration/diagnostics.R")
source("R/integration/plots.R")

input_dir <- file.path("data", "norm_feat")
job_id <- Sys.getenv("SLURM_JOB_ID")
run_id <- if (nzchar(job_id)) {
  paste0("job-", job_id)
} else {
  paste(format(Sys.time(), "%Y%m%d-%H%M%S"), Sys.getpid(), sep = "-")
}
data_dir <- file.path("data", "integration", method, run_id)
output_dir <- file.path("results", "integration", method, run_id)
routes <- "lognorm"
dims <- 1:30
seed <- 1234
n_diagnostic_cells <- 50000
perplexity <- 30

reductions <- c(umap = "pca", switch(method,
  harmony = c(umap_harmony = "harmony"),
  cca = c(umap_cca = "integrated_cca"),
  scvi = c(umap_scvi = "integrated_scvi")
))

# Reserve both run directories without reusing or overwriting an earlier run.
for (directory in c(data_dir, output_dir)) {
  dir.create(dirname(directory), showWarnings = FALSE, recursive = TRUE)
  if (!dir.create(directory, showWarnings = FALSE)) {
    stop("Cannot create new run directory (exists or is not writable): ", directory)
  }
}
message("Integration: ", method, "; route: lognorm; run: ", run_id,
  "\nData: ", data_dir, "\nResults: ", output_dir)

# Downstream integration uses only LogNormalize; existing SCT files are retained.
for (route in routes) {
  input_file <- file.path(input_dir, paste0(route, ".rds"))
  obj <- read_integration_input(input_file, dims = dims)
  metadata <- obj[[]]
  diagnostic_cells <- sample_diagnostic_cells(
    metadata,
    n_cells = n_diagnostic_cells, seed = seed
  )
  write.csv(data.frame(cell = rownames(diagnostic_cells), diagnostic_cells),
    file.path(output_dir, "diagnostic_cells.csv"),
    row.names = FALSE
  )

  if (method == "harmony") {
    obj <- run_harmony_integration(obj, dims = dims, seed = seed)
  }
  if (method == "cca") {
    obj <- run_cca_integration(obj, dims = dims, seed = seed)
  }
  if (method == "scvi") {
    obj <- run_scvi_integration(obj, dims = dims, seed = seed)
  }
  obj <- run_integration_umaps(
    obj, dims = dims, seed = seed, reductions = reductions
  )

  obj@misc$integration <- list(
    input_file = input_file, route = route, batch_col = "sample",
    method = method, run_id = run_id, reductions = reductions,
    dims = dims, seed = seed, diagnostic_cells = rownames(diagnostic_cells),
    diagnostic_sampling = "proportional by sample, largest remainder",
    perplexity = perplexity,
    purpose = "sensitivity analysis against unintegrated PCA",
    session_info = sessionInfo()
  )
  # Keep the integrated checkpoint even if a subsequent diagnostic fails.
  saveRDS(obj, file.path(data_dir, paste0(route, ".rds")))

  # Cross-method UMAP figures are generated separately once both checkpoints exist.
  diagnostics <- compute_integration_metrics(
    obj,
    cells = rownames(diagnostic_cells), dims = dims, perplexity = perplexity,
    reductions = unname(reductions)
  )
  saveRDS(
    diagnostics, file.path(output_dir, paste0("diagnostics_", route, ".rds"))
  )
  write.csv(diagnostics$summary,
    file.path(output_dir, paste0("diagnostics_", route, ".csv")),
    row.names = FALSE
  )
  plot_integration_metrics(
    diagnostics,
    output_dir = output_dir,
    filename = paste0("integration_metrics_", route, ".png")
  )
  rm(obj, diagnostics)
  gc()
}
