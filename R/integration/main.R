# Sample-level correction is a sensitivity analysis, not a replacement for the
# unintegrated reference. Pressure and time are never correction variables.
source("R/integration/io.R")
source("R/integration/harmony.R")
source("R/integration/cca.R")
source("R/integration/reduce.R")
source("R/integration/diagnostics.R")
source("R/norm_feat/plots.R")
source("R/integration/plots.R")

input_dir <- file.path("data", "norm_feat")
data_dir <- file.path("data", "integration")
output_dir <- file.path("results", "integration")
routes <- c("lognorm", "sct")
dims <- 1:30
seed <- 1234
n_diagnostic_cells <- 20000
perplexity <- 30

does_harmony <- TRUE
does_cca <- FALSE

reductions <- c(umap = "pca")
if (does_harmony) reductions["umap_harmony"] <- "harmony"
if (does_cca) reductions["umap_cca"] <- "integrated_cca"

dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
diagnostic_cells <- NULL
reference_metadata <- NULL

# Process routes sequentially to avoid retaining both large objects in memory.
for (route in routes) {
  input_file <- file.path(input_dir, paste0(route, ".rds"))
  obj <- read_integration_input(input_file, dims = dims)
  metadata <- obj[[]]
  if (is.null(diagnostic_cells)) {
    reference_metadata <- metadata
    diagnostic_cells <- sample_diagnostic_cells(
      metadata,
      n_cells = n_diagnostic_cells, seed = seed
    )
    write.csv(data.frame(cell = rownames(diagnostic_cells), diagnostic_cells),
      file.path(output_dir, "diagnostic_cells.csv"),
      row.names = FALSE
    )
  } else {
    stopifnot(
      "Normalization routes must contain the same cells and sample labels" =
        setequal(rownames(metadata), rownames(reference_metadata)),
      identical(
        as.character(metadata[rownames(reference_metadata), "sample"]),
        as.character(reference_metadata$sample)
      )
    )
  }

  if (does_harmony) {
    obj <- run_harmony_integration(obj, dims = dims, seed = seed)
  }
  if (does_cca) {
    obj <- run_cca_integration(obj, dims = dims, seed = seed)
  }
  obj <- run_integration_umaps(
    obj, dims = dims, seed = seed, reductions = reductions
  )

  obj@misc$integration <- list(
    input_file = input_file, route = route, batch_col = "sample",
    reductions = reductions,
    dims = dims, seed = seed, diagnostic_cells = rownames(diagnostic_cells),
    diagnostic_sampling = "proportional by sample, largest remainder",
    perplexity = perplexity,
    purpose = "sensitivity analysis against unintegrated PCA",
    session_info = sessionInfo()
  )
  # Keep the integrated checkpoint even if a subsequent diagnostic fails.
  saveRDS(obj, file.path(data_dir, paste0(route, ".rds")))

  for (group_by in c("sample", "pressure", "time", "pressure_time")) {
    plot_integration_umaps(
      obj,
      group_by = group_by, output_dir = output_dir, seed = seed,
      filename = paste0("umap_", group_by, "_", route, ".png")
    )
  }
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
