args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript R/clustering/main.R <unintegrated|harmony|cca|scvi> <input.rds>")
}

method <- match.arg(tolower(args[1]), c("unintegrated", "harmony", "cca", "scvi"))

source("R/clustering/params.R")
source("R/clustering/io.R")
source("R/clustering/cluster.R")
source("R/clustering/diagnostics.R")

future::plan("sequential")

settings <- clustering_params
settings$rng_kind <- RNGkind()
job_id <- Sys.getenv("SLURM_JOB_ID")
run_id <- if (nzchar(job_id)) paste0("job-", job_id) else
  paste(format(Sys.time(), "%Y%m%d-%H%M%S"), Sys.getpid(), sep = "-")
data_dir <- file.path("data", "clustering", method, run_id)
output_dir <- file.path("results", "clustering", method, run_id)
for (directory in c(data_dir, output_dir)) {
  dir.create(dirname(directory), recursive = TRUE, showWarnings = FALSE)
  if (!dir.create(directory, showWarnings = FALSE)) {
    stop("Cannot create new run directory: ", directory)
  }
}

message("Clustering: ", method, "; run: ", run_id)

obj <- read_clustering_input(args[2], method, clustering_dims(settings, method))
provenance <- clustering_provenance(obj, args[2])
diagnostic_cells <- sample_clustering_cells(obj[[]], settings$n_diagnostic_cells, settings$seed)

obj <- cluster_reduction(obj, method, settings)
obj <- clustering_umap(obj)
obj@misc$clustering$provenance <- provenance
obj@misc$clustering$run_id <- run_id

diagnostics <- evaluate_clustering(obj, diagnostic_cells)

save_clustering_results(obj, diagnostics, provenance, data_dir, output_dir)
