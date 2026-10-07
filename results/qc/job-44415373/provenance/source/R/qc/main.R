# Run the QC cleaning pipeline over the control samples

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

source("R/qc/params.R")
source("R/qc/io.R")
source("R/qc/filter.R")
source("R/qc/reduce.R")
source("R/qc/scrublet.R")
source("R/qc/plots.R")

job_id <- Sys.getenv("SLURM_JOB_ID")

run_id <- if (nzchar(job_id)) paste0("job-", job_id) else
  paste(format(Sys.time(), "%Y%m%d-%H%M%S"), Sys.getpid(), sep = "-")

  labeled_dir <- file.path("data", "qc_labeled_data", run_id)
clean_dir <- file.path("data", "clean_concatenated_data", run_id)
output_dir <- file.path("results", "qc", run_id)

for (directory in c(labeled_dir, clean_dir, output_dir)) {
  dir.create(dirname(directory), recursive = TRUE, showWarnings = FALSE)
  if (!dir.create(directory, showWarnings = FALSE)) {
    stop("Cannot create new run directory: ", directory)
  }
}

message("QC run: ", run_id, "\nLabeled data: ", labeled_dir,
  "\nClean data: ", clean_dir, "\nResults: ", output_dir)

samples <- c(
  "REVO26-C", "REVO26-N4", "REVO26-N10", "REVO26-P4", "REVO26-P10",
  "REVO27-C", "REVO27-N4", "REVO27-N10", "REVO27-P4", "REVO27-P10",
  "REVO29-C", "REVO29-N4", "REVO29-P4", "REVO30-C", "REVO30-N4",
  "REVO30-N10", "REVO30-P10", "REVO30-P4", "REVO31-C", "REVO31-N4",
  "REVO31-N8", "REVO31-P10", "REVO31-P4"
)

# Keep all input cells until the labeled checkpoint and plots are saved.
labeled_objs <- list()

for (sample in samples) {
  mex_dir <- read_filtered_matrix("data/raw_data", sample)
  obj <- build_seurat_obj(mex_dir, sample)
  obj <- inspect_seurat_qc(obj)

  obj <- label_qc_outliers(obj)

  # --- Prepare the temporary reduction and diagnostic UMAP coordinates
  reduced_obj <- preprocess_qc_umap(obj)
  coordinates <- compute_diagnostic_umap(reduced_obj)
  obj <- AddMetaData(obj, coordinates[colnames(obj), , drop = FALSE])

  # --- Label with Scrublet on all input cells; no removal here.
  obj <- detect_scrublet_doublets(obj)
  labeled_objs[[sample]] <- label_retained_cells(obj)
  rm(obj, reduced_obj, coordinates)
  unlink(mex_dir, recursive = TRUE)
  gc()
}

labeled <- concatenate_qc_objects(labeled_objs)
rm(labeled_objs)
gc()
labeled@misc$qc$run_id <- run_id
write_qc_object(labeled, file.path(labeled_dir, "labeled_concatenated.rds"))

# --- These five figures compare all input cells with QC-only retained cells.
plot_qc_comparison(labeled, plot_dir = output_dir)
plot_qc_umaps(labeled, plot_dir = output_dir)

clean <- filter_labeled_cells(labeled)
write_qc_object(clean, file.path(clean_dir, "clean_concatenated_scrublet.rds"))
