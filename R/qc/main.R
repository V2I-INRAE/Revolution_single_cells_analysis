# Run the QC cleaning pipeline over the control samples

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

source("R/qc/params.R")
source("R/qc/io.R")
source("R/qc/filter.R")
source("R/qc/reduce.R")
source("R/qc/doublet_finder.R")
source("R/qc/scrublet.R")
source("R/qc/plots.R")

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
  reduced_obj <- preprocess_for_doublets(obj)
  coordinates <- compute_diagnostic_umap(reduced_obj)
  obj <- AddMetaData(obj, coordinates[colnames(obj), , drop = FALSE])

  # --- Label doublets and copy the calls to the original object
  doublet_obj <- detect_doublets(reduced_obj, sample)
  obj <- add_doublet_calls(obj, doublet_obj)

  # --- Run Scrublet on the same cells; neither method removes cells here.
  obj <- detect_scrublet_doublets(obj)
  labeled_objs[[sample]] <- label_retained_cells(obj)
  rm(obj, reduced_obj, doublet_obj, coordinates)
  unlink(mex_dir, recursive = TRUE)
  gc()
}

labeled <- concatenate_qc_objects(labeled_objs)
rm(labeled_objs)
gc()
write_qc_object(labeled, "data/qc_labeled_data/labeled_concatenated.rds")

# These five figures compare all input cells with QC-only retained cells.
plot_qc_comparison(labeled)
plot_qc_umaps(labeled)

for (method in c("DoubletFinder", "Scrublet")) {
  clean <- filter_labeled_cells(labeled, method)
  write_qc_object(clean, file.path(
    "data/clean_concatenated_data",
    paste0("clean_concatenated_", tolower(method), ".rds")
  ))
  rm(clean)
  gc()
}
