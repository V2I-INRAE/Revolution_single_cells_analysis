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

# Variable placeholder for the plot summaries
qc_summary <- NULL
doublet_summary <- NULL
scrublet_summary <- NULL
clean_doubletfinder_objs <- list()
clean_scrublet_objs <- list()

for (sample in samples) {
  mex_dir <- read_filtered_matrix("data/raw_data", sample)
  obj <- build_seurat_obj(mex_dir, sample)
  obj <- inspect_seurat_qc(obj)

  obj <- label_qc_outliers(obj)
  cells_passing_qc <- colnames(obj)[!obj$qc_outlier]

  # --- Prepare the temporary reduction and diagnostic UMAP coordinates
  reduced_obj <- preprocess_for_doublets(obj)
  coordinates <- compute_diagnostic_umap(reduced_obj)

  # --- Label doublets and copy the calls to the original object
  doublet_obj <- detect_doublets(reduced_obj, sample)
  obj <- add_doublet_calls(obj, doublet_obj)
  sample_summary <- collect_doublet_summary(obj, sample, coordinates)
  doublet_summary <- rbind(doublet_summary, sample_summary)

  # --- Run Scrublet on the same cells and reuse the diagnostic UMAP coordinates
  obj <- detect_scrublet_doublets(obj)
  sample_summary$doublet_score <- obj$doublet_scores
  sample_summary$doublet_class <- ifelse(
    obj$predicted_doublets, "Doublet", "Singlet"
  )
  scrublet_summary <- rbind(scrublet_summary, sample_summary)

  # --- Collect QC data with the diagnostic UMAP coordinates
  qc_summary <- rbind(
    qc_summary,
    collect_qc_summary(
      obj,
      sample,
      "before",
      coordinates = coordinates
    ),
    collect_qc_summary(
      obj,
      sample,
      "after",
      coordinates = coordinates,
      cells = cells_passing_qc
    )
  )

  # --- Filter independently from the same labeled object
  clean_doubletfinder_objs[[sample]] <- filter_labeled_cells(obj, "DoubletFinder")
  clean_scrublet_objs[[sample]] <- filter_labeled_cells(obj, "Scrublet")
}

# --- Generate the QC, doublet and ribosomal comparison plots across all samples
plot_qc_comparison(qc_summary)
plot_doublet_comparison(doublet_summary)
plot_doublet_nfeature_violin(doublet_summary)
plot_doublet_comparison(scrublet_summary, method = "Scrublet", prefix = "scrublet")
plot_doublet_nfeature_violin(scrublet_summary, method = "Scrublet", prefix = "scrublet")
plot_ribo_comparison(qc_summary)

write_concatenated_obj(clean_doubletfinder_objs, "clean_concatenated_doubletfinder.rds")
rm(clean_doubletfinder_objs)
gc()
write_concatenated_obj(clean_scrublet_objs, "clean_concatenated_scrublet.rds")
