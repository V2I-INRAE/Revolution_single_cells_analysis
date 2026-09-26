# Run the QC cleaning pipeline over the control samples

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

source("R/qc/io.R")
source("R/qc/filter.R")
source("R/qc/doublets.R")
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
clean_objs <- list()

for (sample in samples) {
  mex_dir <- read_filtered_matrix("data/raw_data", sample)
  obj <- build_seurat_obj(mex_dir, sample)
  obj <- inspect_seurat_qc(obj)

  # --- Collect QC data before and after outlier labeling for later plotting
  qc_summary <- rbind(qc_summary, collect_qc_summary(obj, sample, "before"))
  obj <- label_qc_outliers(obj)
  cells_passing_qc <- colnames(obj)[!obj$qc_outlier]
  qc_summary <- rbind(
    qc_summary,
    collect_qc_summary(obj, sample, "after", cells_passing_qc)
  )

  # --- Label doublets on the original object and collect data for plotting
  obj_with_doublets <- detect_doublets(obj, sample)
  doublet_summary <- rbind(doublet_summary, obj_with_doublets$summary)

  # --- We filter QC outliers and doublets in one operation
  obj <- filter_labeled_cells(obj_with_doublets$obj)

  # --- We collect the cleaned object for the final concatenation
  clean_objs[[sample]] <- obj
}

# --- Generate the QC and doublet comparison plots across all samples
plot_qc_comparison(qc_summary)
plot_doublet_comparison(doublet_summary)

write_concatenated_obj(clean_objs)
