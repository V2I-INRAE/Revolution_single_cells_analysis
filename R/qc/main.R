# Run the QC cleaning pipeline over the control samples

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

source("R/qc/filter.R")

samples <- c(
  "REVO26-C", "REVO26-N4", "REVO26-N10", "REVO26-P4", "REVO26-P10",
  "REVO27-C", "REVO27-N4", "REVO27-N10", "REVO27-P4", "REVO27-P10",
  "REVO29-C", "REVO29-N4", "REVO29-P4",
  "REVO30-C", "REVO30-N4", "REVO30-N10", "REVO30-P10",
  "REV30-P4",
  "REVO31-C", "REVO31-N4", "REVO31-N8", "REVO31-P10",
  "REV31-P4"
)
for (sample in samples) {
  # Use the function implemented before
  mex_dir <- read_filtered_matrix("raw_data", sample)
  obj <- build_seurat_obj(mex_dir, sample)
  obj <- inspect_seurat_qc(obj)

  # plot before filtering
  plot_qc(obj, "before_filtering", sample)

  # filter_outliers
  obj <- filter_outliers(obj)

  # plot after filtering
  plot_qc(obj, "after_filtering", sample)
}
