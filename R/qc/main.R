# Run the QC cleaning pipeline over the control samples

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

source("R/qc/filter.R")

samples <- c("REVO26-C", "REVO27-C", "REVO29-C", "REVO30-C", "REVO31-C")
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
