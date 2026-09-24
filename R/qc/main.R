# Run the QC cleaning pipeline over the control samples

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

source("R/qc/filter.R")
source("R/qc/doublets.R")

samples <- c(
  "REVO26-C", "REVO26-N4", "REVO26-N10", "REVO26-P4", "REVO26-P10",
  "REVO27-C", "REVO27-N4", "REVO27-N10", "REVO27-P4", "REVO27-P10",
  "REVO29-C", "REVO29-N4", "REVO29-P4", "REVO30-C", "REVO30-N4",
  "REVO30-N10", "REVO30-P10", "REV30-P4", "REVO31-C", "REVO31-N4",
  "REVO31-N8", "REVO31-P10", "REV31-P4"
)

# Variable placeholder for the plot summary
doublet_summary <- NULL

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

  # compute doublets
  obj_with_doublets <- detect_doublets(obj, sample)

  # append the doublet data to the data frame
  doublet_summary <- rbind(doublet_summary, obj_with_doublets$summary)

  # filter the doublets out of the object in the loop (singlets only)
  obj <- subset(obj_with_doublets$obj, doublet_class == "Singlet")
}

# outside the loop: comparison plots to keep track of the doublet
# identification across samples
plot_doublet_comparison(doublet_summary)
