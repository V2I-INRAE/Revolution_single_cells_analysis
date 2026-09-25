# Run the QC cleaning pipeline over the control samples

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

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

  # -- Building the metadata
  obj$sample <- sample
  sample_parts <- strsplit(sample, "-", fixed = TRUE)[[1]]
  condition <- sample_parts[[2]]
  obj$pig <- sample_parts[[1]]
  obj$pressure <- switch(
    substr(condition, 1, 1),
    C = "None",
    P = "positive",
    N = "negative",
    stop("Unknown pressure code in sample ID: ", sample)
  )
  obj$time_point <- if (condition == "C") {
    "T0H"
  } else {
    paste0("T", substring(condition, 2), "H")
  }

  obj <- inspect_seurat_qc(obj)

  # collect the qc data before filtering
  qc_summary <- rbind(qc_summary, collect_qc_summary(obj, sample, "before"))

  # filter_outliers
  obj <- filter_outliers(obj)

  # collect the qc data after filtering
  qc_summary <- rbind(qc_summary, collect_qc_summary(obj, sample, "after"))

  # compute doublets
  obj_with_doublets <- detect_doublets(obj, sample)

  # append the doublet data to the data frame
  doublet_summary <- rbind(doublet_summary, obj_with_doublets$summary)

  # filter the doublets out of the object in the loop (singlets only)
  obj <- subset(obj_with_doublets$obj, doublet_class == "Singlet")

  # collect the cleaned object for the final concatenation
  clean_objs[[sample]] <- obj
}

plot_qc_comparison(qc_summary)
plot_doublet_comparison(doublet_summary)


concat_dir <- file.path("data", "clean_concatenated_data")
dir.create(concat_dir, showWarnings = FALSE, recursive = TRUE)
concatenated <- merge(clean_objs[[1]], y = clean_objs[-1])
saveRDS(concatenated, file.path(concat_dir, "clean_concatenated.rds"))
