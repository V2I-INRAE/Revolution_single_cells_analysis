# Run once from the project root: Rscript app/prepare.R
suppressPackageStartupMessages(library(Seurat))
source("app/config.R")
source("app/R/prepare_data.R")
source("app/R/data.R")
views <- explorer_config(normalizePath("."))
selected <- commandArgs(trailingOnly = TRUE)
if (!length(selected)) selected <- explorer_active_views
if (any(!selected %in% names(views))) {
  stop("Choose views from: ", paste(names(views), collapse = ", "))
}
for (id in selected) {
  config <- views[[id]]
  if (dir.exists(config$cache)) {
    # Reuse only completed caches with the configured reduction and partition.
    open_explorer_data(config)
    message("Keeping prepared view: ", config$label)
  } else {
    prepare_explorer_data(config, normalizePath("."))
  }
}
