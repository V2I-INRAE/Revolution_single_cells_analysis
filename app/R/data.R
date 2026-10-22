# The app reads metadata once and retrieves one gene at a time from sparse HDF5.
open_explorer_data <- function(config) {
  directory <- config$cache
  metadata <- file.path(directory, "metadata.rds")
  if (!file.exists(metadata)) {
    stop("Viewing data not found. From the project root run: Rscript app/prepare.R")
  }
  if (!"Status: completed" %in% readLines(file.path(directory, ".INFO"))) {
    stop("Viewing-data preparation is incomplete. Check ", file.path(directory, ".INFO"))
  }
  data <- readRDS(metadata)
  data$expression <- HDF5Array::TENxMatrix(
    file.path(directory, "expression.h5"), group = "expression")
  embedding_input <- if (is.null(config$embedding_input)) config$input else
    config$embedding_input
  stopifnot(
    identical(rownames(data$expression), data$cells$cell_id),
    identical(colnames(data$expression), data$genes$id),
    identical(data$source$reduction, config$reduction),
    identical(data$source$cluster_column, config$cluster_column),
    identical(basename(dirname(data$source$input)), basename(dirname(config$input))),
    identical(basename(data$source$input), basename(config$input)),
    identical(basename(dirname(data$source$embedding_input)), basename(dirname(embedding_input))),
    identical(basename(data$source$embedding_input), basename(embedding_input))
  )
  data$label <- config$label
  data$parent <- config$parent
  data
}

read_explorer_gene <- function(data, gene) {
  index <- match(gene, data$genes$id)
  if (is.na(index)) stop("Gene not found: ", gene)
  # Only this cell-length vector becomes dense, never the whole matrix.
  values <- as.numeric(data$expression[, index])
  if (anyNA(values) || any(!is.finite(values)) || any(values < 0)) {
    stop("Invalid normalized RNA values for ", gene)
  }
  unavailable <- !data$availability[index, data$cell_layer]
  list(id = gene, label = data$genes$label[index], values = values,
    maximum = max(values), unavailable = unavailable)
}
