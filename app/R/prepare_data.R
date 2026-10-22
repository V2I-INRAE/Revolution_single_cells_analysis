# Preparation runs separately from Shiny and never writes to the Seurat source.
prepare_explorer_data <- function(config, project_dir) {
  if (dir.exists(config$cache)) {
    stop("Cache directory already exists: ", config$cache,
      ". Keep it for reuse, or move it aside before rebuilding.")
  }
  input_info <- file.info(config$input)[, c("size", "mtime")]
  if (is.na(input_info$size)) stop("Input not found: ", config$input)
  embedding_input <- if (is.null(config$embedding_input)) config$input else
    config$embedding_input
  embedding_info <- file.info(embedding_input)[, c("size", "mtime")]
  if (is.na(embedding_info$size)) stop("UMAP input not found: ", embedding_input)
  dir.create(config$cache, recursive = TRUE)
  info_file <- file.path(config$cache, ".INFO")
  started <- Sys.time()
  info <- c("# Local single-cell explorer preparation", "Status: running",
    paste("Started:", started), paste("Input (read-only):", config$input),
    paste("UMAP input (read-only):", embedding_input),
    paste("Source producer:", config$producer),
    paste("Producer: local Rscript app/prepare.R; PID", Sys.getpid()),
    paste("Output:", config$cache),
    paste("Reduction:", config$reduction),
    paste("Clusters:", config$cluster_column),
    "Expression: RNA data layers; existing LogNormalize, scale factor 10000.",
    "No normalization, filtering, refitting, or stochastic steps; seed not applicable.",
    "Storage: HDF5Array::writeTENxMatrix, cells x genes, sparse, double precision.",
    paste("Packages:", paste(vapply(c("Seurat", "Matrix", "HDF5Array"),
      function(package) paste(package, utils::packageVersion(package)), character(1)),
      collapse = "; ")),
    "Gene labels: existing marker_feature_renames; original feature IDs preserved.",
    "Missing sample features: joined as zero, with original availability retained.")
  writeLines(info, info_file)
  completed <- FALSE
  on.exit({
    info[info == "Status: running"] <- if (completed) "Status: completed" else
      "Status: failed/incomplete; see preparation log"
    writeLines(c(info, paste("Finished:", Sys.time())), info_file)
  }, add = TRUE)

  source(file.path(project_dir, "R/utils/save_plots.R"), local = TRUE)
  source(file.path(project_dir, "R/clusters/markers.R"), local = TRUE)
  embedding <- NULL
  if (!identical(embedding_input, config$input)) {
    message("Loading read-only UMAP source: ", embedding_input)
    embedding_obj <- readRDS(embedding_input)
    embedding <- SeuratObject::Embeddings(embedding_obj, reduction = config$reduction)
    rm(embedding_obj)
    gc()
  }
  message("Loading read-only source: ", config$input)
  obj <- readRDS(config$input)
  cells <- obj[[]]
  if (is.null(embedding)) {
    embedding <- SeuratObject::Embeddings(obj, reduction = config$reduction)
  }
  stopifnot(!anyDuplicated(rownames(embedding)),
    all(colnames(obj) %in% rownames(embedding)))
  embedding_cells <- nrow(embedding)
  # Cleaned objects can omit cells or use a different order than the UMAP source.
  embedding <- embedding[colnames(obj), , drop = FALSE]
  stopifnot(
    identical(rownames(cells), colnames(obj)),
    identical(rownames(embedding), colnames(obj)),
    !anyDuplicated(rownames(cells)), !anyNA(embedding),
    all(is.finite(embedding)),
    all(c(config$cluster_column, "sample", "pig", "pressure",
      "time_point") %in% names(cells))
  )
  cells <- data.frame(cell_id = rownames(cells),
    umap_1 = embedding[, 1L], umap_2 = embedding[, 2L],
    cluster = cells[[config$cluster_column]],
    sample = cells$sample, pig = cells$pig, pressure = cells$pressure,
    time_point = cells$time_point, row.names = NULL)
  stopifnot(!anyNA(cells))
  cells$cluster <- droplevels(factor(cells$cluster))
  groups <- levels(cells$cluster)
  stopifnot(length(groups) <= length(viz_palette))
  colors <- list(cluster = stats::setNames(viz_palette[seq_along(groups)], groups))

  layers <- SeuratObject::Layers(obj[["RNA"]], search = "^data($|\\.)")
  if (!length(layers)) stop("Normalized RNA data layers are missing.")
  feature_ids <- rownames(obj[["RNA"]])
  availability <- matrix(FALSE, nrow = length(feature_ids), ncol = length(layers),
    dimnames = list(feature_ids, layers))
  cell_layer <- integer(nrow(cells))
  for (i in seq_along(layers)) {
    layer <- SeuratObject::LayerData(obj, assay = "RNA", layer = layers[i])
    indices <- match(colnames(layer), cells$cell_id)
    stopifnot(!anyNA(indices), all(cell_layer[indices] == 0L))
    cell_layer[indices] <- i
    availability[match(rownames(layer), feature_ids), i] <- TRUE
  }
  rm(layer)
  stopifnot(all(cell_layer > 0L))
  obj <- Seurat::DietSeurat(obj, assays = "RNA", layers = layers,
    dimreducs = NULL, graphs = NULL, misc = FALSE)
  obj <- SeuratObject::JoinLayers(obj, assay = "RNA",
    layers = "^data($|\\.)", new = "data")
  expression <- SeuratObject::LayerData(obj, assay = "RNA", layer = "data")
  stopifnot(identical(colnames(expression), cells$cell_id),
    identical(rownames(expression), feature_ids), !anyDuplicated(feature_ids),
    !anyNA(expression@x), all(is.finite(expression@x)), all(expression@x >= 0))
  genes <- data.frame(id = feature_ids, symbol = feature_ids)
  translated <- match(feature_ids, names(marker_feature_renames))
  genes$symbol[!is.na(translated)] <- unname(marker_feature_renames[
    translated[!is.na(translated)]])
  genes$label <- ifelse(genes$id == genes$symbol, genes$id,
    paste0(genes$symbol, " [", genes$id, "]"))
  message("Prepared ", nrow(cells), " cells and ", nrow(genes), " genes.")

  message("Writing sparse, gene-indexed expression.")
  expression_file <- file.path(config$cache, "expression.h5")
  HDF5Array::writeTENxMatrix(Matrix::t(expression), expression_file,
    group = "expression", level = 6L, verbose = TRUE)
  metadata <- list(cells = cells, genes = genes, colors = colors,
    availability = availability, cell_layer = cell_layer,
    source = list(input = config$input, producer = config$producer,
      embedding_input = embedding_input, embedding_size = embedding_info$size,
      embedding_mtime = embedding_info$mtime,
      reduction = config$reduction, cluster_column = config$cluster_column,
      size = input_info$size, mtime = input_info$mtime))

  # Publish metadata only after checking the saved expression and all ID orders.
  saved <- HDF5Array::TENxMatrix(expression_file, group = "expression")
  stopifnot(identical(rownames(saved), cells$cell_id),
    identical(colnames(saved), genes$id))
  # Compare every feature in bounded blocks; no dense whole-matrix conversion.
  for (first in seq.int(1L, nrow(expression), by = 128L)) {
    selected <- first:min(first + 127L, nrow(expression))
    actual <- as(saved[, selected, drop = FALSE], "dgCMatrix")
    # HDF5 omits explicit zero slots; drop0 changes storage, not matrix values.
    expected <- Matrix::drop0(Matrix::t(expression[selected, , drop = FALSE]))
    stopifnot(isTRUE(all.equal(actual, expected, tolerance = 0)))
  }
  saveRDS(metadata, file.path(config$cache, "metadata.rds"))
  reopened <- readRDS(file.path(config$cache, "metadata.rds"))
  stopifnot(identical(reopened, metadata),
    identical(file.info(config$input)[, c("size", "mtime")], input_info),
    identical(file.info(embedding_input)[, c("size", "mtime")], embedding_info))
  info <- c(info, paste("Cells:", nrow(cells)), paste("Genes:", nrow(genes)),
    paste("UMAP matching: exact barcode; all", nrow(cells), "source cells matched;",
      embedding_cells - nrow(cells), "embedding-only cells not displayed."),
    paste("Cluster counts:", paste(names(table(cells$cluster)),
      as.integer(table(cells$cluster)), sep = "=", collapse = ", ")),
    paste("Translated labels:", sum(!is.na(translated))),
    paste("Expression file bytes:", file.size(expression_file)),
    paste("Elapsed seconds:", round(as.numeric(difftime(Sys.time(), started,
      units = "secs")), 1)),
    "Validation: reopened metadata identical; all saved expression values, cell IDs,",
    "gene IDs and coordinates checked against their sources; RNA/UMAP source sizes/mtimes unchanged.")
  completed <- TRUE
  message("Preparation and saved-data validation passed: ", config$cache)
  invisible(metadata)
}
