# scDblFinder scores and calls for one physical capture, before QC filtering.
# Sourced by main.R; the parameter values live in
# config$qc$doublets$scdblfinder and are recorded in .INFO.
#
# dbr is a fraction, not a percentage. The "bd_table" value derives the rate
# per sample from the BD Rhapsody multiplet table, the same derivation the
# pipeline used with Scrublet; an explicit NULL uses the package's 10X rate
# prior. dbr_sd = 1 removes the rate-deviation penalty in thresholding.
# Neither is an implicit BD Rhapsody default. Scores are not on Scrublet's
# scale: use scDblFinder.class, never the old Scrublet threshold of 0.15.

# Interpolate BD percentages; extend the last segment above the table without a cap.
bd_multiplet_rate <- function(n_cells, table) {
  rate <- approx(table$cells, table$rate, xout = n_cells, rule = 2)$y
  if (n_cells > max(table$cells)) {
    last <- tail(table, 2L)
    rate <- last$rate[2] + (n_cells - last$cells[2]) *
      diff(last$rate) / diff(last$cells)
  }
  return(rate / 100)
}

# Resolve the configured dbr into the adapter's numeric-or-NULL contract:
# "bd_table" derives the per-sample rate from the BD Rhapsody table, a number
# passes through, NULL keeps the package's 10X prior.
resolve_scdblfinder_dbr <- function(n_cells, doublets_config) {
  dbr <- doublets_config$scdblfinder$dbr
  if (identical(dbr, "bd_table")) {
    rate <- bd_multiplet_rate(n_cells, doublets_config$bd_multiplet_table)
    message(sprintf(
      "scDblFinder expected doublet rate: %.4f%% for %d input cells (BD interpolation; upper linear extrapolation, no cap).",
      100 * rate, n_cells
    ))
    return(rate)
  }
  if (!is.null(dbr) && !(is.numeric(dbr) && length(dbr) == 1L)) {
    stop("config$qc$doublets$scdblfinder$dbr must be \"bd_table\", one fraction, or NULL.")
  }
  return(dbr)
}

detect_scdblfinder_doublets <- function(seurat_obj, dbr, dbr_sd, clusters, seed) {
  # Require explicit scientific choices, including an intentional NULL.
  if (missing(dbr) || missing(dbr_sd) || missing(clusters) || missing(seed)) {
    stop("Specify dbr, dbr_sd, clusters and seed explicitly.")
  }
  stopifnot(
    "dbr must be NULL or a single fraction in [0, 1]" =
      is.null(dbr) || (is.numeric(dbr) && length(dbr) == 1L &&
        is.finite(dbr) && dbr >= 0 && dbr <= 1),
    "dbr_sd must be NULL or a single absolute rate tolerance in [0, 1]" =
      is.null(dbr_sd) || (is.numeric(dbr_sd) && length(dbr_sd) == 1L &&
        is.finite(dbr_sd) && dbr_sd >= 0 && dbr_sd <= 1),
    "clusters must be TRUE or FALSE" =
      is.logical(clusters) && length(clusters) == 1L && !is.na(clusters),
    "seed must be one nonnegative integer representable by R" =
      is.numeric(seed) && length(seed) == 1L && is.finite(seed) &&
        seed >= 0 && seed <= .Machine$integer.max && seed %% 1 == 0
  )

  # Do not join a merged Assay5 or assume biological samples are captures.
  # The caller must confirm that this object's cells belong to one capture.
  samples <- unique(seurat_obj$sample)
  count_layers <- SeuratObject::Layers(seurat_obj[["RNA"]], search = "^counts($|\\.)")
  if (length(samples) != 1L || anyNA(samples) || length(count_layers) != 1L) {
    stop("Supply one unmerged sample from one confirmed physical capture.")
  }
  counts <- SeuratObject::LayerData(
    seurat_obj, assay = "RNA", layer = count_layers[[1]]
  )
  stopifnot(
    "RNA counts must cover all cells in object order" =
      identical(colnames(counts), colnames(seurat_obj)),
    "Cell and gene identifiers must be present and unique" =
      !anyNA(dimnames(counts)[[1]]) && !anyNA(dimnames(counts)[[2]]) &&
        !anyDuplicated(rownames(counts)) && !anyDuplicated(colnames(counts)),
    "Counts must be finite and nonnegative" =
      all(is.finite(counts@x)) && all(counts@x >= 0)
  )
  if (any(Matrix::colSums(counts) == 0)) {
    stop("Zero-count cells require an explicit exclusion decision; none were dropped.")
  }

  # Counts-only SCE avoids transferring normalized assays or diagnostic PCA.
  # Expected (fractional) counts are accepted; do not round SoupX output.
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
  set.seed(seed)
  sce <- scDblFinder::scDblFinder(
    sce, samples = NULL, clusters = clusters, dbr = dbr, dbr.sd = dbr_sd,
    returnType = "sce", BPPARAM = BiocParallel::SerialParam(RNGseed = seed)
  )
  stopifnot(
    "scDblFinder must preserve cell and gene identifiers" =
      identical(dimnames(sce), dimnames(counts))
  )
  threshold <- S4Vectors::metadata(sce)$scDblFinder.threshold
  stopifnot("scDblFinder must return one finite threshold" =
    is.numeric(threshold) && length(threshold) == 1L && is.finite(threshold))
  predictions <- as.data.frame(SummarizedExperiment::colData(sce))[
    colnames(seurat_obj), c("scDblFinder.score", "scDblFinder.class"), drop = FALSE
  ]
  stopifnot(
    "Predictions must cover every input barcode" =
      identical(rownames(predictions), colnames(seurat_obj)) && !anyNA(predictions),
    "Scores must be finite and between 0 and 1" =
      all(is.finite(predictions$scDblFinder.score)) &&
        all(predictions$scDblFinder.score >= 0 & predictions$scDblFinder.score <= 1),
    "Classes must be singlet or doublet" =
      all(predictions$scDblFinder.class %in% c("singlet", "doublet"))
  )

  # Retain native columns as well as the existing generic pipeline contract.
  # Preserve the original Seurat counts, metadata, layers and reductions.
  predictions$doublet_scores <- predictions$scDblFinder.score
  predictions$predicted_doublets <- predictions$scDblFinder.class == "doublet"
  predictions$doublet_filter_method <- "scDblFinder"
  seurat_obj <- SeuratObject::AddMetaData(seurat_obj, predictions)
  seurat_obj@misc$qc$scdblfinder <- list(
    sample = samples[[1]], n_input_cells = ncol(counts),
    n_input_genes = nrow(counts), counts_layer = count_layers[[1]],
    parameters = list(dbr = dbr, dbr.sd = dbr_sd, clusters = clusters, seed = seed),
    samples = NULL, BPPARAM = "SerialParam", returnType = "sce",
    threshold = threshold,
    n_doublets = sum(predictions$predicted_doublets),
    versions = vapply(
      c("scDblFinder", "SingleCellExperiment", "BiocParallel", "scater", "xgboost"),
      function(pkg) as.character(utils::packageVersion(pkg)), character(1)
    )
  )
  message(samples[[1]], " / scDblFinder: ", sum(predictions$predicted_doublets),
    " doublets called among ", ncol(counts), " input cells; none removed.")
  return(seurat_obj)
}
