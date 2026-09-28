# Doublet detection with DoubletFinder for BD Rhapsody samples.
#   - protocol precedent: Li et al., Current Protocols 2024
#     (doi 10.1002/cpz1.963) - interpolate BD's table at the dataset cell
#     count, feed the rate to DoubletFinder, and LABEL doublets instead of
#     removing them. We label here; removal happens in main.R, in the loop.


suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
})

# --- Linear interpolation of the BD values.
# As such we can use sample specific doublet rate
bd_multiplet_rate <- function(
  n_cells,
  table = qc_params$doublets$bd_multiplet_table
) {
  rate <- approx(
    table$cells,
    table$rate,
    xout = n_cells,
    rule = 2
  )$y
  if (n_cells > max(table$cells) || n_cells < min(table$cells)) {
    warning(
      "Cell count outside BD table range; rate clamped to ",
      min(table$rate), "-", max(table$rate), "%"
    )
  }
  return(rate / 100)
}

# --- Leverages cell annotations to model the proportion of homotypic doublets.
# --- Here we use the annonated cluster provided by the
estimate_homotypic <- function(seurat_obj, params = qc_params$doublets) {
  cluster_col <- paste0("RNA_snn_res.", params$cluster_resolution)
  stopifnot(cluster_col %in% colnames(seurat_obj[[]]))
  clusters <- seurat_obj[[cluster_col]][[1]]
  modelHomotypic(clusters)
}

# --- pK selection from the DoubletFinder parameter sweep
find_optimal_pk <- function(seurat_obj, params = qc_params$doublets) {
  set.seed(params$seed) # paramSweep samples cells to build artificial doublets
  sweep_res <- paramSweep(seurat_obj, PCs = params$pcs, sct = params$sct)
  sweep_stats <- summarizeSweep(sweep_res, GT = params$ground_truth)
  bcmvn <- find.pK(sweep_stats)
  return(bcmvn)
}

#' Detect doublets in one sample and flag them (no removal)
#'
#' Input: Seurat object processed by preprocess_for_doublets().
#' Output: Processed Seurat object with doublet_score (pANN, 0-1) and
#' doublet_class ("Singlet"/"Doublet") metadata.
detect_doublets <- function(
  seurat_obj,
  sample_id,
  params = qc_params$doublets,
  doublet_rate = params$expected_rate
) {

  n_cells <- ncol(seurat_obj)
  rate <- doublet_rate %||% bd_multiplet_rate(
    n_cells,
    params$bd_multiplet_table
  )
  n_exp <- round(rate * n_cells)

  homotypic <- estimate_homotypic(seurat_obj, params)
  n_exp_adj <- round(n_exp * (1 - homotypic))
  n_exp_adj <- max(n_exp_adj, params$min_expected_doublets)

  bcmvn <- find_optimal_pk(seurat_obj, params)
  # pK comes back as a factor: as.character() first or as.numeric() gives
  # the level index, not the pK value
  pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  stopifnot(!is.na(pK), pK > 0, pK <= 1)

  set.seed(params$seed)
  seurat_obj <- doubletFinder(seurat_obj,
    PCs = params$pcs,
    pN = params$pN,
    pK = pK,
    nExp = n_exp_adj,
    sct = params$sct
  )

  pann_col <- grep("^pANN_", colnames(seurat_obj[[]]), value = TRUE)
  class_col <- grep(
    "^DF.classifications_",
    colnames(seurat_obj[[]]),
    value = TRUE
  )
  stopifnot(length(pann_col) == 1, length(class_col) == 1)
  seurat_obj$doublet_score <- seurat_obj[[pann_col]][[1]]
  seurat_obj$doublet_class <- seurat_obj[[class_col]][[1]]

  n_doublets <- sum(seurat_obj$doublet_class == "Doublet")
  cat("\nDoublet detection:", sample_id, "\n")
  cat(sprintf("  cells                 %d\n", n_cells))
  cat(sprintf("  BD multiplet rate     %.2f%%\n", 100 * rate))
  cat(sprintf("  nExp (before adj.)    %d\n", n_exp))
  cat(sprintf("  homotypic proportion  %.3f\n", homotypic))
  cat(sprintf("  nExp (adjusted)       %d\n", n_exp_adj))
  cat(sprintf("  pK                    %s\n", pK))
  cat(sprintf(
    "  called doublets       %d (%.1f%%)\n",
    n_doublets, 100 * mean(seurat_obj$doublet_class == "Doublet")
  ))

  return(seurat_obj)
}

# --- Copy DoubletFinder calls to the original, unprocessed object
add_doublet_calls <- function(
  seurat_obj,
  doublet_obj
) {
  stopifnot(
    setequal(colnames(doublet_obj), colnames(seurat_obj)),
    all(c("doublet_score", "doublet_class") %in% colnames(doublet_obj[[]]))
  )

  doublet_meta <- data.frame(
    doublet_score = doublet_obj$doublet_score,
    doublet_class = doublet_obj$doublet_class,
    row.names = colnames(doublet_obj)
  )
  stopifnot(
    setequal(rownames(doublet_meta), colnames(seurat_obj)),
    anyDuplicated(rownames(doublet_meta)) == 0
  )
  doublet_meta <- doublet_meta[colnames(seurat_obj), , drop = FALSE]
  seurat_obj <- AddMetaData(seurat_obj, doublet_meta)

  return(seurat_obj)
}

# --- Build per-cell DoubletFinder data for the diagnostic plots
collect_doublet_summary <- function(
  seurat_obj,
  sample_id,
  coordinates
) {
  stopifnot(
    all(c("doublet_score", "doublet_class") %in% colnames(seurat_obj[[]])),
    all(c("cell", "UMAP1", "UMAP2") %in% colnames(coordinates)),
    anyDuplicated(coordinates$cell) == 0
  )
  coordinate_idx <- match(colnames(seurat_obj), coordinates$cell)
  stopifnot(!anyNA(coordinate_idx))

  data.frame(
    sample_id = sample_id,
    cell = colnames(seurat_obj),
    UMAP1 = coordinates$UMAP1[coordinate_idx],
    UMAP2 = coordinates$UMAP2[coordinate_idx],
    doublet_score = seurat_obj$doublet_score,
    doublet_class = seurat_obj$doublet_class,
    nFeature_RNA = seurat_obj$nFeature_RNA,
    nCount_RNA = seurat_obj$nCount_RNA,
    percent.mt = seurat_obj$percent.mt,
    row.names = NULL
  )
}
