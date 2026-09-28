# Doublet detection with DoubletFinder for BD Rhapsody samples.
#   - protocol precedent: Li et al., Current Protocols 2024
#     (doi 10.1002/cpz1.963) - interpolate BD's table at the dataset cell
#     count, feed the rate to DoubletFinder, and LABEL doublets instead of
#     removing them. We label here; removal happens in main.R, in the loop.


suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
  library(ggplot2)
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

# --- Prepare the data for doublets identification
preprocess_for_doublets <- function(
  seurat_obj,
  params = qc_params$doublets
) {
  set.seed(params$seed)
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj,
    nfeatures = params$n_variable_features,
    verbose = FALSE
  )
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(
    seurat_obj,
    npcs = max(params$pcs),
    verbose = FALSE
  )
  seurat_obj <- FindNeighbors(
    seurat_obj,
    dims = params$pcs,
    verbose = FALSE
  )
  seurat_obj <- FindClusters(seurat_obj,
    resolution = params$cluster_resolution,
    algorithm = params$cluster_algorithm,
    random.seed = params$seed,
    verbose = FALSE
  )
  return(seurat_obj)
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
#' Input : Original Seurat object with QC outliers labeled.
#' Output: list(
#' obj = original Seurat object with two additional metadata columns:
#' doublet_score (pANN, 0-1) and doublet_class ("Singlet"/"Doublet");
#' summary = one row per cell (sample, UMAP coordinates, score, class, QC metrics)
#' for plot_doublet_comparison).
#' Doublets are flagged here, the removal happens in main.R
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

  obj <- preprocess_for_doublets(seurat_obj, params)

  homotypic <- estimate_homotypic(obj, params)
  n_exp_adj <- round(n_exp * (1 - homotypic))
  n_exp_adj <- max(n_exp_adj, params$min_expected_doublets)

  bcmvn <- find_optimal_pk(obj, params)
  # pK comes back as a factor: as.character() first or as.numeric() gives
  # the level index, not the pK value
  pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  stopifnot(!is.na(pK), pK > 0, pK <= 1)

  set.seed(params$seed)
  obj <- doubletFinder(obj,
    PCs = params$pcs,
    pN = params$pN,
    pK = pK,
    nExp = n_exp_adj,
    sct = params$sct
  )

  pann_col <- grep("^pANN_", colnames(obj[[]]), value = TRUE)
  class_col <- grep("^DF.classifications_", colnames(obj[[]]), value = TRUE)
  stopifnot(length(pann_col) == 1, length(class_col) == 1)
  obj$doublet_score <- obj[[pann_col]][[1]]
  obj$doublet_class <- obj[[class_col]][[1]]

  n_doublets <- sum(obj$doublet_class == "Doublet")
  cat("\nDoublet detection:", sample_id, "\n")
  cat(sprintf("  cells                 %d\n", n_cells))
  cat(sprintf("  BD multiplet rate     %.2f%%\n", 100 * rate))
  cat(sprintf("  nExp (before adj.)    %d\n", n_exp))
  cat(sprintf("  homotypic proportion  %.3f\n", homotypic))
  cat(sprintf("  nExp (adjusted)       %d\n", n_exp_adj))
  cat(sprintf("  pK                    %s\n", pK))
  cat(sprintf(
    "  called doublets       %d (%.1f%%)\n",
    n_doublets, 100 * mean(obj$doublet_class == "Doublet")
  ))

  # UMAP on the temporary embedding, used only for the comparison plots
  obj <- RunUMAP(
    obj,
    dims = params$pcs,
    seed.use = params$seed,
    verbose = FALSE
  )

  # One row per cell, for the cross-sample comparison plots
  umap <- Embeddings(obj, "umap")
  summary <- data.frame(
    sample_id = sample_id,
    cell = colnames(obj),
    UMAP1 = umap[, 1],
    UMAP2 = umap[, 2],
    doublet_score = obj$doublet_score,
    doublet_class = obj$doublet_class,
    nFeature_RNA = obj$nFeature_RNA,
    nCount_RNA = obj$nCount_RNA,
    percent.mt = obj$percent.mt,
    row.names = NULL
  )

  # --- Join the classifications to the original object by cell barcode.
  # The temporary normalization, embedding and clustering remain only in obj.
  doublet_meta <- data.frame(
    doublet_score = obj$doublet_score,
    doublet_class = obj$doublet_class,
    row.names = colnames(obj)
  )
  stopifnot(
    setequal(rownames(doublet_meta), colnames(seurat_obj)),
    anyDuplicated(rownames(doublet_meta)) == 0
  )
  doublet_meta <- doublet_meta[colnames(seurat_obj), , drop = FALSE]
  seurat_obj <- AddMetaData(seurat_obj, doublet_meta)

  return(list(obj = seurat_obj, summary = summary))
}
