# Doublet detection with DoubletFinder for BD Rhapsody samples.
#
# Method (run per sample, on the Seurat object already filtered of QC outliers):
#   - preliminary embedding (temporary, stripped afterwards):
#     LogNormalize -> HVG (2000) -> ScaleData -> PCA (1:20) -> neighbors ->
#     clusters (Louvain, resolution 0.6)
#   - expected doublet rate: BD Rhapsody microwell multiplet table
#     (BD Rhapsody Instrument User Guide, Doc ID 214062, Poisson-based,
#     keyed on captured cells): 10,000 -> 2.4%, 12,000 -> 2.8%,
#     15,000 -> 3.5%, 17,000 -> 4.0%. The rate is interpolated at the
#     sample's cell count (bd_multiplet_rate). Note: this information is
#     NOT present in the pipeline reports (no "Multiplet" fields).
#   - protocol precedent: Li et al., Current Protocols 2024
#     (doi 10.1002/cpz1.963) - interpolate BD's table at the dataset cell
#     count, feed the rate to DoubletFinder, and LABEL doublets instead of
#     removing them. We label here; removal happens in main.R, in the loop.
#   - homotypic adjustment from the preliminary clusters (no external
#     annotations): prop = sum((cluster size / N)^2); nExp is adjusted by
#     (1 - prop)
#   - DoubletFinder 2.0.6: paramSweep -> summarizeSweep -> find.pK
#     (max BCmetric) -> doubletFinder(pN = 0.25)
#   - output: doublet_score (pANN) and doublet_class metadata on a
#     counts-only object, plus a per-cell summary table for the
#     cross-sample comparison plots (plot_doublet_comparison).

suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
  library(ggplot2)
})

# BD Rhapsody estimated multiplet rate table (Instrument User Guide,
# Doc ID 214062): captured cells -> multiplet rate (%)
bd_multiplet_table <- data.frame(
  cells = c(
    100,
    500,
    1000,
    2000,
    3000,
    4000,
    5000,
    6000,
    7000,
    8000,
    9000,
    10000,
    11000,
    12000,
    13000,
    14000,
    15000,
    16000,
    17000
  ),
  rate = c(
    0.0,
    0.1,
    0.2,
    0.5,
    0.7,
    1.0,
    1.2,
    1.4,
    1.7,
    1.9,
    2.1,
    2.4,
    2.6,
    2.8,
    3.1,
    3.3,
    3.5,
    3.8,
    4.0
  )
)

# --- Linear interpolation of the BD values.
# As such we can use sample specific doublet rate
bd_multiplet_rate <- function(n_cells) {
  rate <- approx(
    bd_multiplet_table$cells,
    bd_multiplet_table$rate,
    xout = n_cells,
    rule = 2
  )$y
  if (n_cells > max(bd_multiplet_table$cells) ||
    n_cells < min(bd_multiplet_table$cells)) {
    warning(
      "Cell count outside BD table range; rate clamped to ",
      min(bd_multiplet_table$rate), "-", max(bd_multiplet_table$rate), "%"
    )
  }
  return(rate / 100)
}

# --- Prepare the data for doublets identification
preprocess_for_doublets <- function(
  seurat_obj,
  nvar = 2000,
  pcs = 1:20,
  resolution = 0.6,
  seed = 1234
) {
  set.seed(seed)
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj,
    nfeatures = nvar,
    verbose = FALSE
  )
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, npcs = max(pcs), verbose = FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, dims = pcs, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj,
    resolution = resolution,
    algorithm = 1,
    random.seed = seed,
    verbose = FALSE
  )
  return(seurat_obj)
}


# --- Leverages cell annotations to model the proportion of homotypic doublets.
# --- Here we use the annonated cluster provided by the
estimate_homotypic <- function(seurat_obj, resolution = 0.6) {
  cluster_col <- paste0("RNA_snn_res.", resolution)
  stopifnot(cluster_col %in% colnames(seurat_obj[[]]))
  clusters <- seurat_obj[[cluster_col]][[1]]
  modelHomotypic(clusters)
}

# --- pK selection from the DoubletFinder parameter sweep
find_optimal_pk <- function(seurat_obj, pcs = 1:20, seed = 1234) {
  set.seed(seed) # paramSweep samples cells to build artificial doublets
  sweep_res <- paramSweep(seurat_obj, PCs = pcs, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  bcmvn <- find.pK(sweep_stats)
  return(bcmvn)
}

#' Detect doublets in one sample and flag them (no removal)
#'
#' Input : Seurat object filtered of QC outliers.
#' Output: list(
#' obj = counts-only Seurat object with two metadata columns:
#' doublet_score (pANN, 0-1) and doublet_class ("Singlet"/"Doublet");
#' summary = one row per cell (sample, UMAP coordinates, score, class, QC metrics)
#' for plot_doublet_comparison).
#' Doublets are flagged here, the removal happens in main.R
detect_doublets <- function(
  seurat_obj,
  sample_id,
  pcs = 1:20,
  pN = 0.25,
  resolution = 0.6, doublet_rate = NULL,
  seed = 1234
) {

  n_cells <- ncol(seurat_obj)
  rate <- doublet_rate %||% bd_multiplet_rate(n_cells)
  n_exp <- round(rate * n_cells)

  obj <- preprocess_for_doublets(seurat_obj,
    pcs = pcs,
    resolution = resolution, seed = seed
  )

  homotypic <- estimate_homotypic(obj, resolution)
  n_exp_adj <- round(n_exp * (1 - homotypic))
  n_exp_adj <- max(n_exp_adj, 1)

  bcmvn <- find_optimal_pk(obj, pcs = pcs, seed = seed)
  # pK comes back as a factor: as.character() first or as.numeric() gives
  # the level index, not the pK value
  pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  stopifnot(!is.na(pK), pK > 0, pK <= 1)

  set.seed(seed)
  obj <- doubletFinder(obj,
    PCs = pcs, pN = pN, pK = pK, nExp = n_exp_adj,
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
  obj <- RunUMAP(obj, dims = pcs, seed.use = seed, verbose = FALSE)

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

  # --- We remove the temporary embedding and just keep counts
  # and doublet metadata only for filtering
  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  clean_obj <- CreateSeuratObject(
    counts = counts,
    project = paste0(sample_id, "_doublets")
  )
  doublet_meta <- data.frame(
    doublet_score = obj$doublet_score,
    doublet_class = obj$doublet_class,
    row.names = colnames(obj)
  )
  clean_obj <- AddMetaData(clean_obj, doublet_meta)

  return(list(obj = clean_obj, summary = summary))
}
