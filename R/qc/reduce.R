# Dimensional reduction used by the QC pipeline

suppressPackageStartupMessages({
  library(Seurat)
})

# --- Prepare the data for doublet identification
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

# --- Compute UMAP coordinates used only by the QC diagnostic plots
compute_diagnostic_umap <- function(
  seurat_obj,
  params = qc_params$doublets
) {
  seurat_obj <- RunUMAP(
    seurat_obj,
    reduction = "pca",
    dims = params$pcs,
    seed.use = params$seed,
    verbose = FALSE
  )
  umap <- Embeddings(seurat_obj, "umap")

  data.frame(
    cell = rownames(umap),
    UMAP1 = umap[, 1],
    UMAP2 = umap[, 2],
    row.names = NULL
  )
}
