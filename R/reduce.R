# Dimensional reduction for the whole pipeline: temporary QC diagnostics,
# production PCA, integration UMAPs and clustering UMAP reuse.

suppressPackageStartupMessages({
  library(Seurat)
})

# --- QC diagnostics -----------------------------------------------------------

# Temporary per-sample PCA used only for QC visualization and Scrublet;
# never alters the saved counts.
preprocess_qc_umap <- function(seurat_obj, params) {
  set.seed(params$seed)
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj,
    nfeatures = params$n_variable_features,
    verbose = FALSE
  )
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, npcs = max(params$pcs), verbose = FALSE)
  seurat_obj
}

# UMAP coordinates used only by the QC diagnostic plots; one per sample,
# never pooled across samples.
compute_diagnostic_umap <- function(seurat_obj, params) {
  seurat_obj <- RunUMAP(seurat_obj,
    reduction = "pca",
    dims = params$pcs,
    seed.use = params$seed,
    verbose = FALSE
  )
  umap <- Embeddings(seurat_obj, "umap")
  data.frame(
    qc_umap_1 = umap[, 1],
    qc_umap_2 = umap[, 2],
    row.names = rownames(umap)
  )
}

# --- Production PCA -------------------------------------------------------------

run_pca_analysis <- function(seurat_obj, n_pcs) {
  message("Running PCA analysis")
  seurat_obj <- RunPCA(
    seurat_obj,
    features = VariableFeatures(seurat_obj),
    npcs = n_pcs,
    verbose = TRUE
  )

  n_features_used <- nrow(Loadings(seurat_obj, reduction = "pca"))
  n_var_features <- length(VariableFeatures(seurat_obj))
  if (n_features_used == n_var_features) {
    message("  PCA loadings verified: ", n_features_used, " variable features used")
  } else {
    message("  WARNING: PCA used ", n_features_used, " features but ",
      n_var_features, " variable features exist. Check feature selection.")
  }

  pca_sdev <- Stdev(seurat_obj, reduction = "pca")
  cumulative <- cumsum(pca_sdev^2 / sum(pca_sdev^2))
  for (i in c(10, 20, 30)) {
    if (i <= length(cumulative)) {
      message(sprintf("  PC1-%d explain %.1f%% of variance", i, 100 * cumulative[i]))
    }
  }
  seurat_obj
}

# --- Integration UMAPs ------------------------------------------------------------

# Unintegrated PCA UMAP beside the integrated-method UMAP, both on the shared
# dims. `reductions` maps UMAP name to source reduction, e.g.
# c(umap = "pca", umap_harmony = "harmony").
run_integration_umaps <- function(seurat_obj, dims, seed, reductions) {
  for (name in names(reductions)) {
    seurat_obj <- RunUMAP(
      seurat_obj,
      reduction = reductions[[name]],
      dims = dims,
      reduction.name = name,
      seed.use = seed,
      verbose = TRUE
    )
  }
  seurat_obj
}

# --- Clustering UMAP ----------------------------------------------------------------

# The clustering stage reuses the UMAP fitted during integration on the same
# shared dims; verify it is present, documented and aligned instead of refitting.
validate_clustering_umap <- function(seurat_obj, umap_name, reduction, dims) {
  stopifnot(
    "Integration UMAP is required for clustering plots" =
      umap_name %in% Reductions(seurat_obj)
  )
  coordinates <- Embeddings(seurat_obj, umap_name)
  stopifnot(
    "Clustering UMAP must be finite, two-dimensional and cell-aligned" =
      identical(rownames(coordinates), colnames(seurat_obj)) &&
        ncol(coordinates) == 2L &&
        all(is.finite(coordinates))
  )
  seurat_obj@misc$clustering$umap <- umap_name
  seurat_obj@misc$clustering$umap_input_reduction <- reduction
  seurat_obj@misc$clustering$umap_dims <- dims
  seurat_obj
}
