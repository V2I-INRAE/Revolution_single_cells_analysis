suppressPackageStartupMessages({
  library(Seurat)
})

run_pca_analysis <- function(
  seurat_obj,
  features = NULL,
  n_pcs = 50,
  verbose = TRUE
) {
  message("Running PCA analysis")

  if (is.null(features)) {
    features <- VariableFeatures(seurat_obj)
    message("  Using variable features: ", length(features))
  } else {
    message("  Using specified features: ", length(features))
  }

  seurat_obj <- RunPCA(
    seurat_obj,
    features = features,
    npcs = n_pcs,
    verbose = verbose
  )

  message("PCA complete")
  message("  PCs computed: ", n_pcs)

  n_features_used <- nrow(Loadings(seurat_obj, reduction = "pca"))
  n_var_features <- length(VariableFeatures(seurat_obj))
  if (n_features_used == n_var_features) {
    message(
      "  PCA loadings verified: ",
      n_features_used,
      " variable features used"
    )
  } else {
    message(
      "  WARNING: PCA used ", n_features_used, " features but ",
      n_var_features, " variable features exist. Check feature selection."
    )
  }

  pca_sdev <- Stdev(seurat_obj, reduction = "pca")
  cumulative <- cumsum(pca_sdev^2 / sum(pca_sdev^2))
  message(sprintf("  PC1-10 explain %.1f%% of variance", 100 * cumulative[10]))
  message(sprintf("  PC1-20 explain %.1f%% of variance", 100 * cumulative[20]))
  message(sprintf("  PC1-30 explain %.1f%% of variance", 100 * cumulative[30]))

  if (verbose) {
    message("\nTop features for PC1-5:")
    print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)
  }

  seurat_obj
}
