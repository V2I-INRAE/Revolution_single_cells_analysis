# LogNormalize, variable-feature selection and scaling.

suppressPackageStartupMessages({
  library(Seurat)
})

run_lognormalize <- function(seurat_obj, norm_config) {
  message("Running LogNormalize normalization")
  message("  Method: ", norm_config$lognorm_method)
  message("  Scale factor: ", norm_config$lognorm_scale_factor)
  seurat_obj <- NormalizeData(
    seurat_obj,
    normalization.method = norm_config$lognorm_method,
    scale.factor = norm_config$lognorm_scale_factor,
    verbose = TRUE
  )
  message("Normalization complete")
  seurat_obj
}

find_hvgs <- function(seurat_obj, norm_config) {
  message("Finding highly variable genes")
  message("  Method: ", norm_config$hvg_selection_method)
  message("  Number of features: ", norm_config$n_variable_features)
  seurat_obj <- FindVariableFeatures(
    seurat_obj,
    selection.method = norm_config$hvg_selection_method,
    nfeatures = norm_config$n_variable_features,
    verbose = TRUE
  )
  hvgs <- VariableFeatures(seurat_obj)
  message("  Variable features identified: ", length(hvgs))
  message("  Top 10 HVGs: ", paste(head(hvgs, 10), collapse = ", "))
  seurat_obj
}

scale_data <- function(seurat_obj, norm_config) {
  message("Scaling data")
  features <- VariableFeatures(seurat_obj)
  message("  Scaling specified features: ", length(features))
  message("  Variables to regress: ",
    paste(norm_config$vars_to_regress, collapse = ", "))
  missing_vars <- setdiff(norm_config$vars_to_regress, colnames(seurat_obj@meta.data))
  if (length(missing_vars) > 0L) {
    stop("Variables not found in metadata: ",
      paste(missing_vars, collapse = ", "))
  }
  seurat_obj <- ScaleData(
    seurat_obj,
    features = features,
    vars.to.regress = norm_config$vars_to_regress,
    verbose = TRUE
  )
  message("Scaling complete")
  seurat_obj
}
