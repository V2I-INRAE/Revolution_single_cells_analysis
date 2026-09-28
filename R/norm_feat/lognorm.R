suppressPackageStartupMessages({
  library(Seurat)
})

run_lognormalize <- function(
  seurat_obj,
  normalization_method = "LogNormalize",
  scale_factor = 10000,
  verbose = TRUE
) {
  message("Running LogNormalize normalization")
  message("  Method: ", normalization_method)
  message("  Scale factor: ", scale_factor)

  # Normalize data
  seurat_obj <- NormalizeData(
    seurat_obj,
    normalization.method = normalization_method,
    scale.factor = scale_factor,
    verbose = verbose
  )

  message("Normalization complete")

  return(seurat_obj)
}

find_hvgs <- function(
  seurat_obj,
  selection_method = "vst",
  n_features = 2000,
  verbose = TRUE
) {
  message("Finding highly variable genes")
  message("  Method: ", selection_method)
  message("  Number of features: ", n_features)

  # Find variable features
  seurat_obj <- FindVariableFeatures(
    seurat_obj,
    selection.method = selection_method,
    nfeatures = n_features,
    verbose = verbose
  )

  # Get variable features
  hvgs <- VariableFeatures(seurat_obj)

  message("  Variable features identified: ", length(hvgs))
  message("  Top 10 HVGs: ", paste(head(hvgs, 10), collapse = ", "))

  return(seurat_obj)
}

scale_data <- function(
  seurat_obj,
  features = NULL,
  vars_to_regress = NULL,
  verbose = TRUE
) {
  message("Scaling data")

  # Use all genes if not specified
  if (is.null(features)) {
    features <- rownames(seurat_obj)
    message("  Scaling all genes: ", length(features))
  } else {
    message("  Scaling specified features: ", length(features))
  }

  # Check if vars_to_regress exist
  if (!is.null(vars_to_regress)) {
    message("  Variables to regress: ", paste(vars_to_regress, collapse = ", "))
    missing_vars <- setdiff(vars_to_regress, colnames(seurat_obj@meta.data))
    if (length(missing_vars) > 0) {
      warning(
        "Variables not found in metadata: ",
        paste(missing_vars, collapse = ", ")
      )
      vars_to_regress <- intersect(
        vars_to_regress, colnames(seurat_obj@meta.data)
      )
    }
  }

  # Scale data
  seurat_obj <- ScaleData(
    seurat_obj,
    features = features,
    vars.to.regress = vars_to_regress,
    verbose = verbose
  )

  message("Scaling complete")

  return(seurat_obj)
}
