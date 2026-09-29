suppressPackageStartupMessages({
  library(Seurat)
})

run_sctransform <- function(
  seurat_obj,
  vars_to_regress,
  n_genes,
  feature_buffer,
  verbose = TRUE
) {

  message("Running SCTransform normalization")
  message("  Variables to regress: ", paste(vars_to_regress, collapse = ", "))
  message("  Variable genes to identify: ", n_genes)

  # Check if vars_to_regress exist in metadata
  missing_vars <- setdiff(vars_to_regress, colnames(seurat_obj@meta.data))
  if (length(missing_vars) > 0) {
    warning(
      "Variables not found in metadata: ", paste(missing_vars, collapse = ", ")
    )
    vars_to_regress <- intersect(
      vars_to_regress, colnames(seurat_obj@meta.data)
    )
  }

  seurat_obj <- SCTransform(
    seurat_obj,
    vars.to.regress = vars_to_regress,
    variable.features.n = n_genes + feature_buffer,
    verbose = verbose
  )

  residual_features <- rownames(SeuratObject::LayerData(
    seurat_obj,
    assay = "SCT",
    layer = "scale.data"
  ))
  variable_features <- head(
    VariableFeatures(seurat_obj)[
      VariableFeatures(seurat_obj) %in% residual_features
    ],
    n_genes
  )
  VariableFeatures(seurat_obj) <- variable_features
  SeuratObject::LayerData(
    seurat_obj, assay = "SCT", layer = "scale.data"
  ) <- SeuratObject::LayerData(
    seurat_obj, assay = "SCT", layer = "scale.data"
  )[variable_features, , drop = FALSE]

  message("SCTransform complete")
  message("  Default assay: ", DefaultAssay(seurat_obj))
  message(
    "  Variable features with stored residuals: ",
    length(variable_features)
  )

  seurat_obj
}
