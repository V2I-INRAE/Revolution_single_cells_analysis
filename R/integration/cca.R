suppressPackageStartupMessages({
  library(Seurat)
})

run_cca_integration <- function(seurat_obj, dims = 1:30, seed = 1234) {
  normalization <- if (inherits(
    seurat_obj[[DefaultAssay(seurat_obj)]], "SCTAssay"
  )) {
    "SCT"
  } else {
    "LogNormalize"
  }
  set.seed(seed)

  IntegrateLayers(
    object = seurat_obj,
    method = CCAIntegration,
    assay = DefaultAssay(seurat_obj),
    normalization.method = normalization,
    features = VariableFeatures(seurat_obj),
    orig.reduction = "pca",
    new.reduction = "integrated_cca",
    dims = dims,
    dims.to.integrate = dims,
    verbose = TRUE
  )
}
