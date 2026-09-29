suppressPackageStartupMessages({
  library(Seurat)
})

run_scvi_integration <- function(
  seurat_obj, dims = 1:30, seed = 1234, max_epochs = NULL
) {
  features <- VariableFeatures(seurat_obj[["RNA"]])
  layers <- Layers(seurat_obj[["RNA"]], search = "^counts($|\\.)")
  stopifnot(
    "scVI requires exactly 3000 LogNormalize RNA HVGs" =
      length(features) == 3000L,
    "Original RNA count layers are required" = length(layers) > 0L
  )
  for (layer in layers) {
    stopifnot(
      "All 3000 HVGs must be present in every RNA count layer" =
        all(features %in% Features(seurat_obj[["RNA"]], layer = layer))
    )
  }

  python <- file.path(getwd(), ".venv-scvi", "bin", "python")
  reticulate::use_python(python, required = TRUE)
  scvi <- reticulate::import("scvi", convert = FALSE)
  set.seed(seed)
  scvi$settings$seed <- as.integer(seed)

  # Select LogNormalize HVGs, but train on original RNA counts, not data.
  # The wrapper derives batches before joining counts; it does not use PCA.
  seurat_obj <- IntegrateLayers(
    object = seurat_obj,
    method = SeuratWrappers::scVIIntegration,
    assay = "RNA",
    features = features,
    layers = layers,
    scale.layer = NULL,
    orig.reduction = NULL,
    conda_env = python,
    new.reduction = "integrated_scvi",
    ndims = length(dims),
    nlayers = 2L,
    gene_likelihood = "nb",
    max_epochs = max_epochs
  )
  # NumPy cell names can retain a 1D array attribute through the wrapper.
  rownames(seurat_obj[["integrated_scvi"]]@cell.embeddings) <-
    as.character(Cells(seurat_obj[["integrated_scvi"]]))
  DefaultAssay(seurat_obj[["integrated_scvi"]]) <- "RNA"
  seurat_obj@misc$scvi <- list(
    count_assay = "RNA", count_layers = layers,
    feature_assay = "RNA", features = features,
    batch_col = "sample", ndims = length(dims), seed = seed,
    nlayers = 2L, gene_likelihood = "nb", max_epochs = max_epochs,
    python = python,
    scvi_version = reticulate::py_to_r(scvi$`__version__`)
  )
  seurat_obj
}
