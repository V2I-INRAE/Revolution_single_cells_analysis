suppressPackageStartupMessages({
  library(Seurat)
})

run_scvi_integration <- function(seurat_obj, params) {
  features <- VariableFeatures(seurat_obj[["RNA"]])
  layers <- Layers(seurat_obj[["RNA"]], search = "^counts($|\\.)")
  stopifnot(
    "scVI requires the selected number of LogNormalize RNA HVGs" =
      length(features) == params$n_features,
    "Original RNA count layers are required" = length(layers) > 0L
  )
  for (layer in layers) {
    stopifnot(
      "All selected HVGs must be present in every RNA count layer" =
        all(features %in% Features(seurat_obj[["RNA"]], layer = layer))
    )
  }

  python <- file.path(getwd(), ".venv-scvi", "bin", "python")
  reticulate::use_python(python, required = TRUE)
  scvi <- reticulate::import("scvi", convert = FALSE)
  set.seed(params$seed)
  scvi$settings$seed <- as.integer(params$seed)

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
    ndims = params$n_latent,
    nlayers = params$n_layers,
    gene_likelihood = params$gene_likelihood,
    max_epochs = params$max_epochs
  )

  rownames(seurat_obj[["integrated_scvi"]]@cell.embeddings) <-
    as.character(Cells(seurat_obj[["integrated_scvi"]]))
  DefaultAssay(seurat_obj[["integrated_scvi"]]) <- "RNA"
  seurat_obj@misc$scvi <- list(
    count_assay = "RNA", count_layers = layers,
    feature_assay = "RNA", features = features,
    batch_col = "sample", ndims = params$n_latent, seed = params$seed,
    nlayers = params$n_layers, gene_likelihood = params$gene_likelihood,
    max_epochs = params$max_epochs,
    python = python,
    scvi_version = reticulate::py_to_r(scvi$`__version__`)
  )
  seurat_obj
}
