suppressPackageStartupMessages({
  library(Seurat)
})

run_integration_umaps <- function(
  seurat_obj, dims, scvi_dims, seed = 1234,
  reductions = c(
    umap = "pca", umap_harmony = "harmony",
    umap_cca = "integrated_cca"
  )
) {
  for (name in names(reductions)) {
    seurat_obj <- RunUMAP(
      seurat_obj,
      reduction = reductions[[name]],
      dims = if (name == "umap_scvi") scvi_dims else dims,
      reduction.name = name,
      seed.use = seed,
      verbose = TRUE
    )
  }
  seurat_obj
}

# Shared by integration and clustering; never refit an existing embedding.
run_saved_tsne <- function(obj, reduction, dims, name, seed = 1234, perplexity = 30) {
  settings <- list(reduction = reduction, dims = dims, seed = seed, perplexity = perplexity)
  if (!name %in% Reductions(obj)) {
    obj <- RunTSNE(obj, reduction = reduction, dims = dims,
      reduction.name = name, seed.use = seed, perplexity = perplexity)
    # RunTSNE's command log is overwritten by subsequent t-SNE fits and omits perplexity.
    obj[[name]]@misc$tsne_settings <- settings
  }
  stopifnot("Saved t-SNE must document matching input and settings" =
    identical(obj[[name]]@misc$tsne_settings, settings))
  coordinates <- Embeddings(obj, name)
  stopifnot("t-SNE must contain finite, aligned coordinates for every cell" =
    identical(rownames(coordinates), colnames(obj)) &&
      ncol(coordinates) == 2L && all(is.finite(coordinates)))
  obj
}

run_integration_tsnes <- function(obj, config) {
  for (umap_name in names(config$reductions)) {
    reduction <- config$reductions[[umap_name]]
    obj <- run_saved_tsne(obj, reduction = reduction,
      dims = if (reduction == "integrated_scvi") config$scvi_dims else config$dims,
      name = sub("^umap", "tsne", umap_name), seed = config$seed)
  }
  obj
}
