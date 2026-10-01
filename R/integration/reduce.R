suppressPackageStartupMessages({
  library(Seurat)
})

run_integration_umaps <- function(
  seurat_obj, dims = 1:30, scvi_dims = dims, seed = 1234,
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
