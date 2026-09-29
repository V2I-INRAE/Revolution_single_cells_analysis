suppressPackageStartupMessages({
  library(Seurat)
})

run_integration_umaps <- function(
  seurat_obj, dims = 1:30, seed = 1234,
  reductions = c(
    umap = "pca", umap_harmony = "harmony",
    umap_cca = "integrated_cca"
  )
) {
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
