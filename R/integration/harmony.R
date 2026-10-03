suppressPackageStartupMessages({
  library(Seurat)
})

run_harmony_integration <- function(seurat_obj, dims, seed = 1234) {
  set.seed(seed)
  # Direct Harmony honors dims.use; Seurat's layer wrapper uses all input PCs.
  harmony::RunHarmony(
    object = seurat_obj,
    group.by.vars = "sample",
    reduction.use = "pca",
    dims.use = dims,
    reduction.save = "harmony",
    project.dim = FALSE,
    max_iter = 10,
    verbose = TRUE
  )
}
