suppressPackageStartupMessages({
  library(Seurat)
})

read_integration_input <- function(input_file, dims = 1:30) {
  obj <- readRDS(input_file)
  setup_integration(obj, dims = dims)
}

setup_integration <- function(seurat_obj, dims = 1:30) {
  metadata <- seurat_obj[[]]
  assay <- seurat_obj[[DefaultAssay(seurat_obj)]]
  pca <- Embeddings(seurat_obj, "pca")
  features <- VariableFeatures(assay)
  stopifnot(
    "Expected dimensions 1 through the selected PC count" =
      identical(as.integer(dims), seq_len(max(dims))),
    anyDuplicated(colnames(seurat_obj)) == 0L,
    identical(rownames(metadata), colnames(seurat_obj)),
    identical(rownames(pca), colnames(seurat_obj)),
    max(dims) <= ncol(pca),
    all(is.finite(pca[, dims, drop = FALSE])),
    DefaultAssay(seurat_obj[["pca"]]) == DefaultAssay(seurat_obj),
    "sample" %in% names(metadata),
    !anyNA(metadata$sample),
    all(metadata$sample != ""),
    length(unique(metadata$sample)) > 1L,
    length(features) > 0L,
    all(features %in% rownames(Loadings(seurat_obj, "pca"))),
    all(features %in% rownames(SeuratObject::LayerData(
      assay, layer = "scale.data"
    ))),
    length(SeuratObject::Layers(
      seurat_obj[["RNA"]], search = "^counts($|\\.)"
    )) > 0L
  )

  if (inherits(assay, "SCTAssay")) {
    groups <- lapply(levels(assay), function(model) {
      stopifnot(all(features %in% rownames(
        SCTResults(assay, slot = "feature.attributes", model = model)
      )))
      rownames(SCTResults(assay, slot = "cell.attributes", model = model))
    })
  } else {
    stopifnot(inherits(assay, "Assay5"))
    layers <- SeuratObject::Layers(assay, search = "^data($|\\.)")
    groups <- lapply(layers, function(layer) {
      stopifnot(all(features %in% SeuratObject::Features(assay, layer = layer)))
      Cells(assay, layer = layer)
    })
  }
  cells <- unlist(groups, use.names = FALSE)
  stopifnot(
    "Layers/models must partition all cells" =
      !anyDuplicated(cells) && setequal(cells, colnames(seurat_obj))
  )
  samples <- vapply(groups, function(ids) {
    sample <- unique(as.character(metadata[ids, "sample"]))
    stopifnot("Each layer/model must contain one sample" = length(sample) == 1L)
    sample
  }, character(1))
  stopifnot(
    "Expected one layer/model per sample" =
      !anyDuplicated(samples) && setequal(samples, metadata$sample)
  )

  message(
    "Integration input: ", ncol(seurat_obj), " cells, ",
    length(samples), " samples, assay ", DefaultAssay(seurat_obj)
  )
  # Keep sample layers/models, normalized values, selected features and PCA.
  seurat_obj
}
