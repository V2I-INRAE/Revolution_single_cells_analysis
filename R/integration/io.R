suppressPackageStartupMessages({
  library(Seurat)
})

read_integration_input <- function(input_file, dims = 1:30) {
  obj <- readRDS(input_file)
  setup_integration(obj, dims = dims)
}

read_integration_comparison <- function(scvi_file, harmony_file) {
  obj <- readRDS(scvi_file)
  stopifnot("Expected a LogNormalize scVI checkpoint" =
    identical(obj@misc$integration$method, "scvi") &&
      identical(obj@misc$integration$route, "lognorm"))

  obj <- DietSeurat(obj, assays = "RNA", layers = "counts",
    dimreducs = c("pca", "umap", "umap_scvi"), graphs = NULL)
  harmony <- readRDS(harmony_file)
  stopifnot("Expected a LogNormalize Harmony checkpoint" =
    identical(harmony@misc$integration$method, "harmony") &&
      identical(harmony@misc$integration$route, "lognorm"))
  cells <- colnames(obj)
  stopifnot("Integration checkpoints must contain the same unique cells" =
    !anyDuplicated(cells) && !anyDuplicated(colnames(harmony)) &&
      setequal(cells, colnames(harmony)))
  for (field in c("sample", "pressure", "time_point")) {
    stopifnot("Integration checkpoints require complete grouping metadata" =
      field %in% colnames(obj[[]]) && field %in% colnames(harmony[[]]) &&
        !anyNA(obj[[field]]) && !anyNA(harmony[[field]]))
    stopifnot("Integration checkpoint metadata disagree" =
      identical(as.character(obj[[]][cells, field]),
        as.character(harmony[[]][cells, field])))
  }
  for (reduction in c("pca", "umap")) {
    stopifnot("Integration checkpoints have different baseline coordinates" =
      identical(Embeddings(obj, reduction)[cells, , drop = FALSE],
        Embeddings(harmony, reduction)[cells, , drop = FALSE]))
  }
  stopifnot("Integration checkpoints have different PCA features or loadings" =
    identical(VariableFeatures(obj), VariableFeatures(harmony)) &&
      identical(Loadings(obj, "pca"), Loadings(harmony, "pca")))
  obj[["umap_harmony"]] <- harmony[["umap_harmony"]]
  obj@misc$integration_comparison <- data.frame(
    method = c("scVI", "Harmony"),
    input_file = normalizePath(c(scvi_file, harmony_file), mustWork = TRUE),
    run_id = c(obj@misc$integration$run_id, harmony@misc$integration$run_id)
  )
  obj
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

  seurat_obj
}
