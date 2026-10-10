suppressPackageStartupMessages(library(Seurat))

read_marker_input <- function(source_file) {
  obj <- readRDS(source_file)
  DefaultAssay(obj) <- "RNA"
  layers <- Layers(obj[["RNA"]], search = "^data($|\\.)")
  # Keep every cell and gene; only discard components unused by marker discovery.
  obj <- DietSeurat(
    obj,
    assays = "RNA",
    layers = layers,
    dimreducs = NULL,
    graphs = NULL,
    misc = FALSE
  )

  JoinLayers(obj, assay = "RNA", layers = "data")
}
