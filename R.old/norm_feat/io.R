suppressPackageStartupMessages({
  library(Seurat)
})

read_normalization_input <- function(input_file) {
  obj <- readRDS(input_file)
  stopifnot(
    inherits(obj, "Seurat"),
    "RNA" %in% Assays(obj),
    anyDuplicated(colnames(obj)) == 0,
    identical(rownames(obj[[]]), colnames(obj))
  )

  DefaultAssay(obj) <- "RNA"
  rna_layers <- Layers(obj[["RNA"]])
  stopifnot(
    "RNA counts are required" = any(grepl("^counts($|\\.)", rna_layers)),
    "Expected counts-only RNA input" = !any(
      grepl("^(data|scale\\.data)($|\\.)", rna_layers)
    )
  )

  # Retain the input layer structure; do not join or split samples here.
  return(obj)
}
