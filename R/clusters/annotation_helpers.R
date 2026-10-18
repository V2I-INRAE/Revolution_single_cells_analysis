# Small in-memory operations on an existing Seurat object.
# These functions neither detect doublets nor infer cell identities. They do
# not save files, record provenance, or recompute/remove fitted analyses.
# The caller selects exact barcodes and handles those responsibilities.

# Remove previously identified doublets, preserving the surviving cell order.
remove_doublets <- function(object, cells) {
  if (!length(cells)) return(object)
  if (!all(cells %in% colnames(object))) stop("Some doublet barcodes are absent from the object.")
  retained <- colnames(object)[!colnames(object) %in% cells]
  subset(object, cells = retained)
}

# Move selected cells between operational compartments stored in $subset.
# source/destination are compartment labels, NOT numeric cluster identifiers.
reassign_cells <- function(object, cells, source, destination) {
  if (!length(cells)) return(object)
  if (length(source) != 1L || length(destination) != 1L ||
      anyNA(c(source, destination)) || any(!nzchar(c(source, destination)))) {
    stop("source and destination must each be one nonempty compartment label.")
  }
  at <- match(cells, colnames(object))
  if (anyNA(at)) stop("Some reassignment barcodes are absent from the object.")
  labels <- as.character(object$subset)
  if (anyNA(labels[at]) || any(labels[at] != source)) {
    stop("Selected cells do not all belong to the source compartment: ", source)
  }
  labels[at] <- destination
  names(labels) <- colnames(object)
  object$subset <- labels # character replacement avoids old factor-level NAs
  SeuratObject::Idents(object) <- "subset"
  object
}
