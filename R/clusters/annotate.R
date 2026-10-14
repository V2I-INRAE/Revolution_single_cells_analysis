# Manual cluster annotation, based on the Seurat skill's annotate_clusters_manual().
# Callers supply the approved labels and own saving the object and run records.

#' Manually annotate clusters with cell type or subset labels
#'
#' @param seurat_obj Seurat object with clusters.
#' @param annotations Named vector: cluster IDs -> labels.
#' @param cluster_col Cluster metadata column (e.g. "harmony_snn_res.0.7").
#'   NULL uses the active identities.
#' @param new_col_name Name for the annotation column (default: "celltype").
#' @return Seurat object with annotations as metadata and active identities.
annotate_clusters_manual <- function(
  seurat_obj,
  annotations,
  cluster_col = NULL,
  new_col_name = "celltype"
) {
  message("Applying manual cell type annotations")

  if (!is.null(cluster_col)) {
    if (!cluster_col %in% colnames(seurat_obj[[]])) {
      stop("Cluster column not found: ", cluster_col)
    }
    current_clusters <- as.character(seurat_obj[[cluster_col, drop = TRUE]])
    message("  Using cluster column: ", cluster_col)
  } else {
    current_clusters <- as.character(SeuratObject::Idents(seurat_obj))
  }

  cell_types <- annotations[current_clusters]

  unmapped <- unique(current_clusters[is.na(cell_types)])
  if (length(unmapped) > 0) {
    warning("Unmapped clusters: ", paste(unmapped, collapse = ", "))
    cell_types[is.na(cell_types)] <- paste0("Cluster_", current_clusters[is.na(cell_types)])
  }

  names(cell_types) <- colnames(seurat_obj)
  seurat_obj[[new_col_name]] <- cell_types
  SeuratObject::Idents(seurat_obj) <- new_col_name

  message("Annotation complete:")
  print(table(seurat_obj[[new_col_name]]))

  return(seurat_obj)
}
