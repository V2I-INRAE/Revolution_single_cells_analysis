# Manual cluster annotation, based on the Seurat skill's annotate_clusters_manual().
# Source this file to reuse the function, or run it to apply coarse_subsets_44542861:
# Rscript R/clusters/annotate.R <clustered.rds> <annotated_clustered.rds>

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

# Run only when called with Rscript, not when sourced by another script.
if (sys.nframe() == 0L) {
  suppressPackageStartupMessages(library(Seurat))
  source("R/clusters/annotations.R")

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) {
    stop("Usage: Rscript R/clusters/annotate.R <clustered.rds> <annotated_clustered.rds>")
  }
  input <- args[[1L]]
  output <- args[[2L]]
  cluster_col <- "harmony_snn_res.0.7"
  input_info <- file.info(input)[, c("size", "mtime")]

  message("Loading (input is never written): ", input)
  obj <- readRDS(input)
  clusters <- as.character(obj[[cluster_col, drop = TRUE]])
  if (anyDuplicated(names(coarse_subsets_44542861)) ||
      !setequal(unique(clusters), names(coarse_subsets_44542861))) {
    stop("coarse_subsets_44542861 must assign each cluster exactly once.")
  }

  obj <- annotate_clusters_manual(obj, coarse_subsets_44542861,
    cluster_col = cluster_col, new_col_name = "subset")

  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, output)
  message("Saved: ", output, " (", file.size(output), " bytes)")

  # Reopen the actual output; metadata includes cell IDs and cluster assignments.
  reopened <- readRDS(output)
  stopifnot(
    identical(reopened[[]], obj[[]]),
    identical(SeuratObject::Idents(reopened), SeuratObject::Idents(obj)),
    identical(reopened@assays, obj@assays),
    identical(reopened@reductions, obj@reductions),
    identical(reopened@graphs, obj@graphs),
    identical(file.info(input)[, c("size", "mtime")], input_info)
  )
  message("Reopen checks passed: ", ncol(reopened),
    " cells; metadata, identities, assays/layers/counts, reductions and graphs preserved.")
  message("Subset annotation completed: ", output)
}
