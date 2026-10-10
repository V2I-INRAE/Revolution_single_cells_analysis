# Sample-batch integration: harmony (default) or CCA, behind one method switch.

suppressPackageStartupMessages({
  library(Seurat)
})

# Minimal handoff contract between normalization and integration: PCA aligned
# and finite on the shared dims, complete sample labels, HVGs consistent with
# the PCA loadings, per-sample count layers intact.
validate_integration_input <- function(seurat_obj, dims) {
  metadata <- seurat_obj[[]]
  pca <- Embeddings(seurat_obj, "pca")
  features <- VariableFeatures(seurat_obj)
  stopifnot(
    identical(as.integer(dims), seq_len(max(dims))),
    max(dims) <= ncol(pca),
    all(is.finite(pca[, dims, drop = FALSE])),
    identical(rownames(pca), colnames(seurat_obj)),
    "sample" %in% names(metadata),
    !anyNA(metadata$sample),
    length(unique(metadata$sample)) > 1L,
    length(features) > 0L,
    all(features %in% rownames(Loadings(seurat_obj, "pca"))),
    length(Layers(seurat_obj[["RNA"]], search = "^counts($|\\.)")) > 0L
  )
  message(
    "Integration input: ", ncol(seurat_obj), " cells, ",
    length(unique(metadata$sample)), " samples, assay ", DefaultAssay(seurat_obj)
  )
  invisible(TRUE)
}

harmony_convergence <- function(messages, max_iter) {
  messages <- trimws(messages)
  converged <- grep(
    "^Harmony converged after [0-9]+ iterations$",
    messages, value = TRUE
  )
  if (length(converged)) {
    return(tail(converged, 1L))
  }
  if (any(grepl(paste0("^Harmony ", max_iter, "/", max_iter), messages))) {
    return(paste0(
      "Reached the ", max_iter, "-iteration limit without a convergence message"
    ))
  }
  "Convergence not captured from Harmony output"
}

integrate_object <- function(seurat_obj, integration_config, dims, seed) {
  validate_integration_input(seurat_obj, dims)
  method <- integration_config$method

  if (method == "harmony") {
    set.seed(seed)
    # Harmony reports convergence via message(); record the lines while
    # letting them print to the normal execution log.
    harmony_log <- character()
    seurat_obj <- withCallingHandlers(
      harmony::RunHarmony(
        object = seurat_obj,
        group.by.vars = "sample",
        reduction.use = "pca",
        dims.use = dims,
        reduction.save = "harmony",
        project.dim = FALSE,
        max_iter = integration_config$harmony$max_iter,
        verbose = TRUE
      ),
      message = function(m) {
        harmony_log <<- c(harmony_log, conditionMessage(m))
      }
    )
    convergence <- harmony_convergence(
      harmony_log, integration_config$harmony$max_iter
    )
  } else if (method == "cca") {
    set.seed(seed)
    seurat_obj <- IntegrateLayers(
      object = seurat_obj,
      method = CCAIntegration,
      assay = DefaultAssay(seurat_obj),
      normalization.method = "LogNormalize",
      features = VariableFeatures(seurat_obj),
      orig.reduction = "pca",
      new.reduction = "integrated_cca",
      dims = dims,
      dims.to.integrate = dims,
      verbose = TRUE
    )
    convergence <- NULL
  } else {
    stop("Unknown integration method: ", method)
  }

  seurat_obj@misc$integration <- list(
    method = method,
    batch_col = "sample",
    dims = dims,
    seed = seed,
    harmony_max_iter = if (method == "harmony") integration_config$harmony$max_iter else NULL,
    harmony_convergence = convergence,
    purpose = "sample correction ahead of clustering"
  )
  seurat_obj
}

# UMAP name/reduction pairs for the chosen method: the unintegrated PCA UMAP
# always sits beside the method UMAP.
integration_reductions <- function(method) {
  c(
    umap = "pca",
    switch(method,
      harmony = c(umap_harmony = "harmony"),
      cca = c(umap_cca = "integrated_cca"),
      stop("Unknown integration method: ", method)
    )
  )
}
