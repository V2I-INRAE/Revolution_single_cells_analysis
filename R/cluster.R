# Leiden clustering on the chosen integrated reduction.

suppressPackageStartupMessages({
  library(Seurat)
})

cluster_object <- function(seurat_obj, method, clustering_config, dims) {
  reduction <- switch(method,
    harmony = "harmony",
    cca = "integrated_cca",
    stop("Unknown clustering reduction for method: ", method)
  )
  graphs <- paste0(method, c("_nn", "_snn"))
  set.seed(clustering_config$seed)
  seurat_obj <- FindNeighbors(
    seurat_obj,
    reduction = reduction,
    dims = dims,
    k.param = clustering_config$k_param,
    prune.SNN = clustering_config$prune_snn,
    nn.method = clustering_config$nn_method,
    annoy.metric = clustering_config$distance,
    n.trees = clustering_config$n_trees,
    graph.name = graphs
  )
  seurat_obj <- FindClusters(
    seurat_obj,
    graph.name = graphs[2],
    resolution = clustering_config$resolutions,
    algorithm = clustering_config$algorithm,
    leiden_method = clustering_config$leiden_method,
    leiden_objective_function = clustering_config$leiden_objective_function,
    random.seed = clustering_config$seed,
    n.iter = clustering_config$n_iter
  )
  partitions <- data.frame(
    resolution = clustering_config$resolutions,
    column = paste0(graphs[2], "_res.", clustering_config$resolutions)
  )
  seurat_obj@misc$clustering <- list(
    method = method,
    reduction = reduction,
    settings = clustering_config,
    dims = dims,
    partitions = partitions,
    graph = graphs[2],
    note = "All resolutions retained; active identity is the last, not a selected optimum"
  )
  seurat_obj
}

# Cluster counts per resolution, for logs and the .INFO stage record.
cluster_counts <- function(seurat_obj) {
  partitions <- seurat_obj@misc$clustering$partitions
  counts <- vapply(
    partitions$column,
    function(column) length(unique(seurat_obj[[]][[column]])),
    integer(1)
  )
  setNames(counts, partitions$resolution)
}

write_cluster_assignments <- function(seurat_obj, path) {
  partitions <- seurat_obj@misc$clustering$partitions
  assignments <- seurat_obj[[]][, partitions$column, drop = FALSE]
  utils::write.csv(
    data.frame(cell = rownames(assignments), assignments),
    path, row.names = FALSE
  )
  invisible(path)
}
