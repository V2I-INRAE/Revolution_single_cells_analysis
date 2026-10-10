suppressPackageStartupMessages(library(Seurat))

cluster_reduction <- function(obj, method, settings) {
  reduction <- switch(method,
    unintegrated = "pca",
    harmony = "harmony",
    cca = "integrated_cca"
  )
  graphs <- paste0(method, c("_nn", "_snn"))
  set.seed(settings$seed)
  obj <- FindNeighbors(obj,
    reduction = reduction, dims = settings$dims,
    k.param = settings$k_param, prune.SNN = settings$prune_snn,
    nn.method = settings$nn_method, annoy.metric = settings$distance,
    n.trees = settings$n_trees,
    graph.name = graphs
  )
  obj <- FindClusters(obj,
    graph.name = graphs[2],
    resolution = settings$resolutions, algorithm = settings$algorithm,
    leiden_method = settings$leiden_method,
    leiden_objective_function = settings$leiden_objective_function,
    random.seed = settings$seed, n.iter = settings$n_iter
  )
  partitions <- data.frame(
    resolution = settings$resolutions,
    column = paste0(graphs[2], "_res.", settings$resolutions)
  )
  obj@misc$clustering <- list(
    method = method, reduction = reduction,
    settings = settings, partitions = partitions,
    algorithm = settings$algorithm, graph = graphs[2], nn_method = settings$nn_method,
    leiden_method = settings$leiden_method,
    leiden_objective_function = settings$leiden_objective_function,
    distance = settings$distance, n_trees = settings$n_trees,
    n_iter = settings$n_iter,
    note = "All resolutions retained; active identity is the last, not a selected optimum"
  )
  obj
}

clustering_umap <- function(obj) {
  config <- obj@misc$clustering
  dims <- config$settings$dims
  name <- switch(config$method,
    unintegrated = "umap",
    harmony = "umap_harmony",
    cca = "umap_cca"
  )
  if (!name %in% Reductions(obj)) {
    obj <- RunUMAP(obj,
      reduction = config$reduction,
      dims = dims, reduction.name = name,
      n.neighbors = 30, min.dist = 0.3, metric = "cosine",
      seed.use = config$settings$seed
    )
  }
  commands <- Filter(
    function(command) identical(command@params$reduction.name, name),
    obj@commands
  )
  stopifnot(
    "Saved UMAP must document the selected reduction and dimensions" =
      length(commands) == 1L &&
        identical(commands[[1]]@params$reduction, config$reduction) &&
        identical(commands[[1]]@params$dims, dims)
  )
  coordinates <- Embeddings(obj, name)
  stopifnot(
    "UMAP must contain finite, aligned coordinates for every cell" =
      identical(rownames(coordinates), colnames(obj)) &&
        ncol(coordinates) == 2L && all(is.finite(coordinates))
  )
  obj@misc$clustering$umap <- name
  obj@misc$clustering$umap_settings <- commands[[1]]@params
  obj
}
