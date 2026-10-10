# Shared UMAP-only exploration; callers own input selection, plots and run records.
umap_exploration_params <- list(
  dims = 1:20,
  umap.method = "uwot",
  metric = "cosine",
  n.neighbors = 30L,
  n.epochs = 200L,
  min.dist = 0.3,
  spread = 1,
  repulsion.strength = 1,
  seed.use = 1234L
)

umap_exploration_tests <- local({
  tests <- list(
    min.dist_0.5 = list(min.dist = 0.5),
    min.dist_0.7 = list(min.dist = 0.7),
    min.dist_1 = list(min.dist = 1),
    n.neighbors_50 = list(n.neighbors = 50L),
    n.neighbors_75 = list(n.neighbors = 75L),
    repulsion.strength_1.5 = list(repulsion.strength = 1.5),
    repulsion.strength_2 = list(repulsion.strength = 2)
  )
  combinations <- expand.grid(
    min.dist = c(0.7, 1), n.neighbors = c(50L, 75L),
    repulsion.strength = c(1.5, 2)
  )
  for (i in seq_len(nrow(combinations))) {
    settings <- as.list(combinations[i, ])
    name <- paste(paste(names(settings), unlist(settings), sep = "_"), collapse = "_")
    tests[[name]] <- settings
  }
  tests
})

run_umap_exploration <- function(obj, params, reduction_name) {
  stopifnot("Input reduction must exist" = params$reduction %in% Seurat::Reductions(obj))
  coordinates <- Seurat::Embeddings(obj, params$reduction)
  stopifnot("Input coordinates must be finite and aligned with cells" =
    identical(rownames(coordinates), colnames(obj)) &&
      max(params$dims) <= ncol(coordinates) &&
      all(is.finite(coordinates[, params$dims, drop = FALSE])))
  if (reduction_name %in% Seurat::Reductions(obj)) {
    stop("Reduction already exists; refusing to overwrite: ", reduction_name)
  }
  obj <- Seurat::RunUMAP(
    object = obj, reduction.name = reduction_name,
    reduction.key = paste0(gsub("[^[:alnum:]]", "", reduction_name), "_"),
    reduction = params$reduction, dims = params$dims,
    umap.method = params$umap.method, metric = params$metric,
    n.neighbors = params$n.neighbors, n.epochs = params$n.epochs,
    min.dist = params$min.dist, spread = params$spread,
    repulsion.strength = params$repulsion.strength, seed.use = params$seed.use,
    verbose = TRUE
  )
  coordinates <- Seurat::Embeddings(obj, reduction_name)
  stopifnot("Tested UMAP must be finite, two-dimensional and aligned with cells" =
    identical(rownames(coordinates), colnames(obj)) &&
      ncol(coordinates) == 2L && all(is.finite(coordinates)))
  obj
}
