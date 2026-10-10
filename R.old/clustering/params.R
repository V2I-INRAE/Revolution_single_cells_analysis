clustering_params <- list(
  dims = 1:20,
  seed = 1234,
  k_param = 100,
  nn_method = "annoy",
  distance = "euclidean",
  n_trees = 50,
  prune_snn = 1 / 15,
  algorithm = 4, # Leiden
  leiden_method = "leidenbase",
  leiden_objective_function = "modularity",
  resolutions = c(0.1, 0.2, 0.3, 0.4),
  n_iter = 10
)
