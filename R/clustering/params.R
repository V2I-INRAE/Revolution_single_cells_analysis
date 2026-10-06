clustering_params <- list(
  dims = 1:20,
  scvi_dim = 1:20,
  seed = 1234,
  k_param = 100,
  nn_method = "annoy",
  distance = "euclidean",
  n_trees = 50,
  prune_snn = 1 / 15,
  algorithm = 1, #c'est louvain
  resolutions = c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1),
  n_start = 10,
  n_iter = 10,
  n_diagnostic_cells = 50000, # c'est pour subset les cellules pour diagnostic de cluster
  silhouette_distance = "euclidean"
)
