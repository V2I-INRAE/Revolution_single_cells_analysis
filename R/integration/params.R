# Shared PCA/Harmony/CCA dimensions for integration, UMAPs and diagnostics.
integration_params <- list(
  dims = 1:20
)

harmony_params <- list(
  max_iter = 30L
)

# UMAP-only exploration on an existing Harmony checkpoint.
umap_exploration_params <- list(
  reduction = "harmony",
  dims = 1:20,
  umap.method = "uwot",
  metric = "cosine",
  n.neighbors = 30L,
  n.epochs = 200L,
  min.dist = 0.5,
  spread = 1,
  repulsion.strength = 1,
  seed.use = 1234L
)

# scVI trains its own latent representation from RNA counts, not PCA.
scvi_params <- list(
  n_features = 3000L,
  n_latent = 20L,
  n_layers = 2L,
  gene_likelihood = "nb",
  max_epochs = NULL,
  seed = 1234L
)
