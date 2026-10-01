# scVI model and training parameters; baseline PCA and UMAP keep their own dims.
scvi_params <- list(
  n_features = 3000L,
  n_latent = 20L,
  n_layers = 2L,
  gene_likelihood = "nb",
  max_epochs = NULL,
  seed = 1234L
)
