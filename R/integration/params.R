# Shared PCA/Harmony/CCA dimensions for integration, UMAPs and diagnostics.
integration_params <- list(
  dims = 1:20
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
