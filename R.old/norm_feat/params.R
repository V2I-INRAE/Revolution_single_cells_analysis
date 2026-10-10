# Parameters controlling normalization, variable feature selection and PCA

norm_feat_params <- list(
  normalization = list(
    n_variable_features = 3000,
    vars_to_regress = c("percent.mt", "nFeature_RNA"),
    lognorm_method = "LogNormalize",
    lognorm_scale_factor = 10000,
    hvg_selection_method = "vst"
  ),
  pca = list(
    n_pcs = 50,
    seed = 1234
  )
)
