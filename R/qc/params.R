# Parameters controlling QC metric identification, filtering and doublet calls

qc_params <- list(
  metrics = list(
    mitochondrial_genes = c(
      "ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB",
      "ND1", "ND2", "ND3", "ND4", "ND5", "ND6"
    ),
    mitochondrial_pattern = "^MT-",
    ribosomal_pattern = "^RP[LS]",
    min_ribo_percent = 2.5
  ),
  filtering = list(
    min_features = 200,
    max_features_cap = 5000,
    feature_mad_multiplier = 4,
    min_log10_genes_per_umi = 0.8,
    max_mito_percent_cap = 20,
    mito_mad_multiplier = 4,
    min_cells_per_feature = 3
  ),
  doublets = list(
    bd_multiplet_table = data.frame(
      cells = c(
        100, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000,
        9000, 10000, 11000, 12000, 13000, 14000, 15000, 16000, 17000
      ),
      rate = c(
        0.0, 0.1, 0.2, 0.5, 0.7, 1.0, 1.2, 1.4, 1.7, 1.9,
        2.1, 2.4, 2.6, 2.8, 3.1, 3.3, 3.5, 3.8, 4.0
      )
    ),
    n_variable_features = 2000,
    pcs = 1:20,
    cluster_resolution = 0.6,
    cluster_algorithm = 1,
    pN = 0.25,
    expected_rate = NULL,
    seed = 1234,
    sct = FALSE,
    ground_truth = FALSE,
    min_expected_doublets = 1
  )
)
