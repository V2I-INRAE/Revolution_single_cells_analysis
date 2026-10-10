# Parameters controlling QC metric identification, filtering and doublet calls

qc_params <- list(
  metrics = list(
    mitochondrial_genes = c(
      "ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB",
      "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"
    ),
    mitochondrial_pattern = "^MT-",
    ribosomal_pattern = "^RP[LS]"
  ),
  filtering = list(
    min_features = 300,
    max_features = 4000,
    min_log10_genes_per_umi = 0.8,
    max_mito_percent = 10,
    min_cells_per_feature = 10
  ),
  doublets = list(
    # BD Rhapsody Instrument User Guide, Doc ID 214062 Rev. 3.0, pp. 80–81:
    # captured cells on retrieved beads, not HT loaded cells per lane.
    # https://www.bdbiosciences.com/content/dam/bdb/marketing-documents/BD-Rhapsody-Single-Cell-Analysis-System-Instrument.pdf
    bd_multiplet_table = data.frame(
      cells = c(
        100, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000,
        9000, 10000, 11000, 12000, 13000, 14000, 15000, 16000, 17000,
        18000, 19000, 20000
      ),
      rate = c(
        0.0, 0.1, 0.2, 0.5, 0.7, 1.0, 1.2, 1.4, 1.7, 1.9,
        2.1, 2.4, 2.6, 2.8, 3.1, 3.3, 3.5, 3.8, 4.0,
        4.2, 4.5, 4.7
      )
    ),
    n_variable_features = 2000,
    pcs = 1:20,
    seed = 1234
  )
)
