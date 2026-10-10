# Unified configuration for the single-entry pipeline (QC -> normalization ->
# integration -> clustering). This file IS the run definition: edit it, launch
# scripts/sbatch_pipeline.sh, and the resolved values are recorded in .INFO.

config <- list(
  run = list(
    # Raw count source used for every sample: rhapsody | soupx | cellbender.
    input_source = "soupx",
    samples = c(
      "REVO26-C", "REVO26-N4", "REVO26-N10", "REVO26-P4", "REVO26-P10",
      "REVO27-C", "REVO27-N4", "REVO27-N10", "REVO27-P4", "REVO27-P10",
      "REVO29-C", "REVO29-N4", "REVO29-P4", "REVO30-C", "REVO30-N4",
      "REVO30-N10", "REVO30-P10", "REVO30-P4", "REVO31-C", "REVO31-N4",
      "REVO31-N8", "REVO31-P10", "REVO31-P4"
    ),
    seed = 1234
  ),

  qc = list(
    metrics = list(
      # Pig mitochondrial gene names lack the MT- prefix and are renamed on load.
      mitochondrial_genes = c(
        "ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB",
        "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"
      ),
      mitochondrial_pattern = "^MT-",
      ribosomal_pattern = "^RP[LS]"
    ),
    # Adaptive QC thresholds: per-sample median + k * MAD (R's mad, constant
    # 1.4826), capped by fixed ceilings; cells AT a limit pass, strictly above
    # is excluded. Ribosomal content keeps a fixed floor: cells strictly above
    # it are retained. Complexity keeps its fixed strict minimum.
    filtering = list(
      min_features = 200,
      max_features_cap = 5000,
      feature_mads = 3,
      max_mito_cap = 20,
      mito_mads = 3,
      min_ribo_percent = 2,
      min_log10_genes_per_umi = 0.8,
      min_cells_per_feature = 10
    ),
    doublets = list(
      # BD Rhapsody Instrument User Guide, Doc ID 214062 Rev. 3.0, pp. 80-81:
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
      # Temporary per-sample PCA/UMAP used for QC diagnostics only.
      n_variable_features = 2000,
      pcs = 1:20,
      seed = 1234,
      # Scrublet calls on all input cells, before any filtering.
      # seed 0 is recorded for provenance only: the pinned scrubletR wrapper
      # fixes the Python random_state to 0 and exposes no seed argument.
      scrublet = list(
        min_counts = 3L,
        n_prin_comps = 30L,
        threshold = 0.15,
        seed = 0L
      )
    )
  ),

  norm = list(
    n_variable_features = 3000,
    vars_to_regress = c("percent.mt", "nFeature_RNA"),
    lognorm_method = "LogNormalize",
    lognorm_scale_factor = 10000,
    hvg_selection_method = "vst",
    n_pcs = 50
  ),

  # Integration route: harmony (default) | cca.
  integration = list(
    method = "harmony",
    harmony = list(max_iter = 10L)
  ),

  clustering = list(
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
)

# Dimensions shared by integration, all UMAPs and clustering: one knob.
config$dims <- 1:20

# Startup validation: the pipeline is fixed, so reject impossible
# configurations before any expensive work.
validate_config <- function(config) {
  stopifnot(
    "run$input_source must be rhapsody, soupx or cellbender" =
      config$run$input_source %in% c("rhapsody", "soupx", "cellbender"),
    "run$samples must be nonempty and unique" =
      length(config$run$samples) > 0L && !anyDuplicated(config$run$samples),
    "integration$method must be harmony or cca" =
      config$integration$method %in% c("harmony", "cca"),
    "dims must be the contiguous range 1..n" =
      identical(as.integer(config$dims), seq_len(max(config$dims))),
    "dims cannot exceed the number of computed PCs" =
      max(config$dims) <= config$norm$n_pcs,
    "clustering resolutions must be positive and unique" =
      length(config$clustering$resolutions) > 0L &&
        all(config$clustering$resolutions > 0) &&
        !anyDuplicated(config$clustering$resolutions),
    "vars_to_regress must exist in QC metadata" =
      all(config$norm$vars_to_regress %in%
        c("percent.mt", "percent.ribo", "nFeature_RNA", "nCount_RNA",
          "log10GenesPerUMI"))
  )
  invisible(TRUE)
}
