# Unified configuration for the single-entry pipeline (QC -> normalization ->
# integration -> clustering). This file IS the run definition: edit it, launch
# scripts/sbatch_pipeline.sh, and the resolved values are recorded in .INFO.

# Feature panels for the marker stage live in their own file.
source("R/clusters/markers.R")

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
      ribosomal_pattern = "^RP[LS]",
      # Hemoglobin genes by explicit list: a ^HB pattern would also match
      # HBEGF, HBP1 and HBS1L, which are not globins.
      hemoglobin_genes = c("HBB", "HBE1", "HBM", "HBQ1", "HBZ")
    ),
    # Adaptive QC thresholds: per-sample median + k * MAD (R's mad, constant
    # 1.4826), capped by fixed ceilings; cells AT a limit pass, strictly above
    # is excluded. Ribosomal content keeps a fixed floor: cells strictly above
    # it are retained. Complexity keeps its fixed strict minimum. Hemoglobin
    # keeps a fixed ceiling: cells AT or above it are excluded (deliberately
    # inclusive, unlike every other ceiling).
    filtering = list(
      min_features = 200,
      max_features_cap = 5000,
      feature_mads = 4,
      max_mito_cap = 20,
      mito_mads = 4,
      max_hb_percent = 0.5,
      min_ribo_percent = 0,
      min_log10_genes_per_umi = 0.8,
      min_cells_per_feature = 10
    ),
    # Feature names dropped from the matrix after cell filtering, exactly as
    # they appear in the rownames at that point (the mitochondrial genes are
    # MT- prefixed by the on-load rename). From the top-20 genes by per-cell
    # percentage of counts in the 5-sample QC pilot: the mitochondrial and
    # ribosomal entries of that top 20, plus two Ensembl IDs kept as bare
    # feature names. Matched exactly by membership, see remove_listed_genes().
    remove_genes = c(
      "MT-ATP6", "MT-COX1", "MT-COX2", "MT-COX3", "MT-CYTB",
      "MT-ND1", "MT-ND2", "MT-ND4",
      "RPL37", "RPS12", "RPS24", "RPS8",
      "ENSSSCG00000052671", "ENSSSCG00000048856"
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
      # Doublet detection with scDblFinder on each sample's counts, all cells,
      # before any filtering. Parameter values below are PROPOSED and await
      # owner approval; alternatives are listed in the run record:
      # - dbr: "bd_table" derives the per-sample rate from the BD table above
      #   (the derivation used in the Scrublet era); numeric = fixed fraction;
      #   NULL = the package's 10X rate prior.
      # - dbr_sd = 1 removes the rate-deviation penalty in thresholding.
      # - clusters = FALSE classifies without internal clustering (this stage
      #   has no cluster structure yet); TRUE would cluster internally.
      # Scores are scDblFinder's, not Scrublet's: calls come from
      # scDblFinder.class, never the old Scrublet score threshold.
      scdblfinder = list(
        dbr = "bd_table",
        dbr_sd = 1,
        clusters = FALSE,
        seed = 1234L
      )
    )
  ),

  norm = list(
    n_variable_features = 3000,
    # S.Score and G2M.Score are added by the cell-cycle step before scaling.
    vars_to_regress = c("percent.mt", "nFeature_RNA", "S.Score", "G2M.Score"),
    lognorm_method = "LogNormalize",
    lognorm_scale_factor = 10000,
    hvg_selection_method = "vst",
    n_pcs = 50
  ),

  # Cell-cycle scoring (normalization stage): cells above their sample's
  # cutoffs are flagged (cc_tail) and removed before HVG selection, so
  # proliferating tail cells cannot drive features, PCs or clusters. Cutoffs
  # from the qc-explore-all23 exploration, approved on the 2026-10-11 transfer
  # test at production QC. The rename map puts assayable pig Ensembl IDs
  # under their human cc.genes symbols; MLF1IP, FAM64A and HN1 are absent
  # from the BD panel and are logged no-ops.
  cell_cycle = list(
    rename_map = c(
      ENSSSCG00000015785 = "MLF1IP",
      ENSSSCG00000031741 = "RRM2",
      ENSSSCG00000026302 = "MKI67",
      ENSSSCG00000008677 = "TACC3",
      ENSSSCG00000017893 = "FAM64A",
      ENSSSCG00000017213 = "HN1",
      ENSSSCG00000007493 = "AURKA"
    ),
    ctrl = 43,
    seed = 1234,
    tail_cutoffs = data.frame(
      sample = c(
        "REVO26-C", "REVO26-N4", "REVO26-N10", "REVO26-P4", "REVO26-P10",
        "REVO27-C", "REVO27-N4", "REVO27-N10", "REVO27-P4", "REVO27-P10",
        "REVO29-C", "REVO29-N4", "REVO29-P4", "REVO30-C", "REVO30-N4",
        "REVO30-N10", "REVO30-P10", "REVO30-P4", "REVO31-C", "REVO31-N4",
        "REVO31-N8", "REVO31-P10", "REVO31-P4"
      ),
      cut_S = c(
        0.258, 0.352, 0.324, 0.257, 0.296,
        0.380, 0.380, 0.380, 0.380, 0.380,
        0.371, 0.371, 0.371, 0.320, 0.320,
        0.320, 0.320, 0.320, 0.380, 0.380,
        0.380, 0.380, 0.380
      ),
      cut_G2M = c(
        0.248, 0.248, 0.248, 0.248, 0.248,
        0.295, 0.295, 0.295, 0.295, 0.295,
        0.236, 0.236, 0.236, 0.250, 0.250,
        0.250, 0.250, 0.250, 0.250, 0.250,
        0.250, 0.250, 0.250
      )
    )
  ),

  # Integration route: harmony (default) | cca.
  integration = list(
    method = "harmony",
    harmony = list(max_iter = 10L),
    # Shared by the unintegrated PCA and integrated-method UMAPs.
    umap = list(min.dist = 0.3, n.neighbors = 30L)
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
    resolutions = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7),
    n_iter = 10
  ),

  # Marker discovery and plotting share R/clusters/find_markers.R. The input
  # object and resolution arrive as launcher arguments; results go to
  # results/<job>/clusters_plots/markers/<resolution>/.
  markers = list(
    cluster_prefix = "harmony_snn_res.",
    reduction = "umap_harmony",
    test_use = "wilcox",
    only_pos = TRUE,
    min_pct = 0.25,
    logfc_threshold = 0.25,
    adjusted_p_cutoff = 0.05,
    # Panels and ID-only feature translations defined in R/clusters/markers.R
    # (sourced above); visualization panels, not cluster annotations.
    feature_panels = marker_feature_panels,
    # Applied to the loaded object BEFORE FindAllMarkers so tables and
    # figures share gene names; collision-checked via rename_cc_features().
    feature_renames = marker_feature_renames
  )
)

# Dimensions shared by integration, all UMAPs and clustering: one knob.
config$dims <- 1:30

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
    "UMAP min.dist must be finite and nonnegative" =
      length(config$integration$umap$min.dist) == 1L &&
        is.numeric(config$integration$umap$min.dist) &&
        is.finite(config$integration$umap$min.dist) &&
        config$integration$umap$min.dist >= 0,
    "UMAP n.neighbors must be an integer of at least 2" =
      length(config$integration$umap$n.neighbors) == 1L &&
        is.numeric(config$integration$umap$n.neighbors) &&
        is.finite(config$integration$umap$n.neighbors) &&
        config$integration$umap$n.neighbors >= 2 &&
        config$integration$umap$n.neighbors %% 1 == 0,
    "dims must be the contiguous range 1..n" =
      identical(as.integer(config$dims), seq_len(max(config$dims))),
    "dims cannot exceed the number of computed PCs" =
      max(config$dims) <= config$norm$n_pcs,
    "clustering resolutions must be positive and unique" =
      length(config$clustering$resolutions) > 0L &&
        all(config$clustering$resolutions > 0) &&
        !anyDuplicated(config$clustering$resolutions),
    "vars_to_regress must exist in QC metadata or be cell-cycle scores" =
      all(config$norm$vars_to_regress %in%
        c("percent.mt", "percent.ribo", "percent.hb", "nFeature_RNA",
          "nCount_RNA", "log10GenesPerUMI", "S.Score", "G2M.Score")),
    "cell_cycle$tail_cutoffs must match run$samples exactly" =
      !anyNA(config$cell_cycle$tail_cutoffs$sample) &&
      identical(sort(as.character(config$cell_cycle$tail_cutoffs$sample)),
        sort(config$run$samples)),
    "cell_cycle cutoffs must be finite and positive" =
      all(is.finite(c(config$cell_cycle$tail_cutoffs$cut_S,
        config$cell_cycle$tail_cutoffs$cut_G2M))) &&
      all(c(config$cell_cycle$tail_cutoffs$cut_S,
        config$cell_cycle$tail_cutoffs$cut_G2M) > 0),
    "cell_cycle$rename_map must have unique sources and targets" =
      !anyDuplicated(c(names(config$cell_cycle$rename_map),
        unname(config$cell_cycle$rename_map))),
    "qc$doublets$scdblfinder$dbr must be \"bd_table\", NULL, or one fraction in [0, 1]" =
      identical(config$qc$doublets$scdblfinder$dbr, "bd_table") ||
        is.null(config$qc$doublets$scdblfinder$dbr) ||
        (is.numeric(config$qc$doublets$scdblfinder$dbr) &&
          length(config$qc$doublets$scdblfinder$dbr) == 1L &&
          is.finite(config$qc$doublets$scdblfinder$dbr) &&
          config$qc$doublets$scdblfinder$dbr >= 0 &&
          config$qc$doublets$scdblfinder$dbr <= 1),
    "qc$doublets$scdblfinder$dbr_sd must be NULL or one absolute rate tolerance in [0, 1]" =
      is.null(config$qc$doublets$scdblfinder$dbr_sd) ||
        (is.numeric(config$qc$doublets$scdblfinder$dbr_sd) &&
          length(config$qc$doublets$scdblfinder$dbr_sd) == 1L &&
          is.finite(config$qc$doublets$scdblfinder$dbr_sd) &&
          config$qc$doublets$scdblfinder$dbr_sd >= 0 &&
          config$qc$doublets$scdblfinder$dbr_sd <= 1),
    "qc$doublets$scdblfinder$clusters must be TRUE or FALSE" =
      isTRUE(config$qc$doublets$scdblfinder$clusters) ||
        isFALSE(config$qc$doublets$scdblfinder$clusters),
    "qc$doublets$scdblfinder$seed must be one nonnegative integer" =
      is.numeric(config$qc$doublets$scdblfinder$seed) &&
        length(config$qc$doublets$scdblfinder$seed) == 1L &&
        is.finite(config$qc$doublets$scdblfinder$seed) &&
        config$qc$doublets$scdblfinder$seed >= 0 &&
        config$qc$doublets$scdblfinder$seed %% 1 == 0
  )
  invisible(TRUE)
}
