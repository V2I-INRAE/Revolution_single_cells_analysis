clustering_params <- list(
  dims = 1:30,
  scvi_dim = 1:20,
  seed = 1234,
  k_param = 100,
  nn_method = "annoy",
  distance = "euclidean",
  n_trees = 50,
  prune_snn = 1 / 15,
  algorithm = 1, #c'est louvain
  resolutions = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3),
  n_start = 10,
  n_iter = 10,
  n_diagnostic_cells = 50000, # c'est pour subset les cellules pour diagnostic de cluster
  silhouette_distance = "euclidean"
)

marker_params <- list(
  resolution = 0.1,
  panels = list(
    Stromal = list(
      `Adventitial fibroblast` = c("PI16", "DPT", "COL14A1"), # X
      `Alveolar/interstitial fibroblast` = c("WNT2", "NPNT", "COL1A1"), # X
      `Activated fibroblast / myofibroblast` = c("COL3A1", "ACTA2", "CTHRC1"), # P; injury-associated
      Pericyte = c("PDGFRB", "RGS5", "CSPG4"), # X
      `Pleural mesothelial cell` = c("MSLN", "WT1", "UPK3B") # X; only if pleura sampled
    ),
    Mesenchymal = list(
      `Airway or vascular smooth-muscle cell` = c("ACTA2", "TAGLN", "MYH11") # X
    ),
    Epithelial = list(
      `Alveolar type 1` = c("AGER", "AQP5", "PDPN"), # P; porcine AQP5 support
      `Alveolar type 2` = c("SFTPC", "SFTPB", "ABCA3"), # P; porcine SFTPC support
      `Basal airway cell` = c("KRT5", "TP63", "KRT14"), # X
      `Club/secretory airway cell` = c("SCGB1A1", "SCGB3A2", "BPIFB1"), # X
      `Ciliated airway cell` = c("FOXJ1", "PIFO", "TPPP3"), # X
      `Goblet cell` = c("MUC5AC", "SPDEF", "AGR2"), # X
      `Pulmonary neuroendocrine cell` = c("CHGA", "ASCL1", "SYP") # X; rare
    ),
    Endothelial = list(
      `Capillary/aerocyte-like cell` = c("CA4", "EMCN", "RGCC"), # X
      `Arterial cell` = c("GJA5", "EFNB2", "SOX17"), # X
      `Venous cell` = c("ACKR1", "NR2F2", "VWF"), # X
      `Lymphatic cell` = c("PROX1", "LYVE1", "FLT4") # X
    ),
    Immune = list(
      `Alveolar macrophage` = c("FABP4", "MARCO", "PPARG"), # P for population; panel X
      `Interstitial macrophage` = c("C1QA", "CD163", "MRC1"), # P for population; panel X
      `Recruited/inflammatory monocyte` = c("FCN1", "S100A8", "LYZ"), # P during infection; panel X
      Neutrophil = c("S100A8", "S100A9", "CSF3R"), # P during infection; panel X
      `Conventional dendritic cell, cDC2-like` = c("CD1C", "FCER1A", "CLEC10A"), # X
      `Plasmacytoid dendritic cell` = c("TCF4", "IL3RA", "IRF7"), # P during infection; panel X
      `CD4 T cell` = c("CD3D", "CD4", "IL7R"), # P for T cells; panel X
      `Cytotoxic CD8 T cell` = c("CD8A", "CD8B", "GZMB"), # P; GZMB depends on activation
      `Natural killer cell` = c("NKG7", "KLRD1", "GNLY"), # X
      `B cell` = c("CD79A", "MS4A1", "CD79B"), # P for population; panel X
      `Plasma cell` = c("JCHAIN", "MZB1", "XBP1"), # P after respiratory infection; panel X
      `Mast cell` = c("KIT", "TPSAB1", "CPA3") # X; check pig gene annotation
    )
  )
)
