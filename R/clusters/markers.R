# Feature panels for the marker stage (R/clusters/find_markers.R).
# Sourced by R/config.R; output is <Panel>/<gene>.png and <Panel>/<gene>.svg.
#
# Selected 2026-10-13 from job-44493014, resolution 0.7:
# results/job-44493014/clusters_plots/markers/0.7/markers_top25_plus_lineage_by_pct.csv.
# Each panel's first five genes follow that CSV's pct.1 order; the next five
# have the highest avg_log2FC among remaining genes (ties: p_val_adj, gene).
# This is the top-25-plus-lineage candidate pool, not all tested genes.
# Owner-approved exceptions: LYVE1 in Lymphatic_EC; FCN1, GPNMB and
# SIGLEC1 in Macrophages_AM_like. Both panels retain ten genes.
# Owner-approved screenshot additions (2026-10-14, annotation-thread request,
# remade in this thread; see .cluster_markers "Screenshot additions"):
# Monocytes_inflammatory +SOD2/CSF1R/C5AR1/TGFB1/CX3CR1;
# Macrophages_AM_like +PON2/TLR4/CD68; Macrophages_C1Q_high +MAF/CCL2/C5AR1/
# TLR4/CSF1R/LAP3/MARCO/IDO1; cDC1 +CD86/CCR7/CD1E; pDC +CLEC12A/CCR2/RUNX2.
#
# One ten-gene panel per cluster, including unresolved mixed cluster 22.
# Cycling_state removed by owner request. These ranked lists replace manual
# panel pruning. Panels are visualization aids, not cluster annotations or
# evidence of lineage specificity. No expression threshold is applied.
# Coverage and naming notes live in .cluster_markers beside the input CSV.
# ID-only feature translations below apply to the in-memory object before
# FindAllMarkers, so tables and plots share names; the saved object is unchanged.

marker_feature_panels <- list(
  # ---- Myeloid ------------------------------------------------------------
  Monocytes_inflammatory = c(
    "MAFB", "TBXAS1", "F13A1", "FCN1", "ADGRE4 (review)",
    "SAA3", "IL1B", "PLBD1", "AREG", "EREG",
    # Owner-approved screenshot Mo-C/MP additions 2026-10-13:
    "SOD2", # cl 9 log2FC 2.24/pct.1 0.787
    "CSF1R", # cl 9 log2FC 4.47/pct.1 0.640; also cl 13
    "C5AR1", # cl 9 pct.1 0.757; also cl 13 (shared signal, FCN1 precedent)
    "TGFB1", # screenshot MP-family; significant in cl 9 here (pct.1 0.690)
    "CX3CR1" # screenshot Mo-C; absent from object (warn-skip)
  ),
  Macrophages_AM_like = c(
    "LGALS3", "APOE", "S100A12", "FCN2", "PTGR1",
    "CIDEC", "FCGR1A", "FCN1", "GPNMB", "SIGLEC1",
    # Owner-approved screenshot-MP additions 2026-10-13:
    "PON2", # cl 3 log2FC 2.21/pct.1 0.513
    "TLR4", # cl 3 pct.1 0.917 (screenshot TLR4hi); also cl 13
    "CD68" # cl 3 log2FC 4.11/pct.1 0.870
  ),
  Macrophages_C1Q_high = c(
    "C1QA", "C1QC", "C1QB", "FCGR2B", "CCL14",
    "LOC100524972", "CCL8", "LOC100511343", "LOC100624077", "LOC100525275 (review)",
    # Owner-approved screenshot Mo-C/MP additions 2026-10-13:
    "MAF", # cl 13 log2FC 2.48/pct.1 0.795
    "CCL2", # cl 13 log2FC 4.05/pct.1 0.693
    "C5AR1", # cl 13 pct.1 0.885; also cl 9
    "TLR4", # cl 13 pct.1 0.717; also cl 3
    "CSF1R", # cl 13 log2FC 3.76/pct.1 0.691
    "LAP3", # screenshot MP; cl 13 pct.1 0.753 (broad signal)
    "MARCO", # owner request; absent from object (warn-skip)
    "IDO1" # screenshot MP; absent from object (warn-skip)
  ),
  cDC1 = c(
    "SLA-DQB1", "SLA-DRB1 (review)", "HLA-DRA", "CD83", "FLT3",
    "XCR1", "TIFAB", "SLC6A12", "DNASE1L3", "RAMP1",
    # Owner-approved screenshot Mo-DC additions 2026-10-13:
    "CD86", # cl 18 log2FC 2.88/pct.1 0.805
    "CCR7", # cl 18 log2FC 2.59/pct.1 0.408
    "CD1E" # screenshot cDC2 comparison; absent from object (warn-skip)
  ),
  pDC = c(
    "CST3", "IRF8", "BLNK", "LOC100620627", "TCF4",
    "VPREB1", "PIGE-108A11.3", "TDRD1", "ENSSSCG00000061277", "S100Z",
    # Owner-approved screenshot additions 2026-10-13:
    "CLEC12A", # cl 14 log2FC 5.39/pct.1 0.686
    "CCR2", # cl 14 log2FC 1.90/pct.1 0.333 (screenshot Mo-DC)
    "RUNX2" # screenshot pDC; cl 14 log2FC 1.49/pct.1 0.815
  ),
  Neutrophils = c(
    "LOC100124654 (review)", "CXCL8", "CSF3R", "TGM3", "SDS",
    "LOC100515345", "STEAP4", "GPR84", "KCNJ15", "G0S2"
  ),
  Eosinophils = c(
    "RBBP7", "EMB", "MCTP1", "ALOX5", "RGS18",
    "EPX", "LOC100158011", "LOC100522787", "CCR3", "ENSSSCG00000034853"
  ),

  # ---- Lymphoid -----------------------------------------------------------
  T_CD4_naive_CM = c(
    "ABLIM1", "CD5", "CD28", "IL7R", "TCF7",
    "CERS6", "ICOS", "LEF1", "CD40LG", "TSPAN5"
  ),
  T_CD8_cytotoxic = c(
    "CD2", "CCL5", "CD3D", "LOC102158035", "CD3E",
    "CD6", "CD8A", "CXCR6", "CD5", "GZMB"
  ),
  T_gd_like = c(
    "CD247", "GATA3", "SYTL3", "RHEX", "TCF7",
    "CD163L1", "SOX13", "BLK", "ENSSSCG00000056250", "SH3BP4"
  ),
  NK_like = c(
    "IL2RB", "KLRK1", "GZMA", "UBD", "TXK",
    "NCR1", "GPR37", "KLRB1", "SEMA6D", "LOC100523789"
  ),
  B_cells = c(
    "SMC6", "CD79B", "MS4A1", "RALGPS2", "BANK1",
    "CD19", "CLEC17A", "ITGAD", "PAX5", "ELL3"
  ),
  ASC = c(
    "JCHAIN", "ENSSSCG00000036224", "TXNDC5", "PRDM1", "SEC11C",
    "ENSSSCG00000053053", "ENSSSCG00000058302", "ENSSSCG00000063511",
    "ENSSSCG00000048391", "ENSSSCG00000031037"
  ),

  # ---- Endothelium and stroma ---------------------------------------------
  Endothelial_arterial = c(
    "PECAM1", "CLEC14A", "CLDN5", "SEMA3G", "EFNB2",
    "GJA5", "SLC52A3", "SSTR1", "CX3CL1", "HTR2B"
  ),
  Endothelial_venous = c(
    "PECAM1", "ADGRF5", "MMRN2", "FMO1", "EMCN",
    "ACKR1", "SLC16A9", "FAM167B", "CCDC85A", "CLEC1A"
  ),
  Lymphatic_EC = c(
    "MMRN1", "APOLD1", "PDPN", "CCL21", "MYCT1",
    "HOXD9", "HOXD8", "CHP2", "LYVE1", "FLT4" # LYVE1 replaces SCG3 2026-10-13: canonical lymphatic marker, c12 pct.1 0.772/pct.2 0.02
  ),
  Fibroblasts = c(
    "COL1A1", "COL3A1", "DCN", "COL1A2", "FBLN1",
    "ENSSSCG00000031292", "PI16", "DPT", "LUM", "CHRDL1"
  ),
  Smooth_muscle = c(
    "SPARCL1", "TAGLN", "NDNF", "MYH11", "WIF1",
    "GREM1", "CNN1", "ACTG2", "ENSSSCG00000058753", "ACTA2"
  ),

  # ---- Epithelium ---------------------------------------------------------
  Secretory_epithelium = c(
    "CBR2", "TFF3", "WFDC2", "SULT1E1", "PIGR",
    "SCGB1A1", "SCGB3A1", "BPIFB1", "HSD17B2", "ENSSSCG00000013714"
  ),
  Ciliated = c(
    "DNAH11", "CIMIP1", "RSPH4A", "C1orf87", "RSPH1",
    "SNTN", "UBXN10", "CIBAR2", "CFAP144", "CIMAP1B"
  ),
  Tuft = c(
    "SH2D7", "POU2F3", "RGS21", "RGS13", "AVIL",
    "HMX3", "CYP4A90", "HMX2", "GNG13", "IHH"
  ),
  AT2 = c(
    "SFTPA1", "SFTPC", "SFTPD", "LAMP3", "NAPSA",
    "GJB1", "IYD", "C5", "CHIA", "SLC34A2"
  ),
  AT1 = c(
    "AGER", "HOPX", "LMO7", "PDPN", "ADIRF",
    "LOC100738720", "COL4A3", "COL4A4", "VSIG2", "FLRT3"
  ),

  # ---- Unresolved mixed population ---------------------------------------
  Mixed_IFN_myeloid_T_unresolved = c(
    "LOC100124654 (review)", "SDS", "MARCKS", "HCAR2", "LOC100515345",
    "VNN1", "IFIT2", "IFIT1", "STEAP4", "IFIT3"
  )
)

# Preserve the earlier 19 translations and add selected CSV names from
# results/misc/Ensembl_ID/gene-mapping/rename_map.csv. Review suffixes follow
# symbol_review_required; ambiguous or unnamed features retain their IDs.
# The existing rename_cc_features() rejects targets present in the object.
marker_feature_renames <- c(
  "ENSSSCG00000039214" = "IL1B",
  "ENSSSCG00000029414" = "FCN1",
  "ENSSSCG00000028461" = "SIRPA",
  "ENSSSCG00000006588" = "S100A9",
  "ENSSSCG00000006350" = "FCGR2B",
  "ENSSSCG00000025644" = "CD1.1 (review)",
  "ENSSSCG00000033183" = "MCEMP1",
  "ENSSSCG00000003236" = "SIGLEC5",
  "ENSSSCG00000034914" = "CD163L1",
  "ENSSSCG00000036618" = "FCGR3A",
  "ENSSSCG00000054409" = "CD3G",
  "ENSSSCG00000026210" = "FCRL5",
  "ENSSSCG00000010077" = "IGLL1 (review)",
  "ENSSSCG00000024759" = "CX3CL1",
  "ENSSSCG00000006286" = "SELE",
  "ENSSSCG00000008898" = "HOPX",
  "ENSSSCG00000006525" = "MUC1",
  "ENSSSCG00000002749" = "HP",
  "ENSSSCG00000008115" = "MAL",
  "ENSSSCG00000040016" = "ADGRE4 (review)",
  "ENSSSCG00000003954" = "CFAP144",
  "ENSSSCG00000003754" = "CYP4A90",
  "ENSSSCG00000037735" = "HCAR2",
  "ENSSSCG00000016184" = "LOC100124654 (review)",
  "ENSSSCG00000035296" = "LOC100158011",
  "ENSSSCG00000015405" = "LOC100511343",
  "ENSSSCG00000016183" = "LOC100515345",
  "ENSSSCG00000008769" = "LOC100520753 (review)",
  "ENSSSCG00000026402" = "LOC100522787",
  "ENSSSCG00000000640" = "LOC100523789",
  "ENSSSCG00000013100" = "LOC100524972",
  "ENSSSCG00000051873" = "LOC100525275 (review)",
  "ENSSSCG00000028944" = "LOC100620627",
  "ENSSSCG00000031640" = "LOC100624077",
  "ENSSSCG00000055794" = "LOC100738720",
  "ENSSSCG00000033721" = "LOC102158035",
  "ENSSSCG00000025687" = "PIGE-108A11.3",
  "ENSSSCG00000037358" = "SAA3",
  "ENSSSCG00000001455" = "SLA-DRB1 (review)"
)
