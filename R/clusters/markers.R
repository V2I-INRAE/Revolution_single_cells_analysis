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
# Owner-approved human EVLP additions (2026-10-14): top-50 signatures per
# human subset filtered to symbols present in annotated_clustered.rds, plus
# rename-map rescues (MCEMP1, S100A9, FCN1, IL1B, FCGR3A, FCGR2B, HOPX, MUC1,
# FCRL5, CD3G). 538 entries across 16 panels, tagged "# EVLP <subset> #<rank>".
# Plus exploratory Mast_like_exploratory panel (43 mast-cell genes, rank-ordered)
# added 2026-10-14 at owner request although no mast-cell cluster is annotated.
# Human AM/MoMac/cMo/ncMo have no one-to-one pig equivalent: cMo/ncMo go to
# Monocytes_inflammatory, AM/MoMac to Macrophages_AM_like, complement C1R/C1S/
# C7/CFD plus FCGR3A to Macrophages_C1Q_high, RETN/S100A9 also to Neutrophils.
# cDC2 goes to cDC1; Club_AT2 (transitional) to Secretory_epithelium (MUC1 also
# in AT2); Endothelial/Stromal genes go to both arterial/venous and Fibroblasts.
# POU2F2 placed in B_cells only (ncMo hit likely ambient). Mast-cell genes
# (31, e.g. KIT/MS4A2/HDC) have no panel and stay out; T_gd_like, ASC,
# Lymphatic_EC, Tuft, Smooth_muscle, pDC and Mixed receive none (no fit).
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
    "CX3CR1", # screenshot Mo-C; absent from object (warn-skip)
    # Human EVLP additions 2026-10-14 (ranks in owner xlsx):
    "S100A8", # EVLP cMo #1
    "TIMP1", # EVLP ncMo #1
    "MCEMP1", # EVLP AM #3
    "S100A9", # EVLP cMo #4
    "PPIF", # EVLP cMo #5
    "FGL2", # EVLP ncMo #6
    "AIF1", # EVLP ncMo #7
    "COTL1", # EVLP ncMo #8
    "VCAN", # EVLP cMo #9
    "NLRP3", # EVLP cMo #10
    "NAMPT", # EVLP ncMo #11
    "THBS1", # EVLP cMo #12
    "TNFRSF1B", # EVLP ncMo #12
    "CFP", # EVLP ncMo #13
    "CD14", # EVLP cMo #14
    "RGS2", # EVLP ncMo #14
    "MXD1", # EVLP cMo #15
    "PLAUR", # EVLP cMo #16
    "MPEG1", # EVLP cMo #17
    "GPR183", # EVLP cMo #18
    "CLEC7A", # EVLP ncMo #20
    "WARS1", # EVLP ncMo #21
    "PTGS2", # EVLP cMo #22
    "PDE4B", # EVLP cMo #25
    "LYN", # EVLP ncMo #25
    "PTPRE", # EVLP cMo #26
    "SAMSN1", # EVLP ncMo #26
    "CYRIA", # EVLP ncMo #29
    "SERPINB9", # EVLP cMo #31
    "TNFSF13B", # EVLP ncMo #32
    "IL1R2", # EVLP cMo #34
    "HES4", # EVLP ncMo #34
    "NINJ1", # EVLP cMo #36
    "AP1S2", # EVLP ncMo #36
    "HIF1A", # EVLP cMo #37
    "MTSS1", # EVLP ncMo #37
    "RILPL2", # EVLP cMo #41
    "ZNF703", # EVLP ncMo #41
    "RASSF5", # EVLP ncMo #42
    "PFKFB3", # EVLP cMo #43
    "CD93", # EVLP cMo #44
    "IRAK3", # EVLP cMo #45
    "DMXL2", # EVLP ncMo #45
    "LYZ", # EVLP cMo #46
    "JARID2", # EVLP cMo #48
    "PLEK", # EVLP ncMo #48
    "BCL2A1", # EVLP cMo #49
    "C19orf38", # EVLP ncMo #49
    "ELL2", # EVLP cMo #50
    "VAMP5" # EVLP ncMo #50
  ),
  Macrophages_AM_like = c(
    "LGALS3", "APOE", "S100A12", "FCN2", "PTGR1",
    "CIDEC", "FCGR1A", "FCN1", "GPNMB", "SIGLEC1", #cd169
    # Dotplot top-five additions, cluster 3 (marker run 44561454):
    "LOC100520753 (review)", "PON1",
    # Owner-approved screenshot-MP additions 2026-10-13:
    "PON2", # cl 3 log2FC 2.21/pct.1 0.513
    "TLR4", # cl 3 pct.1 0.917 (screenshot TLR4hi); also cl 13
    "CD68", # cl 3 log2FC 4.11/pct.1 0.870
    # Human EVLP additions 2026-10-14 (ranks in owner xlsx):
    "INHBA", # EVLP AM #1
    "OLR1", # EVLP MoMac #1
    "MRC1", # EVLP MoMac #3
    "CD14", # EVLP MoMac #4
    "VSIG4", # EVLP MoMac #5
    "IFI30", # EVLP MoMac #6
    "PMP22", # EVLP MoMac #7
    "FN1", # EVLP AM #8
    "LYZ", # EVLP MoMac #8
    "PLAUR", # EVLP MoMac #9
    "MS4A7", # EVLP MoMac #10
    "FABP4", # EVLP AM #11
    "ACP5", # EVLP MoMac #11
    "ALOX5AP", # EVLP AM #12
    "CTSB", # EVLP MoMac #13
    "SNX10", # EVLP AM #14
    "FBP1", # EVLP AM #15
    "RGCC", # EVLP AM #17
    "TREM1", # EVLP AM #18
    "FAM89A", # EVLP AM #19
    "LAPTM5", # EVLP MoMac #20
    "RETN", # EVLP MoMac #22
    "GRN", # EVLP AM #23
    "FCGRT", # EVLP MoMac #23
    "SLC7A7", # EVLP MoMac #24
    "CD52", # EVLP AM #25
    "TFRC", # EVLP AM #26
    "CYP27A1", # EVLP AM #28
    "CTSD", # EVLP MoMac #28
    "PPARG", # EVLP AM #29
    "PHLDA1", # EVLP MoMac #29
    "CAPG", # EVLP AM #30
    "CLEC7A", # EVLP MoMac #30
    "CTSC", # EVLP AM #31
    "BHLHE41", # EVLP AM #32
    "CYBB", # EVLP MoMac #32
    "SPI1", # EVLP AM #33
    "CTSS", # EVLP MoMac #33
    "SERPING1", # EVLP AM #34
    "VIM", # EVLP MoMac #34
    "ANPEP", # EVLP MoMac #35
    "NCF2", # EVLP MoMac #37
    "PLA2G7", # EVLP MoMac #38
    "GLDN", # EVLP AM #39
    "CTSA", # EVLP MoMac #39
    "PLXDC2", # EVLP MoMac #40
    "CCDC88A", # EVLP AM #41
    "SLC11A1", # EVLP MoMac #41
    "FHL1", # EVLP AM #42
    "CES1", # EVLP AM #43
    "TYROBP", # EVLP MoMac #43
    "HEXB", # EVLP MoMac #44
    "CTSZ", # EVLP MoMac #45
    "NCEH1", # EVLP AM #46
    "GRINA", # EVLP MoMac #46
    "MSR1", # EVLP AM #47
    "AIF1", # EVLP MoMac #47
    "IGSF6", # EVLP AM #48
    "RNF130", # EVLP MoMac #48
    "PYCARD", # EVLP AM #49
    "KLF4" # EVLP MoMac #49
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
    "IDO1", # screenshot MP; absent from object (warn-skip)
    # Human EVLP additions 2026-10-14 (ranks in owner xlsx):
    "C1R", # EVLP Stromal #10
    "C1S", # EVLP Stromal #12
    "CFD", # EVLP Stromal #13
    "FCGR3A", # EVLP ncMo #33
    "C7" # EVLP Stromal #35
  ),
  cDC1 = c(
    "SLA-DQB1", "SLA-DRB1 (review)", "HLA-DRA", "CD83", "FLT3",
    "XCR1", "TIFAB", "SLC6A12", "DNASE1L3", "RAMP1",
    # Owner-approved screenshot Mo-DC additions 2026-10-13:
    "CD86", # cl 18 log2FC 2.88/pct.1 0.805
    "CCR7", # cl 18 log2FC 2.59/pct.1 0.408
    "CD1E", # screenshot cDC2 comparison; absent from object (warn-skip)
    # Human EVLP additions 2026-10-14 (ranks in owner xlsx):
    "GPR183", # EVLP cDC2 #1
    "NR4A3", # EVLP cDC2 #2
    "PKIB", # EVLP cDC2 #3
    "RGS1", # EVLP cDC2 #6
    "SGK1", # EVLP cDC2 #8
    "TAMALIN", # EVLP cDC2 #9
    "REL", # EVLP cDC2 #13
    "SERPINB9", # EVLP cDC2 #14
    "RGS2", # EVLP cDC2 #15
    "IL1R2", # EVLP cDC2 #16
    "PLEK", # EVLP cDC2 #17
    "CCL17", # EVLP cDC2 #19
    "FGL2", # EVLP cDC2 #21
    "FCER1A", # EVLP cDC2 #23
    "RGS10", # EVLP cDC2 #24
    "JAML", # EVLP cDC2 #25
    "NFKB1", # EVLP cDC2 #27
    "PMAIP1", # EVLP cDC2 #30
    "OGFRL1", # EVLP cDC2 #31
    "PLEKHO1", # EVLP cDC2 #33
    "CIITA", # EVLP cDC2 #35
    "LDLRAD4", # EVLP cDC2 #36
    "MEF2C", # EVLP cDC2 #37
    "SLC8A1", # EVLP cDC2 #40
    "IL13RA1", # EVLP cDC2 #41
    "DUSP4", # EVLP cDC2 #42
    "ENTPD1", # EVLP cDC2 #44
    "PLSCR1", # EVLP cDC2 #45
    "CPVL", # EVLP cDC2 #46
    "DAPP1", # EVLP cDC2 #47
    "CD74", # EVLP cDC2 #48
    "ETV3", # EVLP cDC2 #49
    "PPP1R14A" # EVLP cDC2 #50
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
    "LOC100515345", "STEAP4", "GPR84", "KCNJ15", "G0S2",
    # Human EVLP additions 2026-10-14 (ranks in owner xlsx):
    "S100A9", # EVLP cMo #4
    "RETN" # EVLP MoMac #22
  ),
  Eosinophils = c(
    "RBBP7", "EMB", "MCTP1", "ALOX5", "RGS18",
    "EPX", "LOC100158011", "LOC100522787", "CCR3", "ENSSSCG00000034853"
  ),
  Mast_like_exploratory = c(
    # Human EVLP Mast cells ranks; exploratory panel 2026-10-14, no pig
    # mast-cell cluster annotated. Core identity genes first (ranks 3-10),
    # then remaining ranks in order:
    "HPGDS", # EVLP Mast cells #3
    "GATA2", # EVLP Mast cells #4
    "MS4A2", # EVLP Mast cells #5
    "KIT", # EVLP Mast cells #7
    "IL1RL1", # EVLP Mast cells #9
    "HDC", # EVLP Mast cells #10
    "LTC4S", # EVLP Mast cells #13
    "SLC18A2", # EVLP Mast cells #23
    "PTGS1", # EVLP Mast cells #20
    "HPGD", # EVLP Mast cells #33
    "RGS2", # EVLP Mast cells #6; also Monocytes_inflammatory
    "RGS13", # EVLP Mast cells #8; also Eosinophils
    "CD69", # EVLP Mast cells #11; also T_CD4_naive_CM, T_CD8_cytotoxic
    "RGS1", # EVLP Mast cells #12; also Monocytes_inflammatory, cDC1, T_CD4_naive_CM, T_CD8_cytotoxic
    "BMP2K", # EVLP Mast cells #14
    "CLU", # EVLP Mast cells #15; also Fibroblasts
    "RHEX", # EVLP Mast cells #17; also T_gd_like
    "BATF", # EVLP Mast cells #18; also T_CD4_naive_CM
    "VWA5A", # EVLP Mast cells #19
    "GPR65", # EVLP Mast cells #21
    "GCSAML", # EVLP Mast cells #22
    "ACSL4", # EVLP Mast cells #24
    "FOSB", # EVLP Mast cells #25
    "TMEM233", # EVLP Mast cells #27
    "RAC2", # EVLP Mast cells #28
    "LAT", # EVLP Mast cells #29
    "LMNA", # EVLP Mast cells #31
    "TDRD3", # EVLP Mast cells #32
    "RAB27B", # EVLP Mast cells #34
    "DUSP6", # EVLP Mast cells #35
    "CTTNBP2", # EVLP Mast cells #36
    "TNIK", # EVLP Mast cells #37
    "RGS10", # EVLP Mast cells #38; also Macrophages_AM_like, cDC1
    "SLC26A2", # EVLP Mast cells #39
    "DDIT4", # EVLP Mast cells #40; also T_CD4_naive_CM
    "BIRC3", # EVLP Mast cells #41; also B_cells, T_CD4_naive_CM
    "MAOB", # EVLP Mast cells #42
    "GALC", # EVLP Mast cells #43
    "CALB2", # EVLP Mast cells #45
    "GNPTAB", # EVLP Mast cells #46; also NK_like
    "CTNNBL1", # EVLP Mast cells #48
    "AREG", # EVLP Mast cells #49; also Monocytes_inflammatory
    "ADAM12" # EVLP Mast cells #50
  ),

  # ---- Lymphoid -----------------------------------------------------------
  T_CD4_naive_CM = c(
    "ABLIM1", "CD5", "CD28", "IL7R", "TCF7",
    "CERS6", "ICOS", "LEF1", "CD40LG", "TSPAN5",
    # Human EVLP additions 2026-10-14 (ranks in owner xlsx):
    "LTB", # EVLP CD4_T #2
    "CD69", # EVLP CD4_T #3
    "CXCR4", # EVLP CD4_T #5
    "BCL11B", # EVLP CD4_T #8
    "RGS1", # EVLP CD4_T #11
    "SPOCK2", # EVLP CD4_T #13
    "CD3G", # EVLP CD4_T #15
    "RORA", # EVLP CD4_T #16
    "ZFP36L2", # EVLP CD4_T #17
    "STK4", # EVLP CD4_T #19
    "ITK", # EVLP CD4_T #20
    "IKZF1", # EVLP CD4_T #21
    "RCAN3", # EVLP CD4_T #22
    "CRYBG1", # EVLP CD4_T #24
    "RHOH", # EVLP CD4_T #25
    "GPR171", # EVLP CD4_T #26
    "CNOT6L", # EVLP CD4_T #28
    "BIRC3", # EVLP CD4_T #29
    "DDIT4", # EVLP CD4_T #30
    "BATF", # EVLP CD4_T #32
    "ETS1", # EVLP CD4_T #33
    "SMCHD1", # EVLP CD4_T #34
    "PIK3IP1", # EVLP CD4_T #36
    "ACAP1", # EVLP CD4_T #37
    "FYB1", # EVLP CD4_T #38
    "SOCS1", # EVLP CD4_T #39
    "TRAT1", # EVLP CD4_T #40
    "GPR183", # EVLP CD4_T #41
    "SARAF", # EVLP CD4_T #42
    "SMAP2", # EVLP CD4_T #43
    "CYTIP", # EVLP CD4_T #44
    "LEPROTL1", # EVLP CD4_T #45
    "SYNE2", # EVLP CD4_T #48
    "TSPYL2", # EVLP CD4_T #49
    "ZC3HAV1" # EVLP CD4_T #50
  ),
  T_CD8_cytotoxic = c(
    "CD2", "CCL5", "CD3D", "LOC102158035", "CD3E",
    "CD6", "CD8A", "CXCR6", "CD5", "GZMB",
    # Human EVLP additions 2026-10-14 (ranks in owner xlsx):
    "CXCR4", # EVLP T-CD8 #3
    "NKG7", # EVLP T-CD8 #4
    "CCL4", # EVLP T-CD8 #5
    "PIK3R1", # EVLP T-CD8 #6
    "KLRD1", # EVLP T-CD8 #7
    "CD3G", # EVLP T-CD8 #8
    "CD69", # EVLP T-CD8 #9
    "CST7", # EVLP T-CD8 #10
    "RUNX3", # EVLP T-CD8 #12
    "FYN", # EVLP T-CD8 #13
    "DUSP2", # EVLP T-CD8 #15
    "PRF1", # EVLP T-CD8 #17
    "GZMK", # EVLP T-CD8 #20
    "CNOT6L", # EVLP T-CD8 #22
    "PPP2R5C", # EVLP T-CD8 #24
    "ETS1", # EVLP T-CD8 #25
    "GNLY", # EVLP T-CD8 #29
    "BCL11B", # EVLP T-CD8 #32
    "SLC38A1", # EVLP T-CD8 #33
    "ITGA4", # EVLP T-CD8 #34
    "LCK", # EVLP T-CD8 #35
    "GZMM", # EVLP T-CD8 #36
    "GZMH", # EVLP T-CD8 #37
    "SH2D1A", # EVLP T-CD8 #38
    "ZFP36L2", # EVLP T-CD8 #39
    "RGS1", # EVLP T-CD8 #40
    "ISG20", # EVLP T-CD8 #41
    "CDC42SE2", # EVLP T-CD8 #42
    "PTPRC", # EVLP T-CD8 #43
    "IFNG", # EVLP T-CD8 #44
    "ARL4C", # EVLP T-CD8 #45
    "CTSW", # EVLP T-CD8 #46
    "PTPN22", # EVLP T-CD8 #48
    "PARP8" # EVLP T-CD8 #49
  ),
  T_gd_like = c(
    "CD247", "GATA3", "SYTL3", "RHEX", "TCF7",
    "CD163L1", "SOX13", "BLK", "ENSSSCG00000056250", "SH3BP4"
  ),
  NK_like = c(
    "IL2RB", "KLRK1", "GZMA", "UBD", "TXK",
    "NCR1", "GPR37", "KLRB1", "SEMA6D", "LOC100523789",
    # Human EVLP additions 2026-10-14 (ranks in owner xlsx):
    "GNLY", # EVLP NK #1
    "NKG7", # EVLP NK #3
    "KLRD1", # EVLP NK #4
    "PRF1", # EVLP NK #5
    "KLRF1", # EVLP NK #7
    "CST7", # EVLP NK #10
    "CCL4", # EVLP NK #11
    "CTSW", # EVLP NK #14
    "SPON2", # EVLP NK #16
    "ARL4C", # EVLP NK #17
    "RUNX3", # EVLP NK #20
    "CEP78", # EVLP NK #22
    "GZMM", # EVLP NK #23
    "S1PR5", # EVLP NK #26
    "MYBL1", # EVLP NK #28
    "ABHD17A", # EVLP NK #31
    "PIK3R1", # EVLP NK #33
    "FCGR3A", # EVLP ncMo #33
    "CLIC3", # EVLP NK #34
    "CEMIP2", # EVLP NK #35
    "FYN", # EVLP NK #37
    "PTPN4", # EVLP NK #38
    "GZMH", # EVLP NK #40
    "DUSP2", # EVLP NK #41
    "GNPTAB", # EVLP NK #42
    "CNOT6L", # EVLP NK #43
    "SYNE2", # EVLP NK #44
    "SH2D2A", # EVLP NK #45
    "LIMD2", # EVLP NK #46
    "METRNL", # EVLP NK #47
    "ETS1", # EVLP NK #48
    "SAMD3" # EVLP NK #49
  ),
  B_cells = c(
    "SMC6", "CD79B", "MS4A1", "RALGPS2", "BANK1",
    "CD19", "CLEC17A", "ITGAD", "PAX5", "ELL3",
    # Human EVLP additions 2026-10-14 (ranks in owner xlsx):
    "CD79A", # EVLP B #5
    "IGHM", # EVLP B #9
    "MZB1", # EVLP B #10
    "MEF2C", # EVLP B #12
    "TNFRSF13C", # EVLP B #13
    "LY9", # EVLP B #15
    "HERPUD1", # EVLP B #16
    "BCL11A", # EVLP B #18
    "LTB", # EVLP B #19
    "POU2F2", # EVLP B #20
    "GNG7", # EVLP B #21
    "BIRC3", # EVLP B #22
    "FCMR", # EVLP B #24
    "RHOH", # EVLP B #26
    "TP53INP1", # EVLP B #27
    "ADAM28", # EVLP B #28
    "EZR", # EVLP B #29
    "CXCR5", # EVLP B #30
    "DAPP1", # EVLP B #31
    "NCOA3", # EVLP B #32
    "AFF3", # EVLP B #34
    "CD22", # EVLP B #35
    "CD37", # EVLP B #37
    "RCSD1", # EVLP B #38
    "CD27", # EVLP B #39
    "RUBCNL", # EVLP B #40
    "POU2AF1", # EVLP B #42
    "TNFRSF13B", # EVLP B #43
    "REL", # EVLP B #45
    "FCRL5", # EVLP B #47
    "P2RX5" # EVLP B #48
  ),
  ASC = c(
    "JCHAIN", "ENSSSCG00000036224", "TXNDC5", "PRDM1", "SEC11C",
    "ENSSSCG00000053053", "ENSSSCG00000058302", "ENSSSCG00000063511",
    "ENSSSCG00000048391", "ENSSSCG00000031037"
  ),

  # ---- Endothelium and stroma ---------------------------------------------
  Endothelial_arterial = c(
    "PECAM1", "CLEC14A", "CLDN5", "SEMA3G", "EFNB2",
    "GJA5", "SLC52A3", "SSTR1", "CX3CL1", "HTR2B",
    # Human EVLP additions 2026-10-14 (ranks in owner xlsx):
    "EPAS1", # EVLP Endothelial cells #2
    "CALCRL", # EVLP Endothelial cells #4
    "TIMP3", # EVLP Endothelial cells #6
    "TM4SF1", # EVLP Endothelial cells #7
    "AQP1", # EVLP Endothelial cells #8
    "VWF", # EVLP Endothelial cells #9
    "IGFBP7", # EVLP Endothelial cells #10
    "GNG11", # EVLP Endothelial cells #12
    "HYAL2", # EVLP Endothelial cells #13
    "EGFL7", # EVLP Endothelial cells #14
    "RAMP2", # EVLP Endothelial cells #19
    "TFPI", # EVLP Endothelial cells #20
    "CAV1", # EVLP Endothelial cells #21
    "RAMP3", # EVLP Endothelial cells #22
    "ARHGAP29", # EVLP Endothelial cells #23
    "CDH5", # EVLP Endothelial cells #24
    "BMPR2", # EVLP Endothelial cells #26
    "CAVIN2", # EVLP Endothelial cells #27
    "PALM2AKAP2", # EVLP Endothelial cells #28
    "EMP1", # EVLP Endothelial cells #29
    "SPTBN1", # EVLP Endothelial cells #30
    "CALD1", # EVLP Endothelial cells #31
    "NPDC1", # EVLP Endothelial cells #32
    "PCDH17", # EVLP Endothelial cells #33
    "FLT1", # EVLP Endothelial cells #34
    "APP", # EVLP Endothelial cells #35
    "CAVIN1", # EVLP Endothelial cells #36
    "S1PR1", # EVLP Endothelial cells #37
    "IGFBP4", # EVLP Endothelial cells #38
    "TAMALIN", # EVLP Endothelial cells #39
    "ECSCR", # EVLP Endothelial cells #40
    "SOX18", # EVLP Endothelial cells #41
    "NFIB", # EVLP Endothelial cells #42
    "PTPRB", # EVLP Endothelial cells #43
    "CD93", # EVLP Endothelial cells #45
    "ADGRL4", # EVLP Endothelial cells #46
    "SLCO2A1", # EVLP Endothelial cells #47
    "CRIP2", # EVLP Endothelial cells #48
    "MGP", # EVLP Endothelial cells #49
    "ESAM" # EVLP Endothelial cells #50
  ),
  Endothelial_venous = c(
    "PECAM1", "ADGRF5", "MMRN2", "FMO1", "EMCN",
    "ACKR1", "SLC16A9", "FAM167B", "CCDC85A", "CLEC1A",
    # Human EVLP additions 2026-10-14 (ranks in owner xlsx):
    "EPAS1", # EVLP Endothelial cells #2
    "CALCRL", # EVLP Endothelial cells #4
    "TIMP3", # EVLP Endothelial cells #6
    "TM4SF1", # EVLP Endothelial cells #7
    "AQP1", # EVLP Endothelial cells #8
    "VWF", # EVLP Endothelial cells #9
    "IGFBP7", # EVLP Endothelial cells #10
    "GNG11", # EVLP Endothelial cells #12
    "HYAL2", # EVLP Endothelial cells #13
    "EGFL7", # EVLP Endothelial cells #14
    "RAMP2", # EVLP Endothelial cells #19
    "TFPI", # EVLP Endothelial cells #20
    "CAV1", # EVLP Endothelial cells #21
    "RAMP3", # EVLP Endothelial cells #22
    "ARHGAP29", # EVLP Endothelial cells #23
    "CDH5", # EVLP Endothelial cells #24
    "BMPR2", # EVLP Endothelial cells #26
    "CAVIN2", # EVLP Endothelial cells #27
    "PALM2AKAP2", # EVLP Endothelial cells #28
    "EMP1", # EVLP Endothelial cells #29
    "SPTBN1", # EVLP Endothelial cells #30
    "CALD1", # EVLP Endothelial cells #31
    "NPDC1", # EVLP Endothelial cells #32
    "PCDH17", # EVLP Endothelial cells #33
    "FLT1", # EVLP Endothelial cells #34
    "APP", # EVLP Endothelial cells #35
    "CAVIN1", # EVLP Endothelial cells #36
    "S1PR1", # EVLP Endothelial cells #37
    "IGFBP4", # EVLP Endothelial cells #38
    "TAMALIN", # EVLP Endothelial cells #39
    "ECSCR", # EVLP Endothelial cells #40
    "SOX18", # EVLP Endothelial cells #41
    "NFIB", # EVLP Endothelial cells #42
    "PTPRB", # EVLP Endothelial cells #43
    "CD93", # EVLP Endothelial cells #45
    "ADGRL4", # EVLP Endothelial cells #46
    "SLCO2A1", # EVLP Endothelial cells #47
    "CRIP2", # EVLP Endothelial cells #48
    "MGP", # EVLP Endothelial cells #49
    "ESAM" # EVLP Endothelial cells #50
  ),
  Lymphatic_EC = c(
    "MMRN1", "APOLD1", "PDPN", "CCL21", "MYCT1",
    "HOXD9", "HOXD8", "CHP2", "LYVE1", "FLT4" # LYVE1 replaces SCG3 2026-10-13: canonical lymphatic marker, c12 pct.1 0.772/pct.2 0.02
  ),
  Fibroblasts = c(
    "COL1A1", "COL3A1", "DCN", "COL1A2", "FBLN1",
    "ENSSSCG00000031292", "PI16", "DPT", "LUM", "CHRDL1",
    # Human EVLP additions 2026-10-14 (ranks in owner xlsx):
    "MGP", # EVLP Stromal #2
    "IGFBP7", # EVLP Stromal #6
    "COL6A2", # EVLP Stromal #8
    "CALD1", # EVLP Stromal #9
    "C1R", # EVLP Stromal #10
    "TIMP3", # EVLP Stromal #11
    "C1S", # EVLP Stromal #12
    "CFD", # EVLP Stromal #13
    "RARRES2", # EVLP Stromal #14
    "CCDC80", # EVLP Stromal #15
    "PLAC9", # EVLP Stromal #16
    "SERPINF1", # EVLP Stromal #17
    "MMP2", # EVLP Stromal #18
    "COL6A1", # EVLP Stromal #19
    "IGFBP6", # EVLP Stromal #20
    "AKAP12", # EVLP Stromal #21
    "CLU", # EVLP Stromal #22
    "MFAP4", # EVLP Stromal #23
    "FSTL1", # EVLP Stromal #24
    "FBN1", # EVLP Stromal #25
    "AEBP1", # EVLP Stromal #26
    "TIMP1", # EVLP Stromal #27
    "NNMT", # EVLP Stromal #28
    "CCN1", # EVLP Stromal #29
    "SELENOM", # EVLP Stromal #30
    "A2M", # EVLP Stromal #31
    "CRISPLD2", # EVLP Stromal #32
    "PRRX1", # EVLP Stromal #33
    "SPARC", # EVLP Stromal #34
    "C7", # EVLP Stromal #35
    "CRYAB", # EVLP Stromal #36
    "CCN2", # EVLP Stromal #39
    "PMP22", # EVLP Stromal #40
    "CCN5", # EVLP Stromal #43
    "COL6A3", # EVLP Stromal #44
    "TPM1", # EVLP Stromal #45
    "PCOLCE", # EVLP Stromal #46
    "IGFBP4", # EVLP Stromal #47
    "LTBP4", # EVLP Stromal #48
    "TPM2", # EVLP Stromal #49
    "NBL1" # EVLP Stromal #50
  ),
  Smooth_muscle = c(
    "SPARCL1", "TAGLN", "NDNF", "MYH11", "WIF1",
    "GREM1", "CNN1", "ACTG2", "ENSSSCG00000058753", "ACTA2"
  ),

  # ---- Epithelium ---------------------------------------------------------
  Secretory_epithelium = c(
    "CBR2", "TFF3", "WFDC2", "SULT1E1", "PIGR",
    "SCGB1A1", "SCGB3A1", "BPIFB1", "HSD17B2", "ENSSSCG00000013714",
    # Human EVLP additions 2026-10-14 (ranks in owner xlsx):
    "SCGB3A2", # EVLP Club_AT2 (transitional) #2
    "GDF15", # EVLP Club_AT2 (transitional) #5
    "ELF3", # EVLP Club_AT2 (transitional) #6
    "GPRC5A", # EVLP Club_AT2 (transitional) #7
    "SOX4", # EVLP Club_AT2 (transitional) #8
    "SFTPB", # EVLP Club_AT2 (transitional) #9
    "CP", # EVLP Club_AT2 (transitional) #11
    "TACSTD2", # EVLP Club_AT2 (transitional) #12
    "TCIM", # EVLP Club_AT2 (transitional) #13
    "SLPI", # EVLP Club_AT2 (transitional) #15
    "NCOA7", # EVLP Club_AT2 (transitional) #17
    "KRT18", # EVLP Club_AT2 (transitional) #18
    "TM4SF1", # EVLP Club_AT2 (transitional) #19
    "CLDN4", # EVLP Club_AT2 (transitional) #20
    "KRT7", # EVLP Club_AT2 (transitional) #21
    "AGR2", # EVLP Club_AT2 (transitional) #24
    "KLK11", # EVLP Club_AT2 (transitional) #25
    "KRT19", # EVLP Club_AT2 (transitional) #26
    "KRT8", # EVLP Club_AT2 (transitional) #27
    "EHF", # EVLP Club_AT2 (transitional) #29
    "LCN2", # EVLP Club_AT2 (transitional) #30
    "DSP", # EVLP Club_AT2 (transitional) #31
    "PLPP2", # EVLP Club_AT2 (transitional) #33
    "CLIC6", # EVLP Club_AT2 (transitional) #36
    "MUC1", # EVLP Club_AT2 (transitional) #39
    "AQP5", # EVLP Club_AT2 (transitional) #40
    "ST6GALNAC1", # EVLP Club_AT2 (transitional) #41
    "TSPAN8", # EVLP Club_AT2 (transitional) #43
    "MMP7", # EVLP Club_AT2 (transitional) #44
    "TMEM45A" # EVLP Club_AT2 (transitional) #45
  ),
  Ciliated = c(
    "DNAH11", "CIMIP1", "RSPH4A", "C1orf87", "RSPH1",
    "SNTN", "UBXN10", "CIBAR2", "CFAP144", "CIMAP1B", "PIFO",
    # Human EVLP additions 2026-10-14 (ranks in owner xlsx):
    "CAPS", # EVLP Ciliated #1
    "TSPAN1", # EVLP Ciliated #4
    "TPPP3", # EVLP Ciliated #5
    "LRRIQ1", # EVLP Ciliated #9
    "MORN2", # EVLP Ciliated #11
    "CETN2", # EVLP Ciliated #15
    "DNAAF1", # EVLP Ciliated #16
    "ROPN1L", # EVLP Ciliated #19
    "ZMYND10", # EVLP Ciliated #20
    "SMIM22", # EVLP Ciliated #21
    "CCDC170", # EVLP Ciliated #23
    "CAPSL", # EVLP Ciliated #24
    "EFHC1", # EVLP Ciliated #26
    "IGFBP7", # EVLP Ciliated #27
    "CEP126", # EVLP Ciliated #28
    "FAM229B", # EVLP Ciliated #30
    "RP1", # EVLP Ciliated #31
    "PLPP2", # EVLP Ciliated #32
    "DYNLRB2", # EVLP Ciliated #34
    "SPA17", # EVLP Ciliated #35
    "SLC44A4", # EVLP Ciliated #38
    "STK33", # EVLP Ciliated #39
    "IGFBP2", # EVLP Ciliated #41
    "DTHD1", # EVLP Ciliated #42
    "PPIL6", # EVLP Ciliated #43
    "CCDC78", # EVLP Ciliated #44
    "DNALI1", # EVLP Ciliated #45
    "DNAH5", # EVLP Ciliated #47
    "SPAG17", # EVLP Ciliated #49
    "CDHR3" # EVLP Ciliated #50
  ),
  Tuft = c(
    "SH2D7", "POU2F3", "RGS21", "RGS13", "AVIL",
    "HMX3", "CYP4A90", "HMX2", "GNG13", "IHH"
  ),
  AT2 = c(
    "SFTPA1", "SFTPC", "SFTPD", "LAMP3", "NAPSA",
    "GJB1", "IYD", "C5", "CHIA", "SLC34A2",
    # Human EVLP additions 2026-10-14 (ranks in owner xlsx):
    "SFTPB", # EVLP AT2 #4
    "SLPI", # EVLP AT2 #7
    "MUC1", # EVLP AT2 #8
    "SFTA2", # EVLP AT2 #9
    "S100A14", # EVLP AT2 #13
    "PEBP4", # EVLP AT2 #16
    "TCIM", # EVLP AT2 #17
    "RNASE1", # EVLP AT2 #20
    "LRRK2", # EVLP AT2 #21
    "KRT18", # EVLP AT2 #22
    "DHCR24", # EVLP AT2 #23
    "SLC39A8", # EVLP AT2 #25
    "LPCAT1", # EVLP AT2 #26
    "CTSH", # EVLP AT2 #27
    "MALL", # EVLP AT2 #28
    "EPCAM", # EVLP AT2 #29
    "GPRC5A", # EVLP AT2 #30
    "ELF3", # EVLP AT2 #31
    "NEDD4L", # EVLP AT2 #32
    "KRT19", # EVLP AT2 #33
    "ABCA3", # EVLP AT2 #34
    "CLDN4", # EVLP AT2 #35
    "KRT8", # EVLP AT2 #36
    "C3", # EVLP AT2 #37
    "C16orf89", # EVLP AT2 #38
    "CYB5A", # EVLP AT2 #39
    "FASN", # EVLP AT2 #41
    "MAL2", # EVLP AT2 #44
    "TACSTD2", # EVLP AT2 #45
    "MTUS1", # EVLP AT2 #47
    "SDR16C5", # EVLP AT2 #48
    "C4BPA", # EVLP AT2 #49
    "TMC5" # EVLP AT2 #50
  ),
  AT1 = c(
    "AGER", "HOPX", "LMO7", "PDPN", "ADIRF",
    "LOC100738720", "COL4A3", "COL4A4", "VSIG2", "FLRT3",
    # Human EVLP additions 2026-10-14 (ranks in owner xlsx):
    "EMP2", # EVLP AT1 #2
    "CAV1", # EVLP AT1 #3
    "KRT7", # EVLP AT1 #4
    "CAV2", # EVLP AT1 #8
    "KRT19", # EVLP AT1 #9
    "GPRC5A", # EVLP AT1 #10
    "SCEL", # EVLP AT1 #11
    "TSPAN13", # EVLP AT1 #13
    "AQP4", # EVLP AT1 #14
    "TACSTD2", # EVLP AT1 #15
    "CLDN18", # EVLP AT1 #16
    "CLIC5", # EVLP AT1 #17
    "KRT18", # EVLP AT1 #18
    "CLIC3", # EVLP AT1 #19
    "GJA1", # EVLP AT1 #22
    "SLC39A8", # EVLP AT1 #23
    "RTKN2", # EVLP AT1 #24
    "LIMCH1", # EVLP AT1 #26
    "PLS3", # EVLP AT1 #27
    "ANKRD29", # EVLP AT1 #28
    "PALM2AKAP2", # EVLP AT1 #30
    "EPCAM", # EVLP AT1 #31
    "GAS6", # EVLP AT1 #32
    "MYL9", # EVLP AT1 #33
    "CADM1", # EVLP AT1 #34
    "MYO1B", # EVLP AT1 #35
    "KRT8", # EVLP AT1 #36
    "GALNT5", # EVLP AT1 #37
    "DSTN", # EVLP AT1 #38
    "SPTBN1", # EVLP AT1 #39
    "ARHGEF26", # EVLP AT1 #40
    "SFTA2", # EVLP AT1 #41
    "MGLL", # EVLP AT1 #42
    "CKB", # EVLP AT1 #43
    "TIMP3", # EVLP AT1 #44
    "TNFRSF12A", # EVLP AT1 #45
    "WWC2", # EVLP AT1 #46
    "LAMB3", # EVLP AT1 #47
    "NEDD9", # EVLP AT1 #48
    "CD151", # EVLP AT1 #49
    "COL8A1" # EVLP AT1 #50
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
