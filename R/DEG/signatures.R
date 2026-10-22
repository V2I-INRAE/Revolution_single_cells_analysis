# Source-human signatures for future intercluster scoring; no pig mapping or scores.
# Owner-approved reviewed response sets, frozen 2026-10-22.
# This file defines data only. It is deliberately NOT sourced by the active DEG
# pipeline; adding it does not alter existing DE, GSEA, diagnostics or figures.
#
# Each entry has positive and negative character vectors. Positive means higher
# expression in the named response (or the catalog membership for mixed pathway
# sets), not that every member activates the pathway. Negative is used only where
# the source explicitly supplies a downregulated arm (Rooney NRF2).
# Preserve the arms: do not flatten signed signatures to an unsigned gene set.
# No score formula, weights, eligibility cutoff or missing-gene policy is chosen.
#
# MSigDB lists are EXACT human-symbol GMT memberships from 2026.1.Hs, matching the
# existing GSEA release, not the live 2026.2.Hs catalog. Original publication
# spellings are retained for Huang/Rooney; aliases require a later explicit audit.
# Source IDs are human, NOT validated pig orthologues. Pig mapping and coverage
# checks are required before any scoring. Intercluster differences can reflect
# cell identity as well as state; overlapping modules are not independent tests.
#
# MSigDB: CC-BY 4.0 with source-specific exceptions; cite source publications.
# Original paper lists retain their publication attribution/license.
# Verification and source-file SHA256: logs/signature-registry-20261022.INFO.
# membership_sha256 hashes UTF-8: "positive\n", each positive symbol + "\n",
# then "negative\n", each negative symbol + "\n", preserving listed order.
#
# Not included: unresolved ATF4/ISR, endothelial/epithelial injury, dissociation
# lists, GO mechanical annotations or mixed death-regulator catalogs. No genes
# were invented to fill these provisional modules.

signatures <- list(
  HALLMARK_INTERFERON_ALPHA_RESPONSE = list(
    positive = c(
      "ADAR", "B2M", "BATF2", "BST2", "C1S", "CASP1", "CASP8", "CCRL2", "CD47", "CD74",
      "CMPK2", "CMTR1", "CNP", "CSF1", "CXCL10", "CXCL11", "DDX60", "DHX58", "EIF2AK2",
      "ELF1", "EPSTI1", "GBP2", "GBP4", "GMPR", "HELZ2", "HERC6", "HLA-C", "IFI27", "IFI30",
      "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT2", "IFIT3", "IFITM1", "IFITM2", "IFITM3",
      "IL15", "IL4R", "IL7", "IRF1", "IRF2", "IRF7", "IRF9", "ISG15", "ISG20", "LAMP3",
      "LAP3", "LGALS3BP", "LPAR6", "LY6E", "MOV10", "MVB12A", "MX1", "NCOA7", "NMI", "NUB1",
      "OAS1", "OASL", "OGFR", "PARP12", "PARP14", "PARP9", "PLSCR1", "PNPT1", "PROCR",
      "PSMA3", "PSMB8", "PSMB9", "PSME1", "PSME2", "RIPK2", "RNF31", "RSAD2", "RTP4",
      "SAMD9", "SAMD9L", "SELL", "SLC25A28", "SP110", "STAT2", "TAP1", "TDRD7", "TENT5A",
      "TMEM140", "TRAFD1", "TRIM14", "TRIM21", "TRIM25", "TRIM26", "TRIM5", "TXNIP", "UBA7",
      "UBE2L6", "USP18", "WARS1"
    ),
    negative = character()
  ),
  HALLMARK_INTERFERON_GAMMA_RESPONSE = list(
    positive = c(
      "ADAR", "APOL6", "ARID5B", "ARL4A", "AUTS2", "B2M", "BANK1", "BATF2", "BPGM", "BST2",
      "BTG1", "C1R", "C1S", "CASP1", "CASP3", "CASP4", "CASP7", "CASP8", "CCL2", "CCL5",
      "CCL7", "CD274", "CD38", "CD40", "CD69", "CD74", "CD86", "CDKN1A", "CFB", "CFH",
      "CIITA", "CMKLR1", "CMPK2", "CMTR1", "CSF2RB", "CXCL10", "CXCL11", "CXCL9", "DDX60",
      "DHX58", "EIF2AK2", "EIF4E3", "EPSTI1", "FAS", "FCGR1A", "FGL2", "FPR1", "GBP4",
      "GBP6", "GCH1", "GPR18", "GZMA", "HELZ2", "HERC6", "HIF1A", "HLA-A", "HLA-B",
      "HLA-DMA", "HLA-DQA1", "HLA-DRB1", "HLA-G", "ICAM1", "IDO1", "IFI27", "IFI30",
      "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT1", "IFIT2", "IFIT3", "IFITM2", "IFITM3",
      "IFNAR2", "IL10RA", "IL15", "IL15RA", "IL18BP", "IL2RB", "IL4R", "IL6", "IL7", "IRF1",
      "IRF2", "IRF4", "IRF5", "IRF7", "IRF8", "IRF9", "ISG15", "ISG20", "ISOC1", "ITGB7",
      "JAK2", "KLRK1", "LAP3", "LATS2", "LCP2", "LGALS3BP", "LY6E", "LYSMD2", "MARCHF1",
      "MT2A", "MTHFD2", "MVP", "MX1", "MX2", "MYD88", "NAMPT", "NCOA3", "NFKB1", "NFKBIA",
      "NLRC5", "NMI", "NOD1", "NUP93", "OAS2", "OAS3", "OASL", "OGFR", "P2RY14", "PARP12",
      "PARP14", "PDE4B", "PELI1", "PFKP", "PIM1", "PLA2G4A", "PLSCR1", "PML", "PNP",
      "PNPT1", "PSMA2", "PSMA3", "PSMB10", "PSMB2", "PSMB8", "PSMB9", "PSME1", "PSME2",
      "PTGS2", "PTPN1", "PTPN2", "PTPN6", "RAPGEF6", "RBCK1", "RIGI", "RIPK1", "RIPK2",
      "RNF213", "RNF31", "RSAD2", "RTP4", "SAMD9L", "SAMHD1", "SECTM1", "SELP", "SERPING1",
      "SLAMF7", "SLC25A28", "SOCS1", "SOCS3", "SOD2", "SP110", "SPPL2A", "SRI", "SSPN",
      "ST3GAL5", "ST8SIA4", "STAT1", "STAT2", "STAT3", "STAT4", "TAP1", "TAPBP", "TDRD7",
      "TMT1B", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TNFSF10", "TOR1B", "TRAFD1", "TRIM14",
      "TRIM21", "TRIM25", "TRIM26", "TXNIP", "UBE2L6", "UPP1", "USP18", "VAMP5", "VAMP8",
      "VCAM1", "WARS1", "XAF1", "XCL1", "ZBP1", "ZNFX1"
    ),
    negative = character()
  ),
  HALLMARK_HYPOXIA = list(
    positive = c(
      "ACKR3", "ADM", "ADORA2B", "AK4", "AKAP12", "ALDOA", "ALDOB", "ALDOC", "AMPD3",
      "ANGPTL4", "ANKZF1", "ANXA2", "ATF3", "ATP7A", "B3GALT6", "B4GALNT2", "BCAN", "BCL2",
      "BGN", "BHLHE40", "BNIP3L", "BRS3", "BTG1", "CA12", "CASP6", "CAV1", "CAVIN1",
      "CAVIN3", "CCN1", "CCN2", "CCN5", "CCNG2", "CDKN1A", "CDKN1B", "CDKN1C", "CHST2",
      "CHST3", "CITED2", "COL5A1", "CP", "CSRP2", "CXCR4", "DCN", "DDIT3", "DDIT4",
      "DPYSL4", "DTNA", "DUSP1", "EDN2", "EFNA1", "EFNA3", "EGFR", "ENO1", "ENO2", "ENO3",
      "ERO1A", "ERRFI1", "ETS1", "EXT1", "F3", "FAM162A", "FBP1", "FOS", "FOSL2", "FOXO3",
      "GAA", "GALK1", "GAPDH", "GAPDHS", "GBE1", "GCK", "GCNT2", "GLRX", "GPC1", "GPC3",
      "GPC4", "GPI", "GRHPR", "GYS1", "HAS1", "HDLBP", "HEXA", "HK1", "HK2", "HMOX1",
      "HOXB9", "HS3ST1", "HSPA5", "IDS", "IER3", "IGFBP1", "IGFBP3", "IL6", "ILVBL", "INHA",
      "IRS2", "ISG20", "JMJD6", "JUN", "KDELR3", "KDM3A", "KIF5A", "KLF6", "KLF7", "KLHL24",
      "LALBA", "LARGE1", "LDHA", "LDHC", "LOX", "LXN", "MAFF", "MAP3K1", "MIF", "MT1E",
      "MT2A", "MXI1", "MYH9", "NAGK", "NCAN", "NDRG1", "NDST1", "NDST2", "NEDD4L", "NFIL3",
      "NOCT", "NR3C1", "P4HA1", "P4HA2", "PAM", "PCK1", "PDGFB", "PDK1", "PDK3", "PFKFB3",
      "PFKL", "PFKP", "PGAM2", "PGF", "PGK1", "PGM1", "PGM2", "PHKG1", "PIM1", "PKLR",
      "PKP1", "PLAC8", "PLAUR", "PLIN2", "PNRC1", "PPARGC1A", "PPFIA4", "PPP1R15A",
      "PPP1R3C", "PRDX5", "PRKCA", "PYGM", "RBPJ", "RORA", "RRAGD", "S100A4", "SAP30",
      "SCARB1", "SDC2", "SDC3", "SDC4", "SELENBP1", "SERPINE1", "SIAH2", "SLC25A1",
      "SLC2A1", "SLC2A3", "SLC2A5", "SLC37A4", "SLC6A6", "SRPX", "STBD1", "STC1", "STC2",
      "SULT2B1", "TES", "TGFB3", "TGFBI", "TGM2", "TIPARP", "TKTL1", "TMEM45A", "TNFAIP3",
      "TPBG", "TPD52", "TPI1", "TPST2", "UGP2", "VEGFA", "VHL", "VLDLR", "WSB1", "XPNPEP1",
      "ZFP36", "ZNF292"
    ),
    negative = character()
  ),
  HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY = list(
    positive = c(
      "ABCC1", "ATOX1", "CAT", "CDKN2D", "EGLN2", "ERCC2", "FES", "FTL", "G6PD", "GCLC",
      "GCLM", "GLRX", "GLRX2", "GPX3", "GPX4", "GSR", "HHEX", "HMOX2", "IPCEF1", "JUNB",
      "LAMTOR5", "LSP1", "MBP", "MGST1", "MPO", "MSRA", "NDUFA6", "NDUFB4", "NDUFS2",
      "NQO1", "OXSR1", "PDLIM1", "PFKP", "PRDX1", "PRDX2", "PRDX4", "PRDX6", "PRNP", "PTPA",
      "SBNO2", "SCAF4", "SELENOS", "SOD1", "SOD2", "SRXN1", "STK25", "TXN", "TXNRD1",
      "TXNRD2"
    ),
    negative = character()
  ),
  HALLMARK_UNFOLDED_PROTEIN_RESPONSE = list(
    positive = c(
      "ALDH18A1", "ARFGAP1", "ASNS", "ATF3", "ATF4", "ATF6", "ATP6V0D1", "BAG3", "BANF1",
      "CALR", "CCL2", "CEBPB", "CEBPG", "CHAC1", "CKS1B", "CNOT2", "CNOT4", "CNOT6",
      "CXXC1", "DCP1A", "DCP2", "DCTN1", "DDIT4", "DDX10", "DKC1", "DNAJA4", "DNAJB9",
      "DNAJC3", "EDC4", "EDEM1", "EEF2", "EIF2AK3", "EIF2S1", "EIF4A1", "EIF4A2", "EIF4A3",
      "EIF4E", "EIF4EBP1", "EIF4G1", "ERN1", "ERO1A", "EXOC2", "EXOSC1", "EXOSC10",
      "EXOSC2", "EXOSC4", "EXOSC5", "EXOSC9", "FKBP14", "FUS", "GEMIN4", "GOSR2", "H2AX",
      "HERPUD1", "HSP90B1", "HSPA5", "HSPA9", "HYOU1", "IARS1", "IFIT1", "IGFBP1", "IMP3",
      "KDELR3", "KHSRP", "KIF5B", "LSM1", "LSM4", "MTHFD2", "MTREX", "NABP1", "NFYA",
      "NFYB", "NHP2", "NOLC1", "NOP14", "NOP56", "NPM1", "PAIP1", "PARN", "PDIA5", "PDIA6",
      "POP4", "PREB", "PSAT1", "RPS14", "RRP9", "SDAD1", "SEC11A", "SEC31A", "SERP1",
      "SHC1", "SKIC3", "SLC1A4", "SLC30A5", "SLC7A5", "SPCS1", "SPCS3", "SRPRA", "SRPRB",
      "SSR1", "STC2", "TARS1", "TATDN2", "TSPYL2", "TUBB2A", "VEGFA", "WFS1", "WIPI1",
      "XBP1", "XPOT", "YIF1A", "YWHAZ", "ZBTB17"
    ),
    negative = character()
  ),
  HALLMARK_TNFA_SIGNALING_VIA_NFKB = list(
    positive = c(
      "ABCA1", "ACKR3", "AREG", "ATF3", "ATP2B1", "B4GALT1", "B4GALT5", "BCL2A1", "BCL3",
      "BCL6", "BHLHE40", "BIRC2", "BIRC3", "BMP2", "BTG1", "BTG2", "BTG3", "CCL2", "CCL20",
      "CCL4", "CCL5", "CCN1", "CCND1", "CCNL1", "CCRL2", "CD44", "CD69", "CD80", "CD83",
      "CDKN1A", "CEBPB", "CEBPD", "CFLAR", "CLCF1", "CSF1", "CSF2", "CXCL1", "CXCL10",
      "CXCL11", "CXCL2", "CXCL3", "CXCL6", "DENND5A", "DNAJB4", "DRAM1", "DUSP1", "DUSP2",
      "DUSP4", "DUSP5", "EDN1", "EFNA1", "EGR1", "EGR2", "EGR3", "EHD1", "EIF1", "ETS2",
      "F2RL1", "F3", "FJX1", "FOS", "FOSB", "FOSL1", "FOSL2", "FUT4", "G0S2", "GADD45A",
      "GADD45B", "GCH1", "GEM", "GFPT2", "GPR183", "HBEGF", "HES1", "ICAM1", "ICOSLG",
      "ID2", "IER2", "IER3", "IER5", "IFIH1", "IFIT2", "IFNGR2", "IL12B", "IL15RA", "IL18",
      "IL1A", "IL1B", "IL23A", "IL6", "IL6ST", "IL7R", "INHBA", "IRF1", "IRS2", "JAG1",
      "JUN", "JUNB", "KDM6B", "KLF10", "KLF2", "KLF4", "KLF6", "KLF9", "KYNU", "LAMB3",
      "LDLR", "LIF", "LITAF", "MAFF", "MAP2K3", "MAP3K8", "MARCKS", "MCL1", "MSC", "MXD1",
      "MYC", "NAMPT", "NFAT5", "NFE2L2", "NFIL3", "NFKB1", "NFKB2", "NFKBIA", "NFKBIE",
      "NINJ1", "NR4A1", "NR4A2", "NR4A3", "OLR1", "PANX1", "PDE4B", "PDLIM5", "PER1",
      "PFKFB3", "PHLDA1", "PHLDA2", "PLAU", "PLAUR", "PLEK", "PLK2", "PLPP3", "PMEPA1",
      "PNRC1", "PPP1R15A", "PTGER4", "PTGS2", "PTPRE", "PTX3", "RCAN1", "REL", "RELA",
      "RELB", "RHOB", "RIGI", "RIPK2", "RNF19B", "SAT1", "SDC4", "SERPINB2", "SERPINB8",
      "SERPINE1", "SGK1", "SIK1", "SLC16A6", "SLC2A3", "SLC2A6", "SMAD3", "SNN", "SOCS3",
      "SOD2", "SPHK1", "SPSB1", "SQSTM1", "STAT5A", "TANK", "TAP1", "TGIF1", "TIPARP",
      "TLR2", "TNC", "TNF", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TNFAIP8", "TNFRSF9", "TNFSF9",
      "TNIP1", "TNIP2", "TRAF1", "TRIB1", "TRIP10", "TSC22D1", "TUBB2A", "VEGFA", "YRDC",
      "ZBTB10", "ZC3H12A", "ZFP36"
    ),
    negative = character()
  ),
  HALLMARK_INFLAMMATORY_RESPONSE = list(
    positive = c(
      "ABCA1", "ABI1", "ACVR1B", "ACVR2A", "ADGRE1", "ADM", "ADORA2B", "ADRM1", "AHR",
      "APLNR", "AQP9", "ATP2A2", "ATP2B1", "ATP2C1", "AXL", "BDKRB1", "BEST1", "BST2",
      "BTG2", "C3AR1", "C5AR1", "CALCRL", "CCL17", "CCL2", "CCL20", "CCL22", "CCL24",
      "CCL5", "CCL7", "CCR7", "CCRL2", "CD14", "CD40", "CD48", "CD55", "CD69", "CD70",
      "CD82", "CDKN1A", "CHST2", "CLEC5A", "CMKLR1", "CSF1", "CSF3", "CSF3R", "CX3CL1",
      "CXCL10", "CXCL11", "CXCL6", "CXCL8", "CXCL9", "CXCR6", "CYBB", "DCBLD2", "EBI3",
      "EDN1", "EIF2AK2", "EMP3", "EREG", "F3", "FFAR2", "FPR1", "FZD5", "GABBR1", "GCH1",
      "GNA15", "GNAI3", "GP1BA", "GPC3", "GPR132", "GPR183", "HAS2", "HBEGF", "HIF1A",
      "HPN", "HRH1", "ICAM1", "ICAM4", "ICOSLG", "IFITM1", "IFNAR1", "IFNGR2", "IL10",
      "IL10RA", "IL12B", "IL15", "IL15RA", "IL18", "IL18R1", "IL18RAP", "IL1A", "IL1B",
      "IL1R1", "IL2RB", "IL4R", "IL6", "IL7R", "INHBA", "IRAK2", "IRF1", "IRF7", "ITGA5",
      "ITGB3", "ITGB8", "KCNA3", "KCNJ2", "KCNMB2", "KIF1B", "KLF6", "LAMP3", "LCK", "LCP2",
      "LDLR", "LIF", "LPAR1", "LTA", "LY6E", "LYN", "MARCO", "MEFV", "MEP1A", "MET",
      "MMP14", "MSR1", "MXD1", "MYC", "NAMPT", "NDP", "NFKB1", "NFKBIA", "NLRP3", "NMI",
      "NMUR1", "NOD2", "NPFFR2", "OLR1", "OPRK1", "OSM", "OSMR", "P2RX4", "P2RX7", "P2RY2",
      "PCDH7", "PDE4B", "PDPN", "PIK3R5", "PLAUR", "PROK2", "PSEN1", "PTAFR", "PTGER2",
      "PTGER4", "PTGIR", "PTPRE", "PVR", "RAF1", "RASGRP1", "RELA", "RGS1", "RGS16", "RHOG",
      "RIPK2", "RNF144B", "ROS1", "RTP4", "SCARF1", "SCN1B", "SELE", "SELENOS", "SELL",
      "SEMA4D", "SERPINE1", "SGMS2", "SLAMF1", "SLC11A2", "SLC1A2", "SLC28A2", "SLC31A1",
      "SLC31A2", "SLC4A4", "SLC7A1", "SLC7A2", "SPHK1", "SRI", "STAB1", "TACR1", "TACR3",
      "TAPBP", "TIMP1", "TLR1", "TLR2", "TLR3", "TNFAIP6", "TNFRSF1B", "TNFRSF9", "TNFSF10",
      "TNFSF15", "TNFSF9", "TPBG", "VIP"
    ),
    negative = character()
  ),
  HALLMARK_P53_PATHWAY = list(
    positive = c(
      "ABAT", "ABCC5", "ABHD4", "ACVR1B", "ADA", "AEN", "AK1", "ALOX15B", "ANKRA2", "APAF1",
      "APP", "ATF3", "BAIAP2", "BAK1", "BAX", "BLCAP", "BMP2", "BTG1", "BTG2", "CASP1",
      "CCND2", "CCND3", "CCNG1", "CCNK", "CCP110", "CD81", "CD82", "CDH13", "CDK5R1",
      "CDKN1A", "CDKN2A", "CDKN2AIP", "CDKN2B", "CEBPA", "CGRRF1", "CLCA2", "COQ8A",
      "CSRNP2", "CTSD", "CTSF", "CYFIP2", "DCXR", "DDB2", "DDIT3", "DDIT4", "DEF6", "DGKA",
      "DNTTIP2", "DRAM1", "EI24", "ELP1", "EPHA2", "EPHX1", "EPS8L2", "ERCC5", "F2R",
      "FAM162A", "FAS", "FBXW7", "FDXR", "FGF13", "FOS", "FOXO3", "FUCA1", "GADD45A",
      "GLS2", "GM2A", "GPX2", "H1-2", "H2AC25", "H2AJ", "HBEGF", "HDAC3", "HEXIM1", "HINT1",
      "HMOX1", "HRAS", "HSPA4L", "IER3", "IER5", "IFI30", "IL1A", "INHBB", "IP6K2", "IRAG2",
      "IRAK1", "ISCU", "ITGB4", "JAG2", "JUN", "KIF13B", "KLF4", "KLK8", "KRT17", "LDHB",
      "LIF", "MAPKAPK3", "MDM2", "MKNK2", "MXD1", "MXD4", "NDRG1", "NHLH2", "NINJ1", "NOL8",
      "NOTCH1", "NUDT15", "NUPR1", "OSGIN1", "PCNA", "PDGFA", "PERP", "PHLDA3", "PIDD1",
      "PITPNC1", "PLK2", "PLK3", "PLXNB2", "PMM1", "POLH", "POM121", "PPM1D", "PPP1R15A",
      "PRKAB1", "PRMT2", "PROCR", "PTPN14", "PTPRE", "PVT1", "RAB40C", "RACK1", "RAD51C",
      "RAD9A", "RALGDS", "RAP2B", "RB1", "RCHY1", "RETSAT", "RGS16", "RHBDF2", "RNF19B",
      "RPL18", "RPL36", "RPS12", "RPS27L", "RRAD", "RRP8", "RXRA", "S100A10", "S100A4",
      "SAT1", "SDC1", "SEC61A1", "SERPINB5", "SERTAD3", "SESN1", "SFN", "SLC19A2",
      "SLC35D1", "SLC3A2", "SLC7A11", "SOCS1", "SP1", "SPHK1", "ST14", "STEAP3", "STOM",
      "TAP1", "TAX1BP3", "TCHH", "TCN2", "TGFA", "TGFB1", "TM4SF1", "TM7SF3", "TNFSF9",
      "TNNI1", "TOB1", "TP53", "TP63", "TPD52L1", "TPRKB", "TRAF4", "TRAFD1", "TRIAP1",
      "TRIB3", "TSC22D1", "TSPYL2", "TXNIP", "UPP1", "VAMP8", "VDR", "VWA5A", "WRAP73",
      "WWP1", "XPC", "ZBTB16", "ZFP36L1", "ZMAT3", "ZNF365"
    ),
    negative = character()
  ),
  HALLMARK_DNA_REPAIR = list(
    positive = c(
      "AAAS", "ADA", "ADCY6", "ADRM1", "AGO4", "AK1", "AK3", "ALYREF", "APRT", "ARL6IP1",
      "BCAM", "BCAP31", "BOLA2", "BRF2", "CANT1", "CCNO", "CDA", "CETN2", "CLP1", "CMPK2",
      "COX17", "CSTF3", "DAD1", "DCTN4", "DDB1", "DDB2", "DGCR8", "DGUOK", "DUT", "EDF1",
      "EIF1B", "ELL", "ELOA", "ERCC1", "ERCC2", "ERCC3", "ERCC4", "ERCC5", "ERCC8", "FEN1",
      "GMPR2", "GPX4", "GSDME", "GTF2A2", "GTF2B", "GTF2F1", "GTF2H1", "GTF2H3", "GTF2H5",
      "GTF3C5", "GUK1", "HCLS1", "HPRT1", "IMPDH2", "ITPA", "LIG1", "MPC2", "MPG", "MRPL40",
      "NCBP2", "NELFB", "NELFCD", "NELFE", "NFX1", "NME1", "NME3", "NME4", "NPR2", "NT5C",
      "NT5C3A", "NUDT21", "NUDT9", "PCNA", "PDE4B", "PDE6G", "PNP", "POLA1", "POLA2",
      "POLB", "POLD1", "POLD3", "POLD4", "POLE4", "POLH", "POLL", "POLR1C", "POLR1D",
      "POLR1H", "POLR2A", "POLR2C", "POLR2D", "POLR2E", "POLR2F", "POLR2G", "POLR2H",
      "POLR2I", "POLR2J", "POLR2K", "POLR3C", "POLR3GL", "POM121", "PRIM1", "RAD51",
      "RAD52", "RAE1", "RALA", "RBX1", "REV3L", "RFC2", "RFC3", "RFC4", "RFC5", "RNMT",
      "RPA2", "RPA3", "RRM2B", "SAC3D1", "SDCBP", "SEC61A1", "SF3A3", "SMAD5", "SNAPC4",
      "SNAPC5", "SRSF6", "SSRP1", "STX3", "SUPT4H1", "SUPT5H", "SURF1", "TAF10", "TAF12",
      "TAF13", "TAF1C", "TAF6", "TAF9", "TARBP2", "TK2", "TMED2", "TP53", "TSG101", "TYMS",
      "UMPS", "UPF3B", "USP11", "VPS28", "VPS37B", "VPS37D", "XPC", "ZNF707", "ZWINT"
    ),
    negative = character()
  ),
  HALLMARK_APOPTOSIS = list(
    positive = c(
      "ADD1", "AIFM3", "ANKH", "ANXA1", "APP", "ATF3", "AVPR1A", "BAX", "BCAP31", "BCL10",
      "BCL2L1", "BCL2L10", "BCL2L11", "BCL2L2", "BGN", "BID", "BIK", "BIRC3", "BMF", "BMP2",
      "BNIP3L", "BRCA1", "BTG2", "BTG3", "CASP1", "CASP2", "CASP3", "CASP4", "CASP6",
      "CASP7", "CASP8", "CASP9", "CAV1", "CCNA1", "CCND1", "CCND2", "CD14", "CD2", "CD38",
      "CD44", "CD69", "CDC25B", "CDK2", "CDKN1A", "CDKN1B", "CFLAR", "CLU", "CREBBP", "CTH",
      "CTNNB1", "CYLD", "DAP", "DAP3", "DCN", "DDIT3", "DFFA", "DIABLO", "DNAJA1", "DNAJC3",
      "DNM1L", "DPYD", "EBP", "EGR3", "EMP1", "ENO2", "ERBB2", "ERBB3", "EREG", "ETF1",
      "F2", "F2R", "FAS", "FASLG", "FDXR", "FEZ1", "GADD45A", "GADD45B", "GCH1", "GNA15",
      "GPX1", "GPX3", "GPX4", "GSN", "GSR", "GSTM1", "GUCY2D", "H1-0", "HGF", "HMGB2",
      "HMOX1", "HSPB1", "IER3", "IFITM3", "IFNB1", "IFNGR1", "IGF2R", "IGFBP6", "IL18",
      "IL1A", "IL1B", "IL6", "IRF1", "ISG20", "JUN", "KRT18", "LEF1", "LGALS3", "LMNA",
      "LUM", "MADD", "MCL1", "MGMT", "MMP2", "NEDD9", "NEFH", "PAK1", "PDCD4", "PDGFRB",
      "PEA15", "PLAT", "PLCB2", "PLPPR4", "PMAIP1", "PPP2R5B", "PPP3R1", "PPT1", "PRF1",
      "PSEN1", "PSEN2", "PTK2", "RARA", "RELA", "RETSAT", "RHOB", "RHOT2", "RNASEL",
      "ROCK1", "SAT1", "SATB1", "SC5D", "SLC20A1", "SMAD7", "SOD1", "SOD2", "SPTAN1",
      "SQSTM1", "TAP1", "TGFB2", "TGFBR3", "TIMP1", "TIMP2", "TIMP3", "TNF", "TNFRSF12A",
      "TNFSF10", "TOP2A", "TSPO", "TXNIP", "VDAC2", "WEE1", "XIAP"
    ),
    negative = character()
  ),
  BUFFA_HYPOXIA_METAGENE = list(
    positive = c(
      "ACOT7", "ADM", "AK4", "ANKRD37", "ANLN", "BNIP3", "CA9", "CDKN3", "CHCHD2", "CORO1C",
      "CTSV", "DDIT4", "ENO1", "ESRP1", "GAPDH", "GPI", "HILPDA", "HK2", "KIF20A", "KIF4A",
      "LDHA", "LRRC42", "MAD2L2", "MAP7D1", "MCTS1", "MIF", "MRGBP", "MRPL13", "MRPL15",
      "MRPS17", "NDRG1", "P4HA1", "PFKP", "PGAM1", "PGK1", "PNP", "PSMA7", "PSRC1",
      "SEC61G", "SHCBP1", "SLC16A1", "SLC25A32", "SLC2A1", "TPI1", "TUBA1A", "TUBA1C",
      "TUBB6", "UTP11", "VEGFA", "YKT6"
    ),
    negative = character()
  ),
  WINTER_HYPOXIA_UP = list(
    positive = c(
      "ADORA2B", "AK4", "ALDOA", "ANGPTL4", "ANKRD37", "ANKRD9", "ANLN", "B4GALT2", "BCAR1",
      "BMS1", "BNIP3", "C16orf74", "CA12", "CA9", "CDCA4", "CNIH4", "COL4A5", "CORO1C",
      "CXCL8", "DPM2", "EEF1AKMT4", "EIF2S1", "GAPDH", "GEMIN2", "GMFB", "GPN3", "GSS",
      "HAUS2", "HES2", "HILPDA", "HOMER1", "IGF2BP2", "KCTD11", "KRT17", "LDHA", "LDLR",
      "METTL22", "MIF", "MNAT1", "MRGBP", "MRPL14", "MRPS17", "MTX1", "NDRG1", "NDUFA4L2",
      "NME1", "NTMT1", "NUDT15", "P4HA1", "PAWR", "PDZD11", "PEDS1", "PFKFB4", "PGAM1",
      "PGF", "PGK1", "PLAU", "PLEKHG3", "PPARD", "PPP4R1", "PSMA7", "PSMB7", "PSMD2",
      "PTGFRN", "PVR", "PYGL", "RAN", "RHOC", "RNF24", "RNPS1", "RUVBL2", "S100A10",
      "S100A3", "SLC16A1", "SLC2A1", "SLC6A8", "SLCO1B3", "SLIRP", "SNX24", "TANC2",
      "TEAD4", "TFAP2C", "TIMM23", "TMEM30B", "TMTC3", "TNS4", "TPBG", "TPD52L2", "TPI1",
      "TRMT5", "TUBB4B", "VAPB", "VEGFA", "VEZT", "XPO5"
    ),
    negative = character()
  ),
  SINGH_NFE2L2_TARGETS = list(
    positive = c(
      "ABCC1", "ABCC2", "G6PD", "GCLC", "GCLM", "GRK6", "GSR", "GSTM4", "HMOX1", "ME1",
      "NQO1", "PRDX1", "TXN", "TXNRD1"
    ),
    negative = character()
  ),
  CORDENONSI_YAP_CONSERVED_SIGNATURE = list(
    positive = c(
      "AGFG2", "AMOTL2", "ANKRD1", "ASAP1", "AXL", "BICC1", "BIRC5", "CAVIN2", "CCN1",
      "CCN2", "CDC20", "CDKN2C", "CENPF", "COL4A3", "CRIM1", "DAB2", "DDAH1", "DLC1",
      "DUSP1", "DUT", "ECT2", "EMP2", "ETV5", "FGF2", "FLNA", "FSCN1", "FSTL1", "GADD45B",
      "GAS6", "GGH", "GLS", "HEXB", "HMMR", "ITGB2", "ITGB5", "LHFPL6", "MARCKS", "MDFIC",
      "NDRG1", "PDLIM2", "PHGDH", "PMP22", "SCHIP1", "SERPINE1", "SGK1", "SH2D4A", "SHCBP1",
      "SLIT2", "STMN1", "TGFB2", "TGM2", "THBS1", "TK1", "TNNT2", "TNS1", "TOP2A", "TSPAN3"
    ),
    negative = character()
  ),
  HUANG_IFN_FIVE_GENE = list(
    positive = c(
      "IFIT1", "IFIT2", "IFI27", "SIGLEC1", "IFI44L"
    ),
    negative = character()
  ),
  ROONEY_NRF2_BIOMARKER = list(
    positive = c(
      "AKR1B10", "AKR1C1", "MAP2", "AKR1C3", "SPANXA1", "TXNRD1", "SLC7A11", "C1S",
      "AGPAT9", "KIAA1324", "PTGR1", "NQO1", "ABCC2", "SRXN1", "GLA", "SLC6A6", "SPSB1",
      "GCLC", "GCLM", "SLC3A2", "CABYR", "ABHD2", "FDXR", "FADS3", "ADM", "B4GALNT1",
      "RRAGD", "AIFM2", "ABCB6", "TSKU", "ALAS1", "PQLC2", "GSR", "ARHGAP5-AS1", "ASF1A",
      "FECH", "EDC3", "STK16", "TMEM70", "POP4", "CBR1", "ATP6V0B", "SH3BGR", "ZNF598",
      "PELO", "LARP6", "SLC48A1", "TXN", "IL17RC", "TALDO1", "MPP3", "BLOC1S4", "DFFA",
      "LOC729680", "BCAN", "C2CD2", "SLC25A51", "STAMBP", "ZNF746", "CTNS", "CLN3",
      "UBQLN1", "GOLT1B", "POLR3C", "SSU72", "CAPN10", "PRRC2B", "ZDHHC18"
    ),
    negative = c(
      "DHRS3", "PAQR8", "FAM117B", "ZAK", "DOCK4", "PNPLA4", "TMTC4", "LZTFL1", "TBC1D8B",
      "TBL1XR1", "PIK3C3", "ALCAM", "USP16", "SH3PXD2A", "PDPK1", "TP63", "MYH10", "LTBP1",
      "WWOX", "NFIA", "HNF4G", "HNRNPU", "YTHDC1", "PTP4A2", "FAM208A", "RRBP1", "FRAS1",
      "PPIA", "RPS6KA3", "UACA", "TIA1", "IVNS1ABP", "ARHGAP18", "TFPI", "LOC389831",
      "FAM162A", "FAM176A", "ADH5", "LRBA", "MTPAP", "CPS1", "PNISR", "SCARB1", "ALDH7A1",
      "TRIQK", "PODXL", "OAS1", "MARCKS", "ZKSCAN1", "GLUD1", "PCM1", "PDCD4", "HNMT",
      "RPS15A", "HMGA2", "CFI", "KLF3", "ALPK2", "DHTKD1", "TACC1", "ST3GAL5", "PCYOX1",
      "COBLL1", "GPX8", "FLRT3", "NAV2", "PPFIBP2", "PTPRK", "STARD4", "PIK3AP1", "CADM1",
      "ST6GAL1", "CARHSP1", "ACSM3", "NREP"
    )
  )
)

# Metadata uses the same names as signatures; counts are source-human, not pig coverage.
signature_metadata <- list(
  HALLMARK_INTERFERON_ALPHA_RESPONSE = list(
    program = "IFN_type_I_associated",
    role = "core",
    species = "Homo sapiens",
    identifier = "gene symbol (frozen source spelling)",
    source = "MSigDB",
    source_version = "2026.1.Hs",
    source_id = "M5911",
    pmid = "26771021",
    doi = "10.1016/j.cels.2015.12.004",
    source_table = "Hallmark collection",
    source_url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/h.all.v2026.1.Hs.symbols.gmt",
    source_sha256 = "eecaf6dad908334ae885406ec72bdc0646d8917588ed7c219fac92fc5363f596",
    status = "source_verified; pig_mapping_pending",
    note = "IFN-associated response, not ligand concentration or type-I exclusivity; overlaps IFN-gamma.",
    membership_sha256 = "4974a37e8444a03247377556f38a5033ca3fdca20333cf240feb9235d8403245",
    n_positive = 97L,
    n_negative = 0L
  ),
  HALLMARK_INTERFERON_GAMMA_RESPONSE = list(
    program = "IFN_type_II_associated",
    role = "core",
    species = "Homo sapiens",
    identifier = "gene symbol (frozen source spelling)",
    source = "MSigDB",
    source_version = "2026.1.Hs",
    source_id = "M5913",
    pmid = "26771021",
    doi = "10.1016/j.cels.2015.12.004",
    source_table = "Hallmark collection",
    source_url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/h.all.v2026.1.Hs.symbols.gmt",
    source_sha256 = "eecaf6dad908334ae885406ec72bdc0646d8917588ed7c219fac92fc5363f596",
    status = "source_verified; pig_mapping_pending",
    note = "IFN-gamma-associated response; overlapping IFN modules are not independent readouts.",
    membership_sha256 = "85de0a73f3a7247c1651f80ca98c6d4d2057a3c339a6f2aba3f0a1689bc9dd0e",
    n_positive = 200L,
    n_negative = 0L
  ),
  HALLMARK_HYPOXIA = list(
    program = "hypoxia_response",
    role = "core",
    species = "Homo sapiens",
    identifier = "gene symbol (frozen source spelling)",
    source = "MSigDB",
    source_version = "2026.1.Hs",
    source_id = "M5891",
    pmid = "26771021",
    doi = "10.1016/j.cels.2015.12.004",
    source_table = "Hallmark collection",
    source_url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/h.all.v2026.1.Hs.symbols.gmt",
    source_sha256 = "eecaf6dad908334ae885406ec72bdc0646d8917588ed7c219fac92fc5363f596",
    status = "source_verified; pig_mapping_pending",
    note = "Hypoxia-associated transcription, not oxygen concentration.",
    membership_sha256 = "b0cb8022d895de3d1a4a4e9ff0063a9c1c66c30709534efb275c8f685a77c1c4",
    n_positive = 200L,
    n_negative = 0L
  ),
  HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY = list(
    program = "oxidative_stress_response",
    role = "core",
    species = "Homo sapiens",
    identifier = "gene symbol (frozen source spelling)",
    source = "MSigDB",
    source_version = "2026.1.Hs",
    source_id = "M5938",
    pmid = "26771021",
    doi = "10.1016/j.cels.2015.12.004",
    source_table = "Hallmark collection",
    source_url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/h.all.v2026.1.Hs.symbols.gmt",
    source_sha256 = "eecaf6dad908334ae885406ec72bdc0646d8917588ed7c219fac92fc5363f596",
    status = "source_verified; pig_mapping_pending",
    note = "Includes protective antioxidant responses; not a measure of ROS burden or oxidative damage.",
    membership_sha256 = "1044e0ca899ec6c785628c0d2b081fa16bda5f524f30433c1274e82c81be778b",
    n_positive = 49L,
    n_negative = 0L
  ),
  HALLMARK_UNFOLDED_PROTEIN_RESPONSE = list(
    program = "unfolded_protein_response",
    role = "core",
    species = "Homo sapiens",
    identifier = "gene symbol (frozen source spelling)",
    source = "MSigDB",
    source_version = "2026.1.Hs",
    source_id = "M5922",
    pmid = "26771021",
    doi = "10.1016/j.cels.2015.12.004",
    source_table = "Hallmark collection",
    source_url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/h.all.v2026.1.Hs.symbols.gmt",
    source_sha256 = "eecaf6dad908334ae885406ec72bdc0646d8917588ed7c219fac92fc5363f596",
    status = "source_verified; pig_mapping_pending",
    note = "Broad unfolded-protein response, not an ATF4-specific integrated-stress module or ER damage assay.",
    membership_sha256 = "dde36b9c13d69d4b0ac49a50d2a0e42572ad2de28aa743f670aa5fa736ae4e0e",
    n_positive = 113L,
    n_negative = 0L
  ),
  HALLMARK_TNFA_SIGNALING_VIA_NFKB = list(
    program = "TNF_NFKB_response",
    role = "core",
    species = "Homo sapiens",
    identifier = "gene symbol (frozen source spelling)",
    source = "MSigDB",
    source_version = "2026.1.Hs",
    source_id = "M5890",
    pmid = "26771021",
    doi = "10.1016/j.cels.2015.12.004",
    source_table = "Hallmark collection",
    source_url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/h.all.v2026.1.Hs.symbols.gmt",
    source_sha256 = "eecaf6dad908334ae885406ec72bdc0646d8917588ed7c219fac92fc5363f596",
    status = "source_verified; pig_mapping_pending",
    note = "NF-kappaB/TNF-associated response, not TNF-specific causal evidence.",
    membership_sha256 = "df1d9739b55a59056ce8c11fda205d3064d2e7d18cf4bbd3df690ecd43f0d99c",
    n_positive = 200L,
    n_negative = 0L
  ),
  HALLMARK_INFLAMMATORY_RESPONSE = list(
    program = "inflammatory_response",
    role = "core",
    species = "Homo sapiens",
    identifier = "gene symbol (frozen source spelling)",
    source = "MSigDB",
    source_version = "2026.1.Hs",
    source_id = "M5932",
    pmid = "26771021",
    doi = "10.1016/j.cels.2015.12.004",
    source_table = "Hallmark collection",
    source_url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/h.all.v2026.1.Hs.symbols.gmt",
    source_sha256 = "eecaf6dad908334ae885406ec72bdc0646d8917588ed7c219fac92fc5363f596",
    status = "source_verified; pig_mapping_pending",
    note = "Broad inflammatory program; cell identity, handling and composition can contribute.",
    membership_sha256 = "f87251214a8124d52c79d17e7ab3b5edd9c766e5af92727593ec28ee218d7c10",
    n_positive = 200L,
    n_negative = 0L
  ),
  HALLMARK_P53_PATHWAY = list(
    program = "p53_associated_program",
    role = "extension",
    species = "Homo sapiens",
    identifier = "gene symbol (frozen source spelling)",
    source = "MSigDB",
    source_version = "2026.1.Hs",
    source_id = "M5939",
    pmid = "26771021",
    doi = "10.1016/j.cels.2015.12.004",
    source_table = "Hallmark collection",
    source_url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/h.all.v2026.1.Hs.symbols.gmt",
    source_sha256 = "eecaf6dad908334ae885406ec72bdc0646d8917588ed7c219fac92fc5363f596",
    status = "source_verified; pig_mapping_pending",
    note = "Pathway membership; higher expression is not a specific measurement of DNA lesions or p53 activity.",
    membership_sha256 = "40711ce465b3fd2041d12dfbe9952ed882f4fb9af11f561c45725166ce8ba055",
    n_positive = 200L,
    n_negative = 0L
  ),
  HALLMARK_DNA_REPAIR = list(
    program = "DNA_repair_program",
    role = "extension",
    species = "Homo sapiens",
    identifier = "gene symbol (frozen source spelling)",
    source = "MSigDB",
    source_version = "2026.1.Hs",
    source_id = "M5898",
    pmid = "26771021",
    doi = "10.1016/j.cels.2015.12.004",
    source_table = "Hallmark collection",
    source_url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/h.all.v2026.1.Hs.symbols.gmt",
    source_sha256 = "eecaf6dad908334ae885406ec72bdc0646d8917588ed7c219fac92fc5363f596",
    status = "source_verified; pig_mapping_pending",
    note = "Repair and proliferation can contribute; not DNA lesion burden.",
    membership_sha256 = "1806b7e4436234e6e0b558dea3194cd6db8c48538c2ea8aca918d3c94ed02078",
    n_positive = 150L,
    n_negative = 0L
  ),
  HALLMARK_APOPTOSIS = list(
    program = "apoptosis_associated_program",
    role = "extension",
    species = "Homo sapiens",
    identifier = "gene symbol (frozen source spelling)",
    source = "MSigDB",
    source_version = "2026.1.Hs",
    source_id = "M5902",
    pmid = "26771021",
    doi = "10.1016/j.cels.2015.12.004",
    source_table = "Hallmark collection",
    source_url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/h.all.v2026.1.Hs.symbols.gmt",
    source_sha256 = "eecaf6dad908334ae885406ec72bdc0646d8917588ed7c219fac92fc5363f596",
    status = "source_verified; pig_mapping_pending",
    note = "Contains death-associated regulators; transcription does not establish apoptosis execution.",
    membership_sha256 = "9a6a41cf81346af37fcc9e81a1de9c7535f8bfafdfcdb079ee6b8924ca3ca80b",
    n_positive = 161L,
    n_negative = 0L
  ),
  BUFFA_HYPOXIA_METAGENE = list(
    program = "hypoxia_response",
    role = "alternative",
    species = "Homo sapiens",
    identifier = "gene symbol (frozen source spelling)",
    source = "MSigDB",
    source_version = "2026.1.Hs",
    source_id = "M34030",
    pmid = "20087356",
    doi = "10.1038/sj.bjc.6605450",
    source_table = "Table S5 (MSigDB source)",
    source_url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/c2.cgp.v2026.1.Hs.symbols.gmt",
    source_sha256 = "170609526c4a52d89f3ce3719b616b199b3b13d8785dc672a593fabb90e72bf2",
    status = "source_verified; pig_mapping_pending",
    note = "Human cancer-derived hypoxia metagene; frozen catalog symbols, not a reproduced prognostic model.",
    membership_sha256 = "5047d8c1c2ab591b6444df4b3921fa8d1a10b189e32fa61fa523778f80377326",
    n_positive = 50L,
    n_negative = 0L
  ),
  WINTER_HYPOXIA_UP = list(
    program = "hypoxia_response",
    role = "alternative",
    species = "Homo sapiens",
    identifier = "gene symbol (frozen source spelling)",
    source = "MSigDB",
    source_version = "2026.1.Hs",
    source_id = "M5466",
    pmid = "17409455",
    doi = "10.1158/0008-5472.CAN-06-3322",
    source_table = "Table 2S (MSigDB source)",
    source_url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/c2.cgp.v2026.1.Hs.symbols.gmt",
    source_sha256 = "170609526c4a52d89f3ce3719b616b199b3b13d8785dc672a593fabb90e72bf2",
    status = "source_verified; pig_mapping_pending",
    note = "Original publication describes99 genes; frozen MSigDB mapping contains95. Cancer-derived, not validated in pig lung.",
    membership_sha256 = "5d0e6ad2ab9159c48ca2d716109722590e492523162bf17a536a0909d439d3f2",
    n_positive = 95L,
    n_negative = 0L
  ),
  SINGH_NFE2L2_TARGETS = list(
    program = "NRF2_response",
    role = "alternative",
    species = "Homo sapiens",
    identifier = "gene symbol (frozen source spelling)",
    source = "MSigDB",
    source_version = "2026.1.Hs",
    source_id = "M2662",
    pmid = "18829555",
    doi = "10.1158/0008-5472.CAN-08-1401",
    source_table = "Table 1 (MSigDB source)",
    source_url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/c2.cgp.v2026.1.Hs.symbols.gmt",
    source_sha256 = "170609526c4a52d89f3ce3719b616b199b3b13d8785dc672a593fabb90e72bf2",
    status = "source_verified; pig_mapping_pending",
    note = "Targets decrease after NRF2 knockdown in lung cancer models; higher expression is NRF2-associated response, not damage.",
    membership_sha256 = "89c1ce13715142cbfca09e21c9ecd5dd6f35c917dbff2d3318ce30cc4e7dc72d",
    n_positive = 14L,
    n_negative = 0L
  ),
  CORDENONSI_YAP_CONSERVED_SIGNATURE = list(
    program = "YAP_TAZ_associated_program",
    role = "exploratory_alternative",
    species = "Homo sapiens",
    identifier = "gene symbol (frozen source spelling)",
    source = "MSigDB",
    source_version = "2026.1.Hs",
    source_id = "M2871",
    pmid = "22078877",
    doi = "10.1016/j.cell.2011.09.048",
    source_table = "Table S3 (MSigDB source)",
    source_url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/c6.all.v2026.1.Hs.symbols.gmt",
    source_sha256 = "d127b735a03a85fdf90c1ab91d49d3cab2c4c507454f769b85e2df5ce3bb3561",
    status = "source_verified; pig_mapping_pending",
    note = "Breast-cancer/Hippo-context signature; may reflect repair/proliferation/polarity, not a universal mechanical-stretch sensor.",
    membership_sha256 = "6c7d684a0a0b100625efffdd2e91a7328baac4a6905dee5268fb11d0e04e28c3",
    n_positive = 57L,
    n_negative = 0L
  ),
  HUANG_IFN_FIVE_GENE = list(
    program = "IFN_response",
    role = "alternative",
    species = "Homo sapiens",
    identifier = "gene symbol (publication spelling)",
    source = "Huang et al.",
    source_version = "publication_2023",
    source_id = "PMID:38268924",
    pmid = "38268924",
    doi = "10.3389/fimmu.2023.1315602",
    source_table = "Article abstract: five commonly used ISGs",
    source_url = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10806211/fullTextXML",
    source_sha256 = "70eaeca2daeb7df67e068a2f8429e6be20006a5b307d44085329792bde05c0af",
    status = "source_verified; pig_mapping_pending",
    note = "Five-gene blood IFN panel; lineage/time-dependent expression. This stores membership, not the original target/reference expression-ratio scoring formula.",
    membership_sha256 = "f7018c9fb0652d473fdf29284baea959a6233d17d4b30cf72fd7409ee9acf688",
    n_positive = 5L,
    n_negative = 0L
  ),
  ROONEY_NRF2_BIOMARKER = list(
    program = "NRF2_response",
    role = "alternative",
    species = "Homo sapiens",
    identifier = "gene symbol (publication spelling; legacy aliases retained)",
    source = "Rooney et al.",
    source_version = "publication_2020_S1",
    source_id = "PMID:32986742",
    pmid = "32986742",
    doi = "10.1371/journal.pone.0239367",
    source_table = "S1 File; worksheet NRF2 Gene Biomarker; Gene and Fold Change columns",
    source_url = "https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0239367.s001&type=supplementary",
    source_sha256 = "eaf242494e5ded37330b763bc6b3c14ff6113e435a8fd329d085c3d479a80598",
    status = "source_verified; pig_mapping_pending",
    note = "143 genes;68 positive/75 negative by published Fold Change. Legacy symbols (e.g. AGPAT9) retained. Original validation used Running Fisher chemical perturbation comparisons, not UCell or pig lungs.",
    membership_sha256 = "50a43e488b25ae922179a358181e254ced2ff1d361582cbcc9ee145bbc279ea7",
    n_positive = 68L,
    n_negative = 75L
  )
)
