test_use <- "wilcox"
only_pos <- TRUE
min_pct <- 0.25
logfc_threshold <- 0.25
adjusted_p_cutoff <- 0.05

markers <- list(
  Macrophages = c(
    "C1QA",
    "C1QB",
    "CD68",
    "CD163",
    "TREM2",
    "ARG1",
    "MRC1",
    "MS4A7",
    "ACP5",
    "MARCO",
    "CLECL12A",
    "LGALS3"
  ),
  cDC = c(
    "FLT3",
    "FCER1",
    "XCR1",
    "BATF3",
    "CADM1",
    "FCER1A"
  ),
  T_cells = c(
    "CD2",
    "CD3D",
    "CD3E",
    "CD5",
    "CD6",
    "CD8A",
    "CD4",
    "IL7R"
  ),
  B_cells = c(
    "CD79A", "CD79B", "MS4A1", "CD19", "PAX5"
  ),
  ASC = c("JCHAIN", "MZB1", "PRDM1", "IRF4", "XBP1"),
  Neutrophils = c("G0S2", "CXCL8", "SELL", "CSF3R"),
  AT2 = c("SFTPC", "SFTPB", "ABCA3"),
  Epithelial = c("KRT5","EpCAM")

)
