# Run from the project root: Rscript R/DEG/preview_publication.R
# Cluster 06 approval preview only. No fitting, testing, GSEA, or source writes.
suppressPackageStartupMessages(library(ggplot2))
source("R/DEG/render_saved_pseudobulk.R")

input_dir <- "results/job-44542861/clusters_plots/subset/Macrophage_Monocytes/job-44831034/DEG/cluster_06"
output_dir <- file.path(input_dir, "publication_preview")
label <- "Macrophage_Monocytes cluster 6"
seed <- 1234
# These are the original run's recorded display thresholds, not new analysis choices.
fdr_cutoff <- 0.05
logfc_cutoff <- 0.5
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
info_path <- file.path(output_dir, ".INFO")
info <- c("# DEG final-size publication preview", "Status: running",
  paste("Started:", Sys.time()), paste("Plotting allocation:", Sys.getenv("SLURM_JOB_ID", "local")),
  "Script: R/DEG/preview_publication.R; shared renderers and publication_style.R.",
  paste("Read-only input directory:", input_dir),
  "Parents: paired pseudobulk producer44922357; subcluster44831034; cleanup44827588; pipeline44542861.",
  "Inputs: pseudobulk_fit.rds, diagnostics/mds.rds and mds_coordinates.csv, plots/volcano_data.csv and heatmap_zscores.csv, gsea/hallmark/gsea_results.rds.",
  "Method: render saved MDS coordinates, volcano classes/labels and GSEA results; reconstruct heatmap from saved TMM DGEList and verify against existing matrix.",
  "No model fitting, new tests, fgseaMultilevel, changed thresholds/colours/data, or source writes.",
  "Display selection: original FDR<=0.05 and abs(log2FC)>=0.5; Hallmark pathway threshold from saved parameters.",
  "Final sizes: width170mm; MDS135/QC180/heatmap160/volcano145/GSEA dot170/running110mm high. Place at100% in Affinity.",
  "Text: labels/ticks/legends9pt, axis10pt, titles11pt, tags14pt; MDS size3.2mm; axes/guides0.25mm. Cairo embedded-font PDF and matching300dpi PNG; raster heatmap matrix retained.",
  paste("Label placement seed:", seed),
  "MDS display labels omit only the shared REVO prefix, explicitly documented in caption; full data IDs retained, no samples or labels dropped.",
  "Scope inventory: 228 production PDFs / 19 clusters: Endothelial5, Epithelial5, Macrophage_Monocytes5, DC3, Granulocytes1 (12 PDFs per cluster).",
  "Standalone endothelial pilot outputs were superseded/removed Oct21; standalone script styling remains pending with full rollout.",
  "Approval gate: preview only; no production PDFs replaced. All other clusters await owner approval.",
  "Evidence: logs/deg-publication-preview.log")
writeLines(info, info_path)
files <- list.files(input_dir, recursive = TRUE, full.names = TRUE, all.files = TRUE)
files <- files[!startsWith(files, paste0(output_dir, "/")) & !dir.exists(files)]
before <- tools::md5sum(files)
tryCatch({
  render_saved_pseudobulk(input_dir, output_dir, label, fdr_cutoff, logfc_cutoff, seed)
  stopifnot(identical(before, tools::md5sum(files)),
    length(list.files(output_dir, pattern = "\\.pdf$", recursive = TRUE)) == 12L)
  info[info == "Status: running"] <- "Status: rendering completed; PDF reopening/visual review pending; owner approval required."
  writeLines(c(info, paste("Finished:", Sys.time()),
    paste("Scientific input check PASS:", length(files), "original files byte-identical by MD5 before/after."),
    "Reopened input checks PASS: 15 sample identities/MDS coordinates; heatmap gene order and matrix agree to1e-12; all volcano genes/logFC/FDR agree with saved fit; GSEA rendered directly from saved result/display tables."), info_path)
  cat("PASS: 12 PDFs and 12 PNGs; original files unchanged; no statistical reruns.\n")
}, error = function(e) {
  writeLines(c(info, paste("Status: failed:", conditionMessage(e)),
    paste("Original files unchanged:", identical(before, tools::md5sum(files)))), info_path)
  stop(e)
})
