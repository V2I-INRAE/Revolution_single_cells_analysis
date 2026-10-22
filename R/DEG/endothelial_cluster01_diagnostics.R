# Run from the project root: Rscript R/DEG/endothelial_cluster01_diagnostics.R
# Sample diagnostics only. No model fitting, DE tests or sample exclusions.
suppressPackageStartupMessages({
  library(Seurat)
  library(edgeR)
  library(plotthis)
  library(ggplot2)
  library(patchwork)
})
source("R/utils/save_plots.R")
source("R/DEG/pseudobulk_diagnostics.R")

input <- "data/job-44542861/subsets/Endothelial/job-44831031/subclustered.rds"
cluster_col <- "harmony_snn_res.0.3"
output_dir <- "results/job-44542861/clusters_plots/subset/Endothelial/job-44831031/DEG/cluster_01/diagnostics"
seed <- 1234
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
info_path <- file.path(output_dir, ".INFO")
info <- c(
  "# Endothelial cluster 1 pseudobulk sample diagnostics",
  "Status: running",
  paste("Started:", Sys.time()),
  paste("Runner allocation:", Sys.getenv("SLURM_JOB_ID", unset = "none")),
  "Script: R/DEG/endothelial_cluster01_diagnostics.R; shared preparation: R/DEG/prepare_endothelial_cluster01.R.",
  paste("Input (read-only):", input),
  "Parents: subclustering 44831031; cleanup 44827588; pipeline 44542861; SoupX 44409807.",
  paste("Population: cluster 1 at", cluster_col, "; Control/T0H, P4/T4H, N4/T4H; five pigs."),
  "T8H grouped as T10H in memory; late samples outside this pilot. No exclusions within the approved panel.",
  "Aggregate: sum original RNA counts per sample, all eligible cells, no downsampling.",
  "Candidate genes: intersection of all 23 count layers (15902 genes), before selecting pilot samples.",
  "filterByExpr(group=condition) defaults; normLibSizes(method=TMM) after filtering and recalculating library sizes.",
  "MDS: edgeR plotMDS(method=logFC, top=500, gene.selection=pairwise, prior.count=2, dim.plot=1:2).",
  "MDS uses normalized log2 CPM; no animal-effect removal, ComBat, scaling, PCA or Harmony input.",
  "Libraries: counts over the approved genes before filtering; filtered library_size; effective_library_size=library_size*tmm_factor.",
  "Metrics: cells, filtered library size, TMM factor and effective library size. TMM reference 1 is not an exclusion threshold.",
  paste("Plotting: plotthis ScatterPlot and BarPlot; PNG 300 dpi and vector PDF; label-placement seed:", seed),
  "No DE tests, significance calls, volcano plots, heatmaps, pathway enrichment or source-object writes.",
  paste("R:", getRversion(), "; Seurat:", packageVersion("Seurat"), "; edgeR:", packageVersion("edgeR"),
    "; plotthis:", packageVersion("plotthis"), "; ggplot2:", packageVersion("ggplot2"),
    "; ggrepel:", packageVersion("ggrepel"), "; patchwork:", packageVersion("patchwork"))
)
writeLines(info, info_path)
source("R/DEG/prepare_endothelial_cluster01.R", local = TRUE)

sample_info$library_before_filter <- as.numeric(colSums(counts))
sample_info$library_size <- y$samples$lib.size
sample_info$tmm_factor <- y$samples$norm.factors
sample_info$effective_library_size <- sample_info$library_size * sample_info$tmm_factor
mds <- plotMDS(y, top = 500, gene.selection = "pairwise", method = "logFC",
  prior.count = 2, dim.plot = c(1, 2), plot = FALSE)
mds_data <- data.frame(sample_info, MDS1 = mds$x, MDS2 = mds$y)
write.csv(sample_info, file.path(output_dir, "sample_qc.csv"), row.names = FALSE)
write.csv(mds_data, file.path(output_dir, "mds_coordinates.csv"), row.names = FALSE)
write.csv(data.frame(gene = rownames(counts), retained = keep),
  file.path(output_dir, "gene_filter.csv"), row.names = FALSE)
saveRDS(list(counts = counts, sample_info = sample_info, keep = keep, dge = y, mds = mds),
  file.path(output_dir, "pseudobulk_qc.rds"))

render_pseudobulk_diagnostics(mds_data, mds, nrow(y), output_dir,
  "Endothelial cluster 1", seed)
info[info == "Status: running"] <- "Status: diagnostics generated; saved-output and visual validation pending"
writeLines(c(info, paste("Finished:", Sys.time()),
  paste("Aggregated:", ncol(pilot), "cells;", ncol(y), "samples;", length(common_genes),
    "candidate genes;", sum(!keep), "low-expression genes removed;", nrow(y), "genes retained."),
  "Outputs: sample_qc.csv, mds_coordinates.csv, gene_filter.csv, pseudobulk_qc.rds, sample_mds/sample_qc PNG and PDF.",
  "Validation: pending independent saved-object/count checks and rendered-figure inspection."), info_path)
print(sample_info, row.names = FALSE)
