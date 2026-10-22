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

condition_colors <- setNames(viz_palette[1:3], levels(sample_info$condition))
pig_colors <- setNames(viz_palette[1:5], levels(sample_info$pig))
axis_labels <- sprintf("Leading logFC dimension %d (%.1f%%)", 1:2, 100 * mds$var.explained[1:2])
mds_plots <- lapply(c("condition", "pig"), function(group) {
  colors <- if (group == "condition") condition_colors else pig_colors
  ScatterPlot(mds_data, x = "MDS1", y = "MDS2", color_by = group,
    color_name = if (group == "condition") "Condition" else "Pig", palcolor = colors,
    size_by = 3, alpha = 1, border_size = 0.3, legend.position = "bottom",
    theme = "theme_classic", theme_args = list(base_size = 12), combine = FALSE,
    title = if (group == "condition") "Samples by condition" else "Samples by pig",
    xlab = axis_labels[1], ylab = axis_labels[2], seed = seed)[[1]] +
    ggrepel::geom_text_repel(data = mds_data, aes(MDS1, MDS2, label = sample),
      inherit.aes = FALSE, size = 3, seed = seed, max.overlaps = Inf,
      min.segment.length = 0, box.padding = 0.4) +
    scale_x_continuous(expand = expansion(mult = 0.15)) +
    scale_y_continuous(expand = expansion(mult = 0.15)) + coord_equal()
})
mds_figure <- wrap_plots(mds_plots, ncol = 2) + plot_annotation(
  title = "Endothelial cluster 1: pseudobulk sample MDS",
  subtitle = "15 samples from five pigs. Identical coordinates in both panels; animal effects retained.",
  caption = "TMM-normalized expression; leading logFC distances using the top 500 genes per sample pair.",
  tag_levels = "A")

qc_data <- sample_info
qc_data$sample <- factor(qc_data$sample, levels = with(sample_info, sample[order(pig, condition)]))
qc_data$library_millions <- qc_data$library_size / 1e6
qc_data$effective_library_millions <- qc_data$effective_library_size / 1e6
metrics <- c(n_cells = "Cluster 1 cells", library_millions = "Library size (million counts)",
  tmm_factor = "TMM normalization factor", effective_library_millions = "Effective library size (millions)")
metric_plots <- lapply(names(metrics), function(metric) {
  p <- BarPlot(qc_data, x = "sample", y = metric, fill_by = "condition", fill_name = "Condition",
    palcolor = condition_colors, x_text_angle = 60, aspect.ratio = 0.55,
    y_min = 0, y_max = 1.12 * max(qc_data[[metric]]),
    theme = "theme_classic", theme_args = list(base_size = 12), combine = FALSE,
    xlab = "Sample (grouped by pig)", ylab = metrics[[metric]], title = metrics[[metric]], seed = seed)[[1]]
  if (metric == "tmm_factor") p <- p + geom_hline(yintercept = 1, linetype = "dashed", color = "grey40")
  p
})
qc_figure <- wrap_plots(metric_plots, ncol = 2, guides = "collect") +
  plot_annotation(title = "Endothelial cluster 1: sample-level count and normalization diagnostics",
    subtitle = paste("All 15 samples retained. Libraries use the", format(nrow(y), big.mark = ","),
      "genes retained by filterByExpr."),
    caption = "Effective library size = filtered library size x TMM factor. Dashed line marks TMM factor 1, not an exclusion cutoff.",
    tag_levels = "A") & theme(legend.position = "bottom")

save_figure(mds_figure, output_dir, "sample_mds", width = 14, height = 6.5)
ggsave(file.path(output_dir, "sample_mds.pdf"), mds_figure, width = 14, height = 6.5)
save_figure(qc_figure, output_dir, "sample_qc", width = 14, height = 10)
ggsave(file.path(output_dir, "sample_qc.pdf"), qc_figure, width = 14, height = 10)
info[info == "Status: running"] <- "Status: diagnostics generated; saved-output and visual validation pending"
writeLines(c(info, paste("Finished:", Sys.time()),
  paste("Aggregated:", ncol(pilot), "cells;", ncol(y), "samples;", length(common_genes),
    "candidate genes;", sum(!keep), "low-expression genes removed;", nrow(y), "genes retained."),
  "Outputs: sample_qc.csv, mds_coordinates.csv, gene_filter.csv, pseudobulk_qc.rds, sample_mds/sample_qc PNG and PDF.",
  "Validation: pending independent saved-object/count checks and rendered-figure inspection."), info_path)
print(sample_info, row.names = FALSE)
