# Run from the project root: Rscript R/DEG/endothelial_cluster01_plots.R
# Add --volcano-only for the FDR 0.05 / absolute log2FC 0.5 volcano rerun.
# Display the saved paired DE results; no aggregation, model fitting or DE tests.
suppressPackageStartupMessages({
  library(edgeR)
  library(ComplexHeatmap)
  library(circlize)
  library(plotthis)
  library(ggplot2)
})
source("R/utils/save_plots.R")
source("R/DEG/pseudobulk_plots.R")

input_dir <- "results/job-44542861/clusters_plots/subset/Endothelial/job-44831031/DEG/cluster_01"
volcano_only <- "--volcano-only" %in% commandArgs(trailingOnly = TRUE)
output_dir <- file.path(input_dir, "plots",
  if (volcano_only) "volcano_fdr0.05_logfc0.5" else "")
fdr_cutoff <- if (volcano_only) 0.05 else 0.01
logfc_cutoff <- if (volcano_only) 0.5 else 1
mapping_file <- "results/misc/Ensembl_ID/gene-mapping/pig_to_human_orthologues_annotated.csv"
seed <- 1234
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
info_path <- file.path(output_dir, ".INFO")
info <- c(
  if (volcano_only) "# Endothelial cluster 1 relaxed-threshold volcano plots" else
    "# Endothelial cluster 1 pseudobulk heatmap and volcano plots",
  "Status: running",
  paste("Started:", Sys.time()),
  paste("Runner allocation:", Sys.getenv("SLURM_JOB_ID", unset = "none")),
  paste("Script: R/DEG/endothelial_cluster01_plots.R", if (volcano_only) "--volcano-only" else ""),
  paste("Read-only input:", file.path(input_dir, "pseudobulk_fit.rds")),
  paste("Label mapping:", mapping_file,
    "; pig_safe_name only, deduplicated by input_gene_id. Original gene IDs retained; no human-symbol substitution."),
  "Parents: local paired DE 2026-10-21 in allocation 44922357; subclustering 44831031; cleanup 44827588; pipeline 44542861.",
  "Panel: Control/T0H, P4/T4H, N4/T4H; five pigs, 15 sample pseudobulks; cluster 1 at harmony_snn_res.0.3.",
  if (volcano_only) "Scope: volcanoes only; existing heatmap and enrichment outputs are not regenerated." else c(
    sprintf("Heatmap genes: union passing FDR<=%g AND abs(logFC)>=%g in the same contrast; no top-gene cap.",
      fdr_cutoff, logfc_cutoff),
    "Heatmap values: edgeR cpm(saved dge, log=TRUE, prior.count=2), using saved TMM factors; gene-wise z-score across all 15 samples.",
    "Heatmap order: Control, P4, N4; the same five pigs within each block; Euclidean/complete row clustering on unclipped z-scores, no dendrogram reordering or module cuts.",
    "Heatmap display: reversed RColorBrewer RdBu, blue-low/near-white-zero/red-high; color saturation +/-2.5 only; no pig-effect removal or further scaling."),
  "Volcanoes: all tested genes separately for P4-Control, N4-Control, P4-N4; exact x=logFC, y=-log10(existing BH FDR), shared axes, no trimming or jitter.",
  sprintf("Volcano classes: FDR>%g grey; FDR<=%g and abs(logFC)<%g black; FDR<=%g and logFC<=-%g blue; FDR<=%g and logFC>=%g red.",
    fdr_cutoff, fdr_cutoff, logfc_cutoff, fdr_cutoff, logfc_cutoff, fdr_cutoff, logfc_cutoff),
  "Labels: up to five genes per outer direction; order by FDR, PValue, decreasing abs(logFC), gene; display safe pig name where available, otherwise original identifier.",
  paste("Seed for label placement:", seed),
  "Outputs: PNG 300 dpi and PDF; heatmap body rasterized without cell borders, volcano points vector in PDF; source CSVs preserve exact values.",
  "Interpretation: exploratory patterns in DE-selected genes, not independently validated signatures; baseline contrasts combine time/perfusion and pressure-associated changes.",
  "No new DE tests or changes to original result tables. The logFC cutoff is an observed-effect display filter, not an effect-size null test.",
  paste("R:", getRversion(), "; edgeR:", packageVersion("edgeR"),
    "; ComplexHeatmap:", packageVersion("ComplexHeatmap"), "; circlize:", packageVersion("circlize"),
    "; plotthis:", packageVersion("plotthis"), "; ggplot2:", packageVersion("ggplot2"),
    "; ggrepel:", packageVersion("ggrepel"))
)
writeLines(info, info_path)

saved <- readRDS(file.path(input_dir, "pseudobulk_fit.rds"))
rename_map <- unique(read.csv(mapping_file)[, c("input_gene_id", "pig_safe_name")])
rd_bu <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
if (!volcano_only) {
  sample_info <- saved$sample_info
  sample_info <- sample_info[order(sample_info$condition, sample_info$pig), ]
  log2_cpm <- cpm(saved$dge, log = TRUE, prior.count = 2)
  passes <- vapply(saved$results, function(tab) {
    genes <- tab$gene[tab$FDR <= fdr_cutoff & abs(tab$logFC) >= logfc_cutoff]
    rownames(log2_cpm) %in% genes
  }, logical(nrow(log2_cpm)))
  rownames(passes) <- rownames(log2_cpm)
  selected <- rownames(log2_cpm)[rowSums(passes) > 0]
  z_scores <- t(scale(t(log2_cpm[selected, sample_info$sample, drop = FALSE])))
  row_hc <- hclust(dist(z_scores, method = "euclidean"), method = "complete")
  drawn <- render_pseudobulk_heatmap(z_scores, row_hc, sample_info, output_dir,
    "Endothelial cluster 1", fdr_cutoff, logfc_cutoff)
  gene_order <- rownames(z_scores)[row_order(drawn)]
  write.csv(data.frame(gene = gene_order, z_scores[gene_order, ], check.names = FALSE),
    file.path(output_dir, "heatmap_zscores.csv"), row.names = FALSE)
  write.csv(data.frame(gene = gene_order, row_order = seq_along(gene_order),
    passes[gene_order, ], check.names = FALSE),
    file.path(output_dir, "heatmap_genes.csv"), row.names = FALSE)
}

class_colors <- c("Not significant" = "grey75", "Significant, small effect" = "black",
  "Down" = rd_bu[1], "Up" = rd_bu[11])
volcano_data <- lapply(names(saved$results), function(comparison) {
  tab <- saved$results[[comparison]]
  safe_symbol <- rename_map$pig_safe_name[match(tab$gene, rename_map$input_gene_id)]
  tab$display_gene <- ifelse(!is.na(safe_symbol) & nzchar(safe_symbol), safe_symbol, tab$gene)
  tab$contrast <- comparison
  tab$neg_log10_fdr <- -log10(tab$FDR)
  tab$class <- "Not significant"
  tab$class[tab$FDR <= fdr_cutoff] <- "Significant, small effect"
  tab$class[tab$FDR <= fdr_cutoff & tab$logFC <= -logfc_cutoff] <- "Down"
  tab$class[tab$FDR <= fdr_cutoff & tab$logFC >= logfc_cutoff] <- "Up"
  tab$class <- factor(tab$class, levels = names(class_colors))
  tab$label <- FALSE
  ranked <- order(tab$FDR, tab$PValue, -abs(tab$logFC), tab$gene)
  for (direction in c("Down", "Up")) {
    indices <- ranked[tab$class[ranked] == direction]
    tab$label[head(indices, 5)] <- TRUE
  }
  tab
})
names(volcano_data) <- names(saved$results)
plot_data <- do.call(rbind, volcano_data)
write.csv(plot_data, file.path(output_dir, "volcano_data.csv"), row.names = FALSE)
class_counts <- as.data.frame(table(plot_data$contrast, plot_data$class),
  responseName = "n_genes")
names(class_counts)[1:2] <- c("contrast", "class")
write.csv(class_counts, file.path(output_dir, "volcano_summary.csv"), row.names = FALSE)
render_pseudobulk_volcano(plot_data, output_dir, "Endothelial cluster 1",
  fdr_cutoff, logfc_cutoff, seed)

info[info == "Status: running"] <- "Status: plots generated; one saved-output and visual review pending"
writeLines(c(info, paste("Finished:", Sys.time()),
  if (!volcano_only) c(
    paste("Heatmap:", length(gene_order), "genes x", nrow(sample_info), "samples;",
      sum(rowSums(passes) == 0), "tested genes outside the joint-threshold union."),
    "Outputs: heatmap_groups.png/pdf; heatmap_zscores.csv, heatmap_genes.csv."),
  "Outputs: volcano_<contrast>.png/pdf for each of three contrasts; volcano_data.csv, volcano_summary.csv.",
  "Validation: no extended audits or repeated validation requested; saved outputs and visual review reported separately."), info_path)
print(class_counts, row.names = FALSE)
