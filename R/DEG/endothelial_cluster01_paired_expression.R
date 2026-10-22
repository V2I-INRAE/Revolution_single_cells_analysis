# Run from the project root: Rscript R/DEG/endothelial_cluster01_paired_expression.R
# Descriptive paired-expression check of genes already labelled in the volcanoes.
suppressPackageStartupMessages({
  library(edgeR)
  library(plotthis)
  library(ggplot2)
})
source("R/utils/save_plots.R")

input_dir <- "results/job-44542861/clusters_plots/subset/Endothelial/job-44831031/DEG/cluster_01"
output_dir <- file.path(input_dir, "paired_expression")
seed <- 1234
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
info_path <- file.path(output_dir, ".INFO")
info <- c(
  "# Endothelial cluster 1 paired-expression check",
  "Status: running",
  paste("Started:", Sys.time()),
  paste("Producer:", Sys.getenv("SLURM_JOB_ID", unset = "local Rscript")),
  "Script: R/DEG/endothelial_cluster01_paired_expression.R",
  paste("Read-only inputs:", file.path(input_dir, "pseudobulk_fit.rds"),
    "and plots/volcano_data.csv from the FDR<=0.01, abs(logFC)>=1 plot run."),
  "Parents: paired DE 2026-10-21 allocation 44922357; subclustering 44831031; cleanup 44827588; pipeline 44542861.",
  "Panel: cluster 1 at harmony_snn_res.0.3; Control/T0H, P4/T4H, N4/T4H; five paired pigs, all 15 samples.",
  "Selection: reuse label=TRUE from each volcano, five up and five down genes; no new gene-selection threshold.",
  "Expression: edgeR cpm(saved dge, log=TRUE, prior.count=2), saved TMM factors; no z-scores or pig-effect removal.",
  "Pairs: same pig across Control/P4/N4; lines connect conditions, not a time course; free y-axis per gene.",
  "Changes: within-pig log2 CPM difference for P4-Control, N4-Control or P4-N4; direction compared with the existing model logFC.",
  "Within-pig descriptive differences are not refitted model logFC estimates or significance tests.",
  "Post-hoc display of DE-selected genes, not independent validation or a genome-wide stability assessment; no exclusions or new DE/enrichment.",
  paste("Plotting: plotthis LinePlot, existing pig palette; PNG 300 dpi and vector PDF; seed:", seed),
  paste("R:", getRversion(), "; edgeR:", packageVersion("edgeR"),
    "; plotthis:", packageVersion("plotthis"), "; ggplot2:", packageVersion("ggplot2"))
)
writeLines(info, info_path)

saved <- readRDS(file.path(input_dir, "pseudobulk_fit.rds"))
volcano <- read.csv(file.path(input_dir, "plots", "volcano_data.csv"))
selected <- volcano[volcano$label, ]
sample_info <- saved$sample_info
sample_info <- sample_info[order(sample_info$condition, sample_info$pig), ]
pig_colors <- setNames(viz_palette[seq_along(levels(sample_info$pig))], levels(sample_info$pig))
log2_cpm <- cpm(saved$dge, log = TRUE, prior.count = 2)
expression_tables <- change_tables <- list()

for (comparison in names(saved$results)) {
  genes <- selected[selected$contrast == comparison, ]
  # Alternate directions so the left/right columns show up/down model effects.
  genes <- genes[match(as.vector(rbind(genes$gene[genes$class == "Up"],
    genes$gene[genes$class == "Down"])), genes$gene), ]
  pair <- strsplit(comparison, "_vs_", fixed = TRUE)[[1]]
  first <- sample_info[sample_info$condition == pair[1], ]
  second <- sample_info[sample_info$condition == pair[2], ]
  second <- second[match(first$pig, second$pig), ]
  changes <- do.call(rbind, lapply(seq_len(nrow(genes)), function(i) {
    gene <- genes$gene[i]
    difference <- log2_cpm[gene, first$sample] - log2_cpm[gene, second$sample]
    data.frame(contrast = comparison, gene = gene, pig = first$pig,
      first_sample = first$sample, second_sample = second$sample,
      log2_cpm_difference = unname(difference), model_logFC = genes$logFC[i],
      model_FDR = genes$FDR[i], same_direction = sign(difference) == sign(genes$logFC[i]))
  }))
  agreements <- vapply(genes$gene, function(gene) sum(changes$same_direction[changes$gene == gene]), integer(1))
  panel_labels <- sprintf("%s\nModel log2FC %.2f; FDR %.2g\nSame direction: %d/5 pigs",
    genes$gene, genes$logFC, genes$FDR, agreements)
  expression <- do.call(rbind, lapply(seq_len(nrow(genes)), function(i) {
    data.frame(contrast = comparison, gene = genes$gene[i], sample_info,
      log2_cpm = unname(log2_cpm[genes$gene[i], sample_info$sample]),
      panel = factor(panel_labels[i], levels = panel_labels))
  }))
  p <- LinePlot(expression, x = "condition", y = "log2_cpm", group_by = "pig",
    facet_by = "panel", facet_ncol = 2, facet_scales = "free_y",
    palcolor = pig_colors, pt_size = 2, line_width = 0.5, line_alpha = 0.85,
    theme = "theme_classic", theme_args = list(base_size = 11), aspect.ratio = 0.55,
    legend.position = "bottom", legend.direction = "horizontal", combine = FALSE, seed = seed,
    title = paste("Endothelial cluster 1:", gsub("_", " ", comparison)),
    subtitle = "Volcano-labelled genes: higher in the first group (left), lower (right). Five paired pigs.",
    xlab = "Condition", ylab = "TMM-normalized log2 CPM (prior count = 2)")[[1]] +
    labs(caption = "Each colour connects samples from one pig, not a time course. Y-axis varies by gene.\nDirection counts use within-pig log2 CPM differences for the named contrast.\nDE-selected genes: descriptive check, not independent validation.") +
    theme(strip.text = element_text(size = 10), plot.caption = element_text(hjust = 0),
      legend.title = element_blank())
  stem <- paste0("paired_", comparison)
  save_figure(p, output_dir, stem, width = 10, height = 12)
  ggsave(file.path(output_dir, paste0(stem, ".pdf")), p, width = 10, height = 12, bg = "white")
  expression_tables[[comparison]] <- expression[, setdiff(names(expression), "panel")]
  change_tables[[comparison]] <- changes
}
expression_data <- do.call(rbind, expression_tables)
change_data <- do.call(rbind, change_tables)
write.csv(expression_data, file.path(output_dir, "paired_expression.csv"), row.names = FALSE)
write.csv(change_data, file.path(output_dir, "paired_changes.csv"), row.names = FALSE)
info[info == "Status: running"] <- "Status: generated; independent saved-output and visual validation pending"
writeLines(c(info, paste("Finished:", Sys.time()),
  paste("Selected:", nrow(selected), "gene-contrast combinations;", length(unique(selected$gene)), "unique genes."),
  "Outputs: paired_<contrast>.png/pdf for each of three contrasts; paired_expression.csv; paired_changes.csv.",
  "Validation: pending independent saved-value, pairing, rendering and input-preservation checks."), info_path)
print(aggregate(same_direction ~ contrast + gene, change_data, sum), row.names = FALSE)
