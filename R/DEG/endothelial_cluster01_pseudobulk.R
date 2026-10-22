# Run from the project root: Rscript R/DEG/endothelial_cluster01_pseudobulk.R
# Paired endothelial cluster 1 pilot: Control, positive 4H and negative 4H.
suppressPackageStartupMessages({
  library(Seurat)
  library(edgeR)
})

input <- "data/job-44542861/subsets/Endothelial/job-44831031/subclustered.rds"
cluster_col <- "harmony_snn_res.0.3"
output_dir <- "results/job-44542861/clusters_plots/subset/Endothelial/job-44831031/DEG/cluster_01"
fdr_cutoff <- 0.05
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
info_path <- file.path(output_dir, ".INFO")
info <- c(
  "# Endothelial cluster 1 paired pseudobulk pilot",
  "Status: running",
  paste("Started:", Sys.time()),
  paste("Producer:", Sys.getenv("SLURM_JOB_ID", unset = "local Rscript")),
  "Script: R/DEG/endothelial_cluster01_pseudobulk.R",
  paste("Input (read-only):", input),
  "Parents: subclustering 44831031; cleanup 44827588; pipeline 44542861; SoupX 44409807.",
  paste("Population: cluster 1 at", cluster_col),
  "Groups: paired Control (T0H), P4 and N4 (T4H); five pigs, 15 samples.",
  "T8H grouped as T10H in memory; late samples excluded from this pilot.",
  "Counts: sum original unnormalized RNA counts by sample; all eligible cells, no downsampling.",
  "Candidate genes: intersection of all 23 original RNA count layers (approved 15902 genes).",
  "Filter: filterByExpr(group=condition); min.count=10, min.total.count=15, large.n=10, min.prop=0.7.",
  "Normalization: TMM after gene filtering, with library sizes recalculated.",
  "Model: ~ pig + condition; edgeR glmQLFit legacy=FALSE, dispersion=NULL, robust=FALSE (v4 defaults).",
  "Contrasts: P4-Control, N4-Control, P4-N4; tests include both expression directions.",
  paste("Significance: BH FDR <=", fdr_cutoff, "separately per contrast; no fold-change cutoff."),
  "Positive logFC: higher in P4, N4 and P4 respectively; negative logFC: lower.",
  "Baseline comparisons combine time/perfusion and pressure-associated changes, not isolated pressure effects.",
  "Deterministic analysis; random seed not applicable. No normalization/integration of cells or source writes.",
  paste("R:", getRversion(), "; Seurat:", packageVersion("Seurat"), "; edgeR:", packageVersion("edgeR"))
)
writeLines(info, info_path)

source("R/DEG/prepare_endothelial_cluster01.R", local = TRUE)
# edgeR v4 estimates the NB dispersion inside glmQLFit; no estimateDisp step.
fit <- glmQLFit(y, design, legacy = FALSE)

# The pig coefficients and intercept cancel in every within-pig contrast.
contrasts <- matrix(0, nrow = ncol(design), ncol = 3,
  dimnames = list(colnames(design), c("P4_vs_Control", "N4_vs_Control", "P4_vs_N4")))
contrasts["conditionP4", ] <- c(1, 0, 1)
contrasts["conditionN4", ] <- c(0, 1, -1)
results <- lapply(seq_len(ncol(contrasts)), function(i) {
  test <- glmQLFTest(fit, contrast = contrasts[, i])
  tab <- topTags(test, n = Inf, adjust.method = "BH", sort.by = "PValue")$table
  tab <- data.frame(gene = rownames(tab), tab, row.names = NULL)
  tab$direction <- ifelse(tab$logFC > 0, "up", ifelse(tab$logFC < 0, "down", "unchanged"))
  tab
})
names(results) <- colnames(contrasts)

sample_info$library_size <- y$samples$lib.size
sample_info$tmm_factor <- y$samples$norm.factors
write.csv(sample_info, file.path(output_dir, "sample_metadata.csv"), row.names = FALSE)
write.csv(data.frame(gene = rownames(counts), tested = keep),
  file.path(output_dir, "gene_filter.csv"), row.names = FALSE)
for (comparison in names(results)) {
  tab <- results[[comparison]]
  write.csv(tab, file.path(output_dir, paste0(comparison, "_all.csv")), row.names = FALSE)
  write.csv(tab[tab$FDR <= fdr_cutoff, ],
    file.path(output_dir, paste0(comparison, "_significant.csv")), row.names = FALSE)
}
summary_table <- data.frame(
  contrast = names(results),
  tested = vapply(results, nrow, integer(1)),
  up = vapply(results, function(tab) sum(tab$FDR <= fdr_cutoff & tab$logFC > 0), integer(1)),
  down = vapply(results, function(tab) sum(tab$FDR <= fdr_cutoff & tab$logFC < 0), integer(1))
)
write.csv(summary_table, file.path(output_dir, "DEG_summary.csv"), row.names = FALSE)
saveRDS(list(counts = counts, sample_info = sample_info, keep = keep,
  design = design, contrasts = contrasts, dge = y, fit = fit, results = results),
  file.path(output_dir, "pseudobulk_fit.rds"))

info[info == "Status: running"] <- "Status: analysis completed; independent saved-output validation pending"
writeLines(c(info, paste("Finished:", Sys.time()),
  paste("Cells aggregated:", ncol(pilot), "; samples:", ncol(counts)),
  paste("Genes outside approved intersection:", nrow(endothelial) - length(common_genes)),
  paste("Candidate genes:", length(common_genes), "; low-expression genes removed:", sum(!keep),
    "; genes tested:", sum(keep)),
  "Outputs: sample_metadata.csv, gene_filter.csv, three full and significant result tables, DEG_summary.csv, pseudobulk_fit.rds.",
  "Validation: pending independent reopening and checks of saved counts, metadata and results."), info_path)
print(summary_table, row.names = FALSE)
