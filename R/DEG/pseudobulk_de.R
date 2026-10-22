# Sum original RNA counts once per cluster, then reuse the filtered TMM object.
prepare_pseudobulk <- function(object, common_genes, cluster_col, cluster) {
  cells <- colnames(object)[object[[cluster_col, drop = TRUE]] == cluster &
    object$time_point %in% c("T0H", "T4H")]
  population <- subset(object, cells = cells, features = common_genes)
  counts <- Seurat::AggregateExpression(population, assays = "RNA", group.by = "sample",
    return.seurat = FALSE, verbose = FALSE)$RNA
  sample_info <- unique(population[[]][, c("sample", "pig", "pressure", "time_point")])
  sample_info <- sample_info[match(colnames(counts), sample_info$sample), ]
  rownames(sample_info) <- sample_info$sample
  sample_info$pig <- factor(sample_info$pig)
  sample_info$condition <- factor(ifelse(sample_info$time_point == "T0H", "Control",
    ifelse(sample_info$pressure == "positive", "P4", "N4")),
    levels = c("Control", "P4", "N4"))
  sample_info$n_cells <- as.integer(table(population$sample)[sample_info$sample])
  design <- model.matrix(~ pig + condition, data = sample_info)
  y <- edgeR::DGEList(counts = counts, samples = sample_info)
  keep <- edgeR::filterByExpr(y, group = sample_info$condition)
  y <- edgeR::normLibSizes(y[keep, , keep.lib.sizes = FALSE], method = "TMM")
  sample_info$library_size <- y$samples$lib.size
  sample_info$tmm_factor <- y$samples$norm.factors
  sample_info$effective_library_size <- sample_info$library_size * sample_info$tmm_factor
  list(counts = counts, sample_info = sample_info, keep = keep, design = design, dge = y)
}

fit_pseudobulk <- function(pb, output_dir, fdr_cutoff = 0.05) {
  # edgeR v4 estimates the NB dispersion internally; retain the pilot settings.
  pb$fit <- edgeR::glmQLFit(pb$dge, pb$design, legacy = FALSE)
  contrasts <- matrix(0, nrow = ncol(pb$design), ncol = 3,
    dimnames = list(colnames(pb$design), c("P4_vs_Control", "N4_vs_Control", "P4_vs_N4")))
  contrasts["conditionP4", ] <- c(1, 0, 1)
  contrasts["conditionN4", ] <- c(0, 1, -1)
  pb$contrasts <- contrasts
  pb$results <- setNames(lapply(seq_len(ncol(contrasts)), function(i) {
    test <- edgeR::glmQLFTest(pb$fit, contrast = contrasts[, i])
    tab <- edgeR::topTags(test, n = Inf, adjust.method = "BH", sort.by = "PValue")$table
    tab <- data.frame(gene = rownames(tab), tab, row.names = NULL)
    tab$direction <- ifelse(tab$logFC > 0, "up", ifelse(tab$logFC < 0, "down", "unchanged"))
    tab
  }), colnames(contrasts))
  write.csv(pb$sample_info, file.path(output_dir, "sample_metadata.csv"), row.names = FALSE)
  write.csv(data.frame(gene = rownames(pb$counts), tested = pb$keep),
    file.path(output_dir, "gene_filter.csv"), row.names = FALSE)
  for (contrast in names(pb$results)) {
    tab <- pb$results[[contrast]]
    write.csv(tab, file.path(output_dir, paste0(contrast, "_all.csv")), row.names = FALSE)
    write.csv(tab[tab$FDR <= fdr_cutoff, ],
      file.path(output_dir, paste0(contrast, "_significant.csv")), row.names = FALSE)
  }
  pb$summary <- data.frame(contrast = names(pb$results),
    tested = vapply(pb$results, nrow, integer(1)),
    up = vapply(pb$results, function(x) sum(x$FDR <= fdr_cutoff & x$logFC > 0), integer(1)),
    down = vapply(pb$results, function(x) sum(x$FDR <= fdr_cutoff & x$logFC < 0), integer(1)))
  write.csv(pb$summary, file.path(output_dir, "DEG_summary.csv"), row.names = FALSE)
  saveRDS(pb, file.path(output_dir, "pseudobulk_fit.rds"))
  pb
}
