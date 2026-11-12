# Population count aggregation is shared; every comparison is filtered and fitted independently.
pseudobulk_comparisons <- function() {
  data.frame(contrast = c("P4_vs_Control", "N4_vs_Control", "P4_vs_N4", "P4_vs_P10", "N4_vs_N10", "P10_vs_N10"),
    group1 = c("P4", "N4", "P4", "P4", "N4", "P10"),
    group2 = c("Control", "Control", "N4", "P10", "N10", "N10"))
}

prepare_pseudobulk <- function(object, common_genes, cluster_col, cluster, sample_panel,
    min_matched_pigs = 3L) {
  if (min_matched_pigs < 3L) stop("At least three complete pig pairs are required.")
  cells <- colnames(object)[as.character(object[[cluster_col, drop = TRUE]]) == cluster &
    object$sample %in% sample_panel$sample]
  audit <- sample_panel
  audit$n_cells <- as.integer(table(factor(object$sample[match(cells, colnames(object))], levels = audit$sample)))
  audit$passes_cell_filter <- audit$n_cells >= 25L
  audit$exclusion_reason <- ifelse(audit$n_cells == 0L, "absent_cluster_population",
    ifelse(audit$n_cells < 25L, "below_25_cells", ""))
  comparisons <- pseudobulk_comparisons()
  comparisons$matched_pigs <- vapply(seq_len(nrow(comparisons)), function(i) {
    passing <- audit[audit$passes_cell_filter, ]
    length(intersect(as.character(passing$pig[passing$condition == comparisons$group1[i]]),
      as.character(passing$pig[passing$condition == comparisons$group2[i]])))
  }, integer(1))
  population <- list(sample_audit = audit, min_matched_pigs = min_matched_pigs,
    min_cells = 25L, candidate_genes = common_genes, counts = NULL)
  if (!any(comparisons$matched_pigs >= min_matched_pigs)) return(population)
  passing <- audit$sample[audit$passes_cell_filter]
  cells <- cells[object$sample[match(cells, colnames(object))] %in% passing]
  selected <- subset(object, cells = cells, features = common_genes)
  counts <- Seurat::AggregateExpression(selected, assays = "RNA", group.by = "sample",
    return.seurat = FALSE, verbose = FALSE)$RNA
  if (anyDuplicated(colnames(counts)) || !setequal(colnames(counts), passing) ||
      !setequal(rownames(counts), common_genes)) stop("Aggregated count IDs disagree with population/sample IDs.")
  counts <- counts[common_genes, passing, drop = FALSE]
  if (any(!is.finite(counts)) || any(counts < 0 | counts != round(counts)) || any(colSums(counts) <= 0)) {
    stop("Pseudobulk libraries must have nonnegative integer raw counts and positive sizes.")
  }
  population$counts <- counts
  population
}

prepare_paired_comparison <- function(population, contrast) {
  spec <- pseudobulk_comparisons()
  spec <- spec[match(contrast, spec$contrast), ]
  if (anyNA(spec)) stop("Unknown comparison: ", contrast)
  audit <- population$sample_audit
  audit$in_comparison <- audit$condition %in% c(spec$group1, spec$group2)
  passing <- audit[audit$passes_cell_filter & audit$in_comparison, ]
  pairs <- intersect(as.character(passing$pig[passing$condition == spec$group1]),
    as.character(passing$pig[passing$condition == spec$group2]))
  audit$paired_pig <- as.character(audit$pig) %in% pairs
  audit$retained <- audit$passes_cell_filter & audit$in_comparison & audit$paired_pig
  audit$exclusion_reason <- ifelse(!audit$in_comparison, "outside_comparison",
    ifelse(audit$n_cells == 0L, "absent_cluster_population",
      ifelse(!audit$passes_cell_filter, "below_25_cells", ifelse(!audit$paired_pig, "unmatched_pig", ""))))
  samples <- droplevels(audit[audit$retained, ])
  samples$condition <- factor(as.character(samples$condition), levels = c(spec$group2, spec$group1))
  eligible <- length(pairs) >= population$min_matched_pigs
  reason <- if (eligible) "" else if (!all(c(spec$group1, spec$group2) %in% as.character(passing$condition))) {
    "missing_group_after_25_cell_filter"
  } else paste0("fewer_than_", population$min_matched_pigs, "_complete_pig_pairs")
  status <- data.frame(spec, group1_passing_samples = sum(passing$condition == spec$group1),
    group2_passing_samples = sum(passing$condition == spec$group2),
    group1_samples = sum(samples$condition == spec$group1), group2_samples = sum(samples$condition == spec$group2),
    matched_pigs = length(pairs), matched_pig_ids = paste(pairs, collapse = ";"),
    min_cells = population$min_cells, min_matched_pigs = population$min_matched_pigs,
    eligible, contrast_status = if (eligible) "pending" else "skipped", skip_reason = reason)
  pb <- list(status = if (eligible) "eligible" else "no_eligible_contrast", reason = reason,
    min_cells = population$min_cells, min_matched_pigs = population$min_matched_pigs,
    model_scope = "comparison_specific_complete_pairs", sample_audit = audit, sample_info = samples,
    contrast_status = status, results = list(), design_rank = NA_integer_, residual_df = NA_integer_,
    filter_min_group_size = length(pairs), candidate_genes = population$candidate_genes)
  if (!eligible) return(pb)
  if (anyDuplicated(paste(samples$pig, samples$condition)) || !all(table(samples$pig, samples$condition) == 1L)) {
    stop("Each retained pig must contribute exactly one specimen per compared condition.")
  }
  design <- model.matrix(~ pig + condition, data = samples,
    contrasts.arg = list(pig = "contr.treatment", condition = "contr.treatment"))
  pb$design <- design
  pb$design_rank <- qr(design)$rank
  pb$residual_df <- nrow(design) - pb$design_rank
  if (pb$design_rank != ncol(design) || pb$residual_df <= 0L) stop("Invalid complete-pair model rank/residual df.")
  pb$contrasts <- matrix(0, nrow = ncol(design), ncol = 1L,
    dimnames = list(colnames(design), contrast))
  same_pig <- as.character(samples$pig) == pairs[1]
  pb$contrasts[, contrast] <- design[same_pig & samples$condition == spec$group1, ] -
    design[same_pig & samples$condition == spec$group2, ]
  pb$counts <- population$counts[, samples$sample, drop = FALSE]
  y <- edgeR::DGEList(counts = pb$counts, samples = samples)
  keep <- edgeR::filterByExpr(y, group = samples$condition,
    min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7)
  if (!any(keep)) stop("No genes retained by comparison-specific filterByExpr: ", contrast)
  y <- y[keep, , keep.lib.sizes = FALSE]
  if (any(y$samples$lib.size <= 0)) stop("A filtered pseudobulk library has zero counts.")
  y <- edgeR::normLibSizes(y, method = "TMM")
  pb$sample_info$library_size <- y$samples$lib.size
  pb$sample_info$tmm_factor <- y$samples$norm.factors
  pb$sample_info$effective_library_size <- y$samples$lib.size * y$samples$norm.factors
  pb$keep <- keep
  pb$dge <- y
  pb$parameters <- list(min_cells = 25L, min_matched_pigs = population$min_matched_pigs,
    model = "~pig+condition", model_scope = pb$model_scope,
    filterByExpr = list(min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7),
    normalization = "TMM", legacy = FALSE, robust = FALSE, fdr = "BH within population/contrast")
  pb
}

fit_pseudobulk <- function(pb, output_dir, gene_mapping, fdr_cutoff = 0.05) {
  if (pb$status != "eligible" || ncol(pb$contrasts) != 1L) stop("One eligible complete-pair comparison is required.")
  if (any(pb$sample_info$n_cells < 25L) || nlevels(pb$sample_info$pig) < 3L ||
      !all(table(pb$sample_info$pig, pb$sample_info$condition) == 1L)) stop("Pair/cell eligibility violated before fitting.")
  if (qr(pb$design)$rank != ncol(pb$design) || nrow(pb$design) <= ncol(pb$design)) {
    stop("Paired pig + condition design must be full rank with positive residual degrees of freedom.")
  }
  pb$fit <- edgeR::glmQLFit(pb$dge, pb$design, legacy = FALSE, robust = FALSE)
  test <- edgeR::glmQLFTest(pb$fit, contrast = pb$contrasts[, 1L])
  tab <- edgeR::topTags(test, n = Inf, adjust.method = "BH", sort.by = "PValue")$table
  tab <- data.frame(gene = rownames(tab), tab, row.names = NULL)
  tab$direction <- ifelse(tab$logFC > 0, "up", ifelse(tab$logFC < 0, "down", "unchanged"))
  tab <- annotate_pig_genes(tab, gene_mapping)
  contrast <- colnames(pb$contrasts)[1L]
  pb$results <- setNames(list(tab), contrast)
  pb$status <- "completed"
  pb$contrast_status$contrast_status <- "tested"
  write.csv(pb$sample_info, file.path(output_dir, "sample_metadata.csv"), row.names = FALSE)
  write.csv(pb$sample_audit, file.path(output_dir, "sample_audit.csv"), row.names = FALSE)
  write.csv(pb$contrast_status, file.path(output_dir, "contrast_status.csv"), row.names = FALSE)
  write.csv(data.frame(gene = rownames(pb$counts), tested = pb$keep),
    file.path(output_dir, "gene_filter.csv"), row.names = FALSE)
  write.csv(tab, file.path(output_dir, paste0(contrast, "_all.csv")), row.names = FALSE)
  write.csv(tab[tab$FDR <= fdr_cutoff, ], file.path(output_dir, paste0(contrast, "_significant.csv")), row.names = FALSE)
  pb$summary <- data.frame(contrast, tested = nrow(tab),
    up = sum(tab$FDR <= fdr_cutoff & tab$logFC > 0), down = sum(tab$FDR <= fdr_cutoff & tab$logFC < 0))
  write.csv(merge(pb$contrast_status, pb$summary, by = "contrast", sort = FALSE),
    file.path(output_dir, "DEG_summary.csv"), row.names = FALSE)
  saveRDS(pb, file.path(output_dir, "pseudobulk_fit.rds"))
  pb
}
