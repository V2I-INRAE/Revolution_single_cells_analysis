# The driver validates the original observed specimen panel before cluster selection.
# Zero cell counts below are audit metadata only: absent populations never become libraries.
prepare_pseudobulk <- function(object, common_genes, cluster_col, cluster, sample_panel) {
  cells <- colnames(object)[as.character(object[[cluster_col, drop = TRUE]]) == cluster]
  sample_audit <- sample_panel
  sample_audit$n_cells <- as.integer(table(factor(object$sample[match(cells, colnames(object))],
    levels = sample_audit$sample)))
  sample_audit$retained <- sample_audit$n_cells >= 10L
  sample_audit$exclusion_reason <- ifelse(sample_audit$n_cells == 0L, "absent_cluster_population",
    ifelse(sample_audit$n_cells < 10L, "below_10_cells", ""))
  sample_info <- droplevels(sample_audit[sample_audit$retained, ])
  group_counts <- do.call(rbind, lapply(levels(sample_panel$condition), function(group) {
    x <- sample_audit[sample_audit$condition == group, ]
    data.frame(condition = group, panel_samples = nrow(x), observed_samples = sum(x$n_cells > 0L),
      retained_samples = sum(x$retained), retained_pigs = length(unique(x$pig[x$retained])),
      excluded_samples = sum(!x$retained), absent_samples = sum(x$n_cells == 0L),
      below_10_samples = sum(x$n_cells > 0L & x$n_cells < 10L),
      observed_cells = sum(x$n_cells), retained_cells = sum(x$n_cells[x$retained]))
  }))
  contrast_status <- data.frame(
    contrast = c("P4_vs_Control", "N4_vs_Control", "P4_vs_N4", "P4_vs_P10", "N4_vs_N10", "P10_vs_N10"),
    group1 = c("P4", "N4", "P4", "P4", "N4", "P10"),
    group2 = c("Control", "Control", "N4", "P10", "N10", "N10"),
    group1_samples = 0L, group2_samples = 0L, matched_pigs = 0L, matched_pig_ids = "",
    eligible = FALSE, contrast_status = "skipped", skip_reason = "")
  for (i in seq_len(nrow(contrast_status))) {
    first <- as.character(sample_info$pig[sample_info$condition == contrast_status$group1[i]])
    second <- as.character(sample_info$pig[sample_info$condition == contrast_status$group2[i]])
    paired <- intersect(first, second)
    contrast_status$group1_samples[i] <- length(first)
    contrast_status$group2_samples[i] <- length(second)
    contrast_status$matched_pigs[i] <- length(paired)
    contrast_status$matched_pig_ids[i] <- paste(paired, collapse = ";")
    contrast_status$eligible[i] <- length(paired) >= 3L
    contrast_status$contrast_status[i] <- if (length(paired) >= 3L) "pending" else "skipped"
    contrast_status$skip_reason[i] <- if (!length(first) || !length(second)) {
      "missing_group_after_sample_filter"
    } else if (length(paired) < 3L) "fewer_than_3_directly_matched_pigs" else ""
  }
  pb <- list(status = "eligible", reason = "", sample_audit = sample_audit,
    sample_info = sample_info, group_counts = group_counts, contrast_status = contrast_status,
    contrasts = matrix(numeric(), nrow = 0L, ncol = 0L, dimnames = list(NULL, character())),
    results = list(), design_rank = NA_integer_, residual_df = NA_integer_,
    filter_min_group_size = if (nrow(sample_info)) min(as.integer(table(sample_info$condition))) else NA_integer_)
  if (!any(contrast_status$eligible)) {
    pb$status <- "no_eligible_contrasts"
    pb$reason <- "No requested contrast has at least 3 directly matched pigs after sample filtering."
    return(pb)
  }
  # Only observed levels enter the joint model; all passing specimens remain, including sparse late groups.
  design <- model.matrix(~ pig + condition, data = sample_info,
    contrasts.arg = list(pig = "contr.treatment", condition = "contr.treatment"))
  pb$design <- design
  pb$design_rank <- qr(design)$rank
  pb$residual_df <- nrow(design) - pb$design_rank
  if (pb$design_rank != ncol(design) || pb$residual_df <= 0L) {
    pb$status <- "blocked"
    pb$reason <- if (pb$design_rank != ncol(design)) "rank_deficient_joint_design" else "nonpositive_residual_df"
    pending <- pb$contrast_status$eligible
    pb$contrast_status$contrast_status[pending] <- "skipped"
    pb$contrast_status$skip_reason[pending] <- pb$reason
    return(pb)
  }
  supported <- which(contrast_status$eligible)
  pb$contrasts <- matrix(0, nrow = ncol(design), ncol = length(supported),
    dimnames = list(colnames(design), contrast_status$contrast[supported]))
  for (i in supported) {
    # Same-pig design-row differences work even when Control is absent or the reference level changes.
    first <- which(sample_info$condition == contrast_status$group1[i])
    second <- which(sample_info$condition == contrast_status$group2[i])
    paired <- intersect(as.character(sample_info$pig[first]), as.character(sample_info$pig[second]))[1]
    first <- first[as.character(sample_info$pig[first]) == paired]
    second <- second[as.character(sample_info$pig[second]) == paired]
    pb$contrasts[, contrast_status$contrast[i]] <- design[first, ] - design[second, ]
  }
  cells <- cells[object$sample[match(cells, colnames(object))] %in% sample_info$sample]
  population <- subset(object, cells = cells, features = common_genes)
  counts <- Seurat::AggregateExpression(population, assays = "RNA", group.by = "sample",
    return.seurat = FALSE, verbose = FALSE)$RNA
  if (anyDuplicated(colnames(counts)) || !setequal(colnames(counts), sample_info$sample) ||
      !setequal(rownames(counts), common_genes)) {
    stop("Aggregated count IDs do not match the retained specimens/common-layer genes.")
  }
  counts <- counts[common_genes, sample_info$sample, drop = FALSE]
  if (any(!is.finite(counts)) || any(counts < 0 | counts != round(counts)) || any(colSums(counts) <= 0)) {
    stop("Pseudobulk libraries must contain nonnegative integer raw counts and have positive sizes.")
  }
  y <- edgeR::DGEList(counts = counts, samples = sample_info)
  # Keep the approved group-based rule: a retained singleton group reduces its minimum group size to 1.
  keep <- edgeR::filterByExpr(y, group = sample_info$condition)
  if (!any(keep)) stop("No genes retained by filterByExpr; cluster cannot be fitted.")
  y <- y[keep, , keep.lib.sizes = FALSE]
  if (any(y$samples$lib.size <= 0)) stop("A filtered pseudobulk library has zero counts.")
  y <- edgeR::normLibSizes(y, method = "TMM")
  sample_info$library_size <- y$samples$lib.size
  sample_info$tmm_factor <- y$samples$norm.factors
  sample_info$effective_library_size <- sample_info$library_size * sample_info$tmm_factor
  pb$counts <- counts
  pb$sample_info <- sample_info
  pb$keep <- keep
  pb$dge <- y
  pb
}

fit_pseudobulk <- function(pb, output_dir, fdr_cutoff = 0.05) {
  if (pb$status != "eligible" || !ncol(pb$contrasts)) stop("No supported contrasts/model available to fit.")
  if (qr(pb$design)$rank != ncol(pb$design) || nrow(pb$design) <= ncol(pb$design)) {
    stop("Joint pig + condition design must be full rank with positive residual degrees of freedom.")
  }
  # edgeR v4 estimates the NB dispersion internally; retain the approved settings.
  pb$fit <- edgeR::glmQLFit(pb$dge, pb$design, legacy = FALSE, robust = FALSE)
  contrasts <- pb$contrasts
  pb$results <- setNames(lapply(seq_len(ncol(contrasts)), function(i) {
    test <- edgeR::glmQLFTest(pb$fit, contrast = contrasts[, i])
    tab <- edgeR::topTags(test, n = Inf, adjust.method = "BH", sort.by = "PValue")$table
    tab <- data.frame(gene = rownames(tab), tab, row.names = NULL)
    tab$direction <- ifelse(tab$logFC > 0, "up", ifelse(tab$logFC < 0, "down", "unchanged"))
    tab
  }), colnames(contrasts))
  pb$status <- "completed"
  pb$contrast_status$contrast_status[pb$contrast_status$contrast %in% names(pb$results)] <- "tested"
  write.csv(pb$sample_info, file.path(output_dir, "sample_metadata.csv"), row.names = FALSE)
  write.csv(pb$sample_audit, file.path(output_dir, "sample_audit.csv"), row.names = FALSE)
  write.csv(pb$contrast_status, file.path(output_dir, "contrast_status.csv"), row.names = FALSE)
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
  summary <- merge(pb$contrast_status, pb$summary, by = "contrast", all.x = TRUE, sort = FALSE)
  summary <- summary[match(pb$contrast_status$contrast, summary$contrast), , drop = FALSE]
  write.csv(summary, file.path(output_dir, "DEG_summary.csv"), row.names = FALSE)
  saveRDS(pb, file.path(output_dir, "pseudobulk_fit.rds"))
  pb
}
