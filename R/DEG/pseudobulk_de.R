# One pseudobulk population per cluster; every contrast is filtered and fitted on its own complete pig pairs.
pseudobulk_comparisons <- function() {
  data.frame(
    contrast = c(
      "P4_vs_Control",
      "N4_vs_Control",
      "P4_vs_N4",
      "P4_vs_P10",
      "N4_vs_N10",
      "P10_vs_N10"
    ),
    group1 = c("P4", "N4", "P4", "P4", "N4", "P10"),
    group2 = c("Control", "Control", "N4", "P10", "N10", "N10")
  )
}

prepare_pseudobulk <- function(
  object,
  common_genes,
  cluster_col,
  cluster,
  sample_panel
) {
  min_cells <- deg_config$min_cells
  in_cluster <- as.character(object[[cluster_col, drop = TRUE]]) == cluster
  audit <- sample_panel
  audit$n_cells <- as.integer(table(factor(
    object$sample[in_cluster],
    levels = audit$sample
  )))
  audit$passes_cell_filter <- audit$n_cells >= min_cells
  audit$exclusion_reason <- ifelse(
    audit$n_cells == 0L,
    "absent_cluster_population",
    ifelse(!audit$passes_cell_filter, paste0("below_", min_cells, "_cells"), "")
  )
  population <- list(sample_audit = audit, counts = NULL)
  passing <- audit[audit$passes_cell_filter, ]
  comparisons <- pseudobulk_comparisons()
  matched_pigs <- vapply(
    seq_len(nrow(comparisons)),
    function(i) {
      length(intersect(
        as.character(passing$pig[passing$condition == comparisons$group1[i]]),
        as.character(passing$pig[passing$condition == comparisons$group2[i]])
      ))
    },
    integer(1)
  )
  if (!any(matched_pigs >= deg_config$min_matched_pigs)) {
    return(population)
  }
  cells <- colnames(object)[in_cluster & object$sample %in% passing$sample]
  selected <- subset(object, cells = cells, features = common_genes)
  counts <- Seurat::AggregateExpression(
    selected,
    assays = "RNA",
    group.by = "sample",
    return.seurat = FALSE,
    verbose = FALSE
  )$RNA
  population$counts <- counts[common_genes, passing$sample, drop = FALSE]
  population
}

prepare_paired_comparison <- function(population, contrast) {
  min_cells <- deg_config$min_cells
  min_pairs <- deg_config$min_matched_pigs
  spec <- pseudobulk_comparisons()
  spec <- spec[spec$contrast == contrast, ]
  audit <- population$sample_audit
  audit$in_comparison <- audit$condition %in% c(spec$group1, spec$group2)
  passing <- audit[audit$passes_cell_filter & audit$in_comparison, ]
  pairs <- intersect(
    as.character(passing$pig[passing$condition == spec$group1]),
    as.character(passing$pig[passing$condition == spec$group2])
  )
  audit$paired_pig <- as.character(audit$pig) %in% pairs
  audit$retained <- audit$passes_cell_filter &
    audit$in_comparison &
    audit$paired_pig
  audit$exclusion_reason <- ifelse(
    !audit$in_comparison,
    "outside_comparison",
    ifelse(
      audit$n_cells == 0L,
      "absent_cluster_population",
      ifelse(
        !audit$passes_cell_filter,
        paste0("below_", min_cells, "_cells"),
        ifelse(!audit$paired_pig, "unmatched_pig", "")
      )
    )
  )
  samples <- droplevels(audit[audit$retained, ])
  samples$condition <- factor(
    as.character(samples$condition),
    levels = c(spec$group2, spec$group1)
  )
  eligible <- length(pairs) >= min_pairs
  reason <- if (eligible) {
    ""
  } else if (
    !all(c(spec$group1, spec$group2) %in% as.character(passing$condition))
  ) {
    paste0("missing_group_after_", min_cells, "_cell_filter")
  } else {
    paste0("fewer_than_", min_pairs, "_complete_pig_pairs")
  }
  status <- data.frame(
    spec,
    group1_passing_samples = sum(passing$condition == spec$group1),
    group2_passing_samples = sum(passing$condition == spec$group2),
    group1_samples = sum(samples$condition == spec$group1),
    group2_samples = sum(samples$condition == spec$group2),
    matched_pigs = length(pairs),
    matched_pig_ids = paste(pairs, collapse = ";"),
    min_cells = min_cells,
    min_matched_pigs = min_pairs,
    eligible,
    contrast_status = if (eligible) "pending" else "skipped",
    skip_reason = reason
  )
  pb <- list(
    status = if (eligible) "eligible" else "no_eligible_contrast",
    reason = reason,
    contrast = contrast,
    sample_audit = audit,
    sample_info = samples,
    contrast_status = status,
    design_rank = NA_integer_,
    residual_df = NA_integer_,
    filter_min_group_size = length(pairs)
  )
  if (!eligible) {
    return(pb)
  }
  design <- model.matrix(
    ~ pig + condition,
    data = samples,
    contrasts.arg = list(pig = "contr.treatment", condition = "contr.treatment")
  )
  pb$design <- design
  pb$design_rank <- qr(design)$rank
  pb$residual_df <- nrow(design) - pb$design_rank
  # First-minus-second contrast: difference between the two design rows of one pig.
  same_pig <- as.character(samples$pig) == pairs[1]
  pb$contrast_vector <- design[same_pig & samples$condition == spec$group1, ] -
    design[same_pig & samples$condition == spec$group2, ]
  pb$counts <- population$counts[, samples$sample, drop = FALSE]
  y <- edgeR::DGEList(counts = pb$counts, samples = samples)
  keep <- edgeR::filterByExpr(
    y,
    group = samples$condition,
    min.count = deg_config$filter_min_count,
    min.total.count = deg_config$filter_min_total_count,
    large.n = deg_config$filter_large_n,
    min.prop = deg_config$filter_min_prop
  )
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::normLibSizes(y, method = "TMM")
  pb$sample_info$library_size <- y$samples$lib.size
  pb$sample_info$tmm_factor <- y$samples$norm.factors
  pb$sample_info$effective_library_size <- y$samples$lib.size *
    y$samples$norm.factors
  pb$keep <- keep
  pb$dge <- y
  pb
}

fit_pseudobulk <- function(pb, output_dir, gene_mapping) {
  fdr <- deg_config$deg_fdr
  pb$fit <- edgeR::glmQLFit(pb$dge, pb$design, legacy = FALSE, robust = FALSE)
  test <- edgeR::glmQLFTest(pb$fit, contrast = pb$contrast_vector)
  tab <- edgeR::topTags(
    test,
    n = Inf,
    adjust.method = "BH",
    sort.by = "PValue"
  )$table
  tab <- data.frame(gene = rownames(tab), tab, row.names = NULL)
  tab$direction <- ifelse(
    tab$logFC > 0,
    "up",
    ifelse(tab$logFC < 0, "down", "unchanged")
  )
  tab <- annotate_pig_genes(tab, gene_mapping)
  pb$results <- tab
  pb$status <- "completed"
  pb$contrast_status$contrast_status <- "tested"
  contrast <- pb$contrast
  write.csv(
    data.frame(gene = rownames(pb$counts), tested = pb$keep),
    file.path(output_dir, "gene_filter.csv"),
    row.names = FALSE
  )
  write.csv(
    tab,
    file.path(output_dir, paste0(contrast, "_all.csv")),
    row.names = FALSE
  )
  write.csv(
    tab[tab$FDR <= fdr, ],
    file.path(output_dir, paste0(contrast, "_significant.csv")),
    row.names = FALSE
  )
  pb$summary <- data.frame(
    contrast,
    tested = nrow(tab),
    up = sum(tab$FDR <= fdr & tab$logFC > 0),
    down = sum(tab$FDR <= fdr & tab$logFC < 0)
  )
  write.csv(
    merge(pb$contrast_status, pb$summary, by = "contrast", sort = FALSE),
    file.path(output_dir, "DEG_summary.csv"),
    row.names = FALSE
  )
  saveRDS(
    pb[c("contrast", "sample_info", "design", "contrast_vector", "dge", "fit")],
    file.path(output_dir, "pseudobulk_fit.rds")
  )
  pb
}
