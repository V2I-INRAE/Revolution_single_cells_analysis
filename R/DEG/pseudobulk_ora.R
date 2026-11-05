# Over-representation among all FDR-significant genes, against tested genes.
# References use the same pig identifiers and frozen collection as GSEA, but
# GSEA-specific rank exclusions are not automatically applied to the ORA list.
run_pseudobulk_ora <- function(pb, reference, output_dir, label, fdr_cutoff = 0.05) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  mapping <- reference$mapping[match(rownames(pb$dge), reference$mapping$gene), , drop = FALSE]
  if (anyNA(mapping$gene)) stop("ORA mapping does not cover all tested features")
  mapping$ora_reason <- ifelse(mapping$match_method == "ambiguous_pig_identifier", "ambiguous_identifier", "")
  if ("mapping_change" %in% names(mapping)) {
    mapping$ora_reason[mapping$mapping_change == "conflicting_canonical_ID"] <- "conflicting_canonical_ID"
  }
  mapping$ora_included <- !nzchar(mapping$ora_reason)
  write.csv(mapping, file.path(output_dir, "tested_gene_mapping.csv"), row.names = FALSE)
  retained <- mapping[mapping$ora_included, , drop = FALSE]
  # An approved gene-level collapse: one background gene per canonical ID;
  # a gene enters the query when any mapped tested feature is significant.
  collisions <- duplicated(retained$rank_id) | duplicated(retained$rank_id, fromLast = TRUE)
  write.csv(retained[collisions, , drop = FALSE], file.path(output_dir, "collapsed_features.csv"), row.names = FALSE)
  universe <- unique(retained$rank_id)
  if (!length(universe)) stop("No unambiguous tested genes for ORA")
  pathways <- lapply(reference$pathways, intersect, y = universe)
  sizes <- lengths(pathways)
  pathways <- pathways[sizes >= 15L & sizes <= 500L]
  if (!length(pathways)) stop("No ORA pathways have 15–500 tested genes")
  N <- length(universe)
  summary <- list()
  for (contrast in names(pb$results)) {
    tab <- pb$results[[contrast]]
    if (anyDuplicated(tab$gene) || !setequal(tab$gene, mapping$gene)) stop("ORA DEG and tested gene IDs disagree")
    significant <- tab$gene[!is.na(tab$FDR) & tab$FDR <= fdr_cutoff]
    query <- retained$rank_id[retained$gene %in% significant]
    query <- unique(query)
    n <- length(query)
    K <- lengths(pathways)
    overlap <- vapply(pathways, function(x) sum(x %in% query), integer(1))
    p <- if (n) stats::phyper(overlap - 1L, K, N - K, n, lower.tail = FALSE) else rep(1, length(K))
    results <- data.frame(pathway = names(pathways), universe_size = N, pathway_size = K,
      query_size = n, overlap = overlap, fold_enrichment = if (n) (overlap / n) / (K / N) else NA_real_,
      pvalue = p, FDR = stats::p.adjust(p, method = "BH"),
      overlap_genes = vapply(pathways, function(x) paste(retained$gene[retained$rank_id %in% intersect(x, query)], collapse = ";"), character(1)),
      row.names = NULL)
    results <- results[order(results$FDR, results$pvalue, results$pathway), ]
    write.csv(results, file.path(output_dir, paste0(contrast, "_ora.csv")), row.names = FALSE)
    shown <- head(results[results$FDR <= fdr_cutoff & results$overlap > 0, ], 20L)
    fig <- if (nrow(shown)) {
      ggplot2::ggplot(shown, ggplot2::aes(x = fold_enrichment, y = stats::reorder(pathway, fold_enrichment), size = overlap, colour = FDR)) +
        ggplot2::geom_point() + ggplot2::scale_colour_viridis_c(direction = -1) +
        ggplot2::labs(x = "Fold enrichment", y = NULL, title = paste(label, contrast), subtitle = paste(reference$collection, "ORA; FDR <=", fdr_cutoff), size = "Overlap") + ggplot2::theme_bw()
    } else {
      ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0, label = "No significant enriched pathways") + ggplot2::theme_void() + ggplot2::labs(title = paste(label, contrast))
    }
    ggplot2::ggsave(file.path(output_dir, paste0(contrast, "_ora.pdf")), fig, width = 10, height = 8, units = "in", bg = "white")
    summary[[contrast]] <- data.frame(contrast, tested_features = nrow(mapping), excluded_ambiguous = sum(!mapping$ora_included),
      universe = N, collapsed_source_features = sum(collisions), significant_source_genes = length(significant), significant_in_universe = n,
      pathways_tested = nrow(results), pathways_FDR_le_0.05 = sum(results$FDR <= fdr_cutoff))
  }
  summary <- do.call(rbind, summary)
  write.csv(summary, file.path(output_dir, "coverage_summary.csv"), row.names = FALSE)
  saveRDS(list(summary = summary, universe = universe, pathways = pathways, mapping = mapping,
    collection = reference$collection, database_version = reference$database_version,
    reference_files = reference$reference_files, fdr_cutoff = fdr_cutoff), file.path(output_dir, "ora_results.rds"))
  list(status = "completed", summary = summary)
}
