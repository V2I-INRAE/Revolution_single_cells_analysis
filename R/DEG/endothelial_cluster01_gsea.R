# Run from the project root: Rscript R/DEG/endothelial_cluster01_gsea.R
suppressPackageStartupMessages({
  library(dplyr)
  library(fgsea)
  library(plotthis)
  library(ggplot2)
})
source("R/utils/save_plots.R")

input_dir <- "results/job-44542861/clusters_plots/subset/Endothelial/job-44831031/DEG/cluster_01"
output_dir <- file.path(input_dir, "gsea")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
info_path <- file.path(output_dir, ".INFO")
info <- c(
  "# Endothelial cluster 1 separate Hallmark and Reactome GSEA pilot",
  "Status: running", paste("Started:", Sys.time()),
  "Preparation: started 2026-10-21 16:16; Reactome reference and identifier audit saved before analysis.",
  "Script: R/DEG/endothelial_cluster01_gsea.R",
  paste("Producer: local Rscript; allocation", Sys.getenv("SLURM_JOB_ID", unset = "none")),
  "Parents: paired DE allocation 44922357; Endothelial subclustering 44831031; cleanup 44827588; pipeline 44542861.",
  paste("Read-only inputs:", input_dir,
    "P4_vs_N4_all.csv, P4_vs_Control_all.csv, N4_vs_Control_all.csv; ora_top100/ora_results.rds frozen aliases; ora_top100/hallmark_pig.rds; gsea/reactome_pig.rds."),
  "Species: Sus scrofa, taxon 9823; MSigDB human Hallmark and C2/CP:REACTOME projected to pig orthologues using msigdbr. Not native pig pathway curation.",
  "Reference: frozen Hallmark/Reactome snapshots, MSigDB 2026.1.Hs; no new online annotation or concurrent orthologue-run outputs used.",
  "Matching: combined-reference exact pig symbols/Ensembl IDs; fallback to frozen safe_name from existing rename_status=rename records. Canonical pig NCBI ID, or Ensembl ID when NCBI absent; no human-ID matching.",
  "Exclusions: owner approved all eight features in four collision pairs on 2026-10-21: B2M/ENSSSCG00000004687, LDHB/ENSSSCG00000000576, ARFIP2/TRIM3, GNA12/ENSSSCG00000032549. No arbitrary rank aggregation.",
  "Unmatched/ambiguous identifiers: retain as distinct SOURCE:<original feature> ranks, never assign arbitrary gene identities. Ambiguous RABIF/IFIT3 remain unmatched. Unannotated features contribute to the ranked universe but cannot be pathway hits.",
  "Ranks: all remaining tested features, sign(logFC)*sqrt(edgeR QL F), decreasing order with original gene ID as tie-breaker. No gene FDR/logFC filter, downsampling, refitting or jitter.",
  "Method: fgseaMultilevel, minSize=15, maxSize=500 after matching to ranks, gseaParam=1, scoreType=std, eps=0, nproc=1, seed1234 reset for each collection/contrast; other algorithm parameters at package defaults.",
  "Multiple testing: BH independently per collection and contrast across all size-eligible pathways. All results retained; padj is pathway FDR. No correction across collections/contrasts.",
  "Direction: positive NES favours the first group in P4_vs_N4/P4_vs_Control/N4_vs_Control; negative NES favours the second. No claim of exclusive induction, pathway activation or independent validation.",
  "Display: up to eight pathways per NES direction by padj/pval/name, regardless of significance; an asterisk marks pathway FDR<0.01. One native running-score plot per direction, selected by the same rule. All eligible pathways remain in CSV/RDS.",
  "Dependencies: fgsea/msigdbr/babelgene installed locally through renv; renv.lock unchanged.",
  "Validation: pending saved CSV/RDS, independent ranking/ES/BH and deterministic rerun checks, PNG/PDF inspection.",
  "Evidence: logs/endothelial-cluster01-gsea-20261021.log",
  paste("R:", getRversion(), "; fgsea:", packageVersion("fgsea"),
    "; msigdbr:", packageVersion("msigdbr"), "; babelgene:", packageVersion("babelgene"),
    "; plotthis:", packageVersion("plotthis"), "; ggplot2:", packageVersion("ggplot2"))
)
writeLines(info, info_path)

references <- list(
  Hallmark = readRDS(file.path(input_dir, "ora_top100", "hallmark_pig.rds")),
  Reactome = readRDS(file.path(output_dir, "reactome_pig.rds"))
)
stopifnot(all(vapply(references, function(x) all(x$ortholog_taxon_id == 9823), logical(1))))
reference <- bind_rows(references, .id = "collection") %>%
  mutate(pig_gene = if_else(!is.na(ncbi_gene), paste0("NCBI:", ncbi_gene), ensembl_gene))
aliases <- bind_rows(
  reference %>% transmute(alias = gene_symbol, pig_gene),
  reference %>% transmute(alias = ensembl_gene, pig_gene)
) %>% filter(!is.na(alias), alias != "", !is.na(pig_gene)) %>% distinct()
ambiguous <- aliases %>% count(alias) %>% filter(n > 1) %>% pull(alias)
lookup <- aliases %>% filter(!alias %in% ambiguous)
mapping <- readRDS(file.path(input_dir, "ora_top100", "ora_results.rds"))$mapping %>%
  select(gene, safe_name) %>%
  mutate(pig_gene = lookup$pig_gene[match(gene, lookup$alias)],
    match_method = if_else(!is.na(pig_gene), "exact_pig_identifier",
      if_else(gene %in% ambiguous, "ambiguous_pig_identifier", "no_collection_match")))
rescue <- which(is.na(mapping$pig_gene) & !mapping$gene %in% ambiguous & !is.na(mapping$safe_name))
mapping$pig_gene[rescue] <- lookup$pig_gene[match(mapping$safe_name[rescue], lookup$alias)]
mapping$match_method[rescue[!is.na(mapping$pig_gene[rescue])]] <- "existing_safe_name"
excluded_features <- c("B2M", "ENSSSCG00000004687", "LDHB", "ENSSSCG00000000576",
  "ARFIP2", "TRIM3", "GNA12", "ENSSSCG00000032549")
mapping <- mapping %>% mutate(excluded = gene %in% excluded_features,
  rank_id = if_else(!is.na(pig_gene), pig_gene, paste0("SOURCE:", gene)))
stopifnot(!anyDuplicated(mapping$rank_id[!mapping$excluded]))
write.csv(mapping, file.path(output_dir, "tested_gene_mapping.csv"), row.names = FALSE)
write.csv(filter(mapping, excluded), file.path(output_dir, "excluded_features.csv"), row.names = FALSE)

contrasts <- c("P4_vs_N4", "P4_vs_Control", "N4_vs_Control")
rank_tables <- ranks <- list()
for (contrast in contrasts) {
  tab <- read.csv(file.path(input_dir, paste0(contrast, "_all.csv"))) %>%
    left_join(mapping, by = "gene") %>% filter(!excluded) %>%
    mutate(rank_statistic = sign(logFC) * sqrt(F)) %>% arrange(desc(rank_statistic), gene)
  rank_tables[[contrast]] <- tab
  ranks[[contrast]] <- setNames(tab$rank_statistic, tab$rank_id)
  write.csv(tab, file.path(output_dir, paste0(contrast, "_ranked_genes.csv")), row.names = FALSE)
}

results <- pathways <- coverage <- plot_tables <- list()
for (collection in names(references)) {
  directory <- file.path(output_dir, tolower(collection))
  dir.create(directory, showWarnings = FALSE)
  members <- reference %>% filter(.data$collection == .env$collection) %>%
    select(gs_name, pig_gene) %>% distinct()
  pathways[[collection]] <- split(members$pig_gene, members$gs_name)
  matched_size <- vapply(pathways[[collection]], function(x) sum(x %in% names(ranks[[1]])), integer(1))
  cov <- data.frame(pathway = names(matched_size),
    reference_size = lengths(pathways[[collection]]), matched_size = unname(matched_size),
    eligible = matched_size >= 15 & matched_size <= 500, row.names = NULL)
  coverage[[collection]] <- cov
  write.csv(cov, file.path(directory, "pathway_coverage.csv"), row.names = FALSE)

  for (contrast in contrasts) {
    set.seed(1234)
    result <- as.data.frame(fgseaMultilevel(pathways = pathways[[collection]],
      stats = ranks[[contrast]], minSize = 15, maxSize = 500,
      gseaParam = 1, scoreType = "std", eps = 0, nproc = 1)) %>%
      arrange(padj, pval, pathway)
    result$direction <- ifelse(result$NES > 0,
      sub("_vs_.*", "", contrast), sub(".*_vs_", "", contrast))
    result$leading_edge_source_genes <- vapply(result$leadingEdge, function(ids)
      paste(mapping$gene[match(ids, mapping$rank_id)], collapse = ";"), character(1))
    results[[collection]][[contrast]] <- result
    csv <- result
    csv$leadingEdge <- vapply(csv$leadingEdge, paste, collapse = ";", character(1))
    write.csv(csv, file.path(directory, paste0(contrast, "_gsea.csv")), row.names = FALSE)
    cat("\n", collection, contrast, ":", nrow(result), "pathways; FDR<0.01",
      sum(result$padj < 0.01, na.rm = TRUE), "\n")
    print(head(result[, c("pathway", "NES", "padj", "size", "direction")], 10), row.names = FALSE)

    display <- result %>% filter(is.finite(NES), !is.na(padj)) %>%
      group_by(direction) %>% slice_head(n = 8) %>% ungroup() %>% arrange(NES) %>%
      mutate(label = paste0(gsub("_", " ", sub("^(HALLMARK|REACTOME)_", "", pathway)),
        if_else(padj < 0.01, " *", "")), neg_log10_fdr = -log10(padj))
    display$label <- factor(display$label, levels = display$label)
    plot_tables[[collection]][[contrast]] <- display
    write.csv(select(display, -leadingEdge),
      file.path(directory, paste0(contrast, "_displayed_pathways.csv")), row.names = FALSE)
    groups <- strsplit(contrast, "_vs_", fixed = TRUE)[[1]]
    figure <- ScatterPlot(display, x = "NES", y = "label", size_by = "size",
      size_name = "Matched genes", color_by = "neg_log10_fdr", color_name = "-log10(BH FDR)",
      palette = "viridis", lower_cutoff = 0, upper_cutoff = max(display$neg_log10_fdr),
      shape = 21, alpha = 1, border_size = 0.3, aspect.ratio = NULL,
      theme = "theme_classic", theme_args = list(base_size = 11),
      title = paste("Endothelial cluster 1:", collection, "GSEA"),
      subtitle = paste(groups[1], "minus", groups[2], "| Five paired pigs | * pathway FDR < 0.01"),
      xlab = paste("Normalized enrichment score (NES): negative favours", groups[2],
        "; positive favours", groups[1]), ylab = "", combine = FALSE, seed = 1234)[[1]] +
      scale_y_discrete(labels = function(x) vapply(x, function(label)
        paste(strwrap(label, width = 55), collapse = "\n"), character(1))) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
      labs(caption = paste0("Up to eight terms per direction by FDR; unstarred terms may not be significant.\n",
        "Pig orthologue projection of MSigDB ", unique(reference$db_version),
        ". BH across ", nrow(result), " pathways in this collection/contrast.\n",
        length(ranks[[contrast]]), " ranked features; no gene FDR/logFC filtering."))
    save_figure(figure, directory, paste0(contrast, "_gsea_dotplot"), width = 12.5, height = 10)
    ggsave(file.path(directory, paste0(contrast, "_gsea_dotplot.pdf")), figure,
      width = 12.5, height = 10, bg = "white")

    top <- result %>% filter(is.finite(NES), !is.na(padj)) %>%
      group_by(direction) %>% slice_head(n = 1) %>% ungroup()
    running <- lapply(seq_len(nrow(top)), function(i) {
      fgsea::plotEnrichment(pathways[[collection]][[top$pathway[i]]], ranks[[contrast]]) +
        labs(title = paste(strwrap(gsub("_", " ", top$pathway[i]), width = 53), collapse = "\n"),
          subtitle = sprintf("%s minus %s | NES %.2f | BH FDR %.3g | favours %s",
            groups[1], groups[2], top$NES[i], top$padj[i], top$direction[i]),
          x = paste("Rank:", groups[1], "higher (left) to", groups[2], "higher (right)"),
          y = "Running enrichment score") +
        theme_classic(base_size = 11)
    })
    figure <- patchwork::wrap_plots(running, ncol = 2) +
      patchwork::plot_annotation(title = paste(collection, "GSEA:", groups[1], "minus", groups[2]),
        caption = "Lowest-FDR pathway in each direction, regardless of significance. Pig orthologue projection; no claim of pathway activation.")
    save_figure(figure, directory, paste0(contrast, "_running_scores"), width = 14, height = 5)
    ggsave(file.path(directory, paste0(contrast, "_running_scores.pdf")), figure,
      width = 14, height = 5, bg = "white")
  }
}

summary <- bind_rows(lapply(names(results), function(collection) {
  bind_rows(lapply(contrasts, function(contrast) {
    result <- results[[collection]][[contrast]]
    data.frame(collection, contrast, input_features = nrow(mapping),
      excluded_features = sum(mapping$excluded), ranked_features = length(ranks[[contrast]]),
      pathway_matched_features = sum(names(ranks[[contrast]]) %in% unique(unlist(pathways[[collection]]))),
      reference_pathways = length(pathways[[collection]]), tested_pathways = nrow(result),
      positive_FDR_lt_0.01 = sum(result$NES > 0 & result$padj < 0.01, na.rm = TRUE),
      negative_FDR_lt_0.01 = sum(result$NES < 0 & result$padj < 0.01, na.rm = TRUE),
      missing_pvalues = sum(is.na(result$pval)))
  }))
}))
write.csv(summary, file.path(output_dir, "coverage_summary.csv"), row.names = FALSE)
saveRDS(list(results = results, ranks = ranks, rank_tables = rank_tables, mapping = mapping,
  pathways = pathways, coverage = coverage, summary = summary, plot_tables = plot_tables,
  database_version = unique(reference$db_version), taxon = 9823L,
  parameters = list(minSize = 15L, maxSize = 500L, gseaParam = 1, scoreType = "std",
    eps = 0, nproc = 1L, seed = 1234L)), file.path(output_dir, "gsea_results.rds"))
info[info == "Status: running"] <- "Status: outputs generated; independent validation pending"
writeLines(c(info, paste("Finished:", Sys.time()),
  capture.output(print(summary, row.names = FALSE))), info_path)
print(summary, row.names = FALSE)
