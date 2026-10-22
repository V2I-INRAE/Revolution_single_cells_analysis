# Run from the project root: Rscript R/DEG/endothelial_cluster01_gsea_go_bp.R
suppressPackageStartupMessages({
  library(dplyr)
  library(fgsea)
  library(plotthis)
  library(ggplot2)
})
source("R/utils/save_plots.R")
source("R/DEG/publication_style.R")

parent_dir <- "results/job-44542861/clusters_plots/subset/Endothelial/job-44831031/DEG/cluster_01/gsea"
output_dir <- file.path(parent_dir, "go_bp")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
info_path <- file.path(output_dir, ".INFO")
info <- c(
  "# Endothelial cluster 1 GO Biological Process GSEA pilot",
  "Status: running", paste("Started:", Sys.time()),
  "Script: R/DEG/endothelial_cluster01_gsea_go_bp.R",
  paste("Producer: local Rscript; allocation", Sys.getenv("SLURM_JOB_ID", unset = "none")),
  "Parents: saved Hallmark/Reactome GSEA ranks, paired DE 44922357, Endothelial subclustering 44831031, cleanup 44827588, pipeline 44542861.",
  paste("Read-only input:", file.path(parent_dir, "gsea_results.rds")),
  "Reference: msigdbr(db_species=HS, species=Sus scrofa, collection=C5, subcollection=GO:BP), MSigDB 2026.1.Hs, frozen locally as go_bp_pig.rds. Human GO sets projected to pig orthologues, not native pig GO curation.",
  "Identifiers: reuse saved pig NCBI/Ensembl rank IDs unchanged. All 3638 unresolved SOURCE:<feature> IDs remain in the ranking but cannot match GO membership. No new identifier rescue, remapping, exclusions or concurrent orthologue outputs.",
  "Ranks: same 10163 features per contrast, sign(logFC)*sqrt(edgeR QL F), original order and eight approved collision exclusions. No gene FDR/logFC filter, refitting or downsampling.",
  "Method: fgseaMultilevel, matched minSize=15, maxSize=500, gseaParam=1, scoreType=std, eps=0, nproc=1, seed1234 reset per contrast; other parameters at package defaults.",
  "Multiple testing: BH separately across all eligible GO:BP terms in each contrast, separate from Hallmark/Reactome. Retain all results; reporting threshold pathway FDR<0.01.",
  "Direction: positive NES favours the first group; negative NES favours the second. No claim of pathway activation or exclusive induction. GO terms overlap and are not independent signatures.",
  "Display: up to eight terms per direction by FDR/p/name, regardless of significance; asterisk marks FDR<0.01. Lowest-FDR running-score example per direction, regardless of significance.",
  "Validation: pending saved-output, independent ranks/ES/leading-edge/BH and visual checks.",
  "Evidence: logs/endothelial-cluster01-gsea-go-bp-20261021.log",
  paste("R:", getRversion(), "; fgsea:", packageVersion("fgsea"),
    "; msigdbr:", packageVersion("msigdbr"), "; babelgene:", packageVersion("babelgene"),
    "; plotthis:", packageVersion("plotthis"), "; ggplot2:", packageVersion("ggplot2"))
)
writeLines(info, info_path)

parent <- readRDS(file.path(parent_dir, "gsea_results.rds"))
ranks <- parent$ranks
mapping <- parent$mapping
reference_path <- file.path(output_dir, "go_bp_pig.rds")
if (!file.exists(reference_path)) {
  reference <- msigdbr::msigdbr(db_species = "HS", species = "Sus scrofa",
    collection = "C5", subcollection = "GO:BP")
  saveRDS(reference, reference_path)
}
reference <- readRDS(reference_path)
stopifnot(all(reference$ortholog_taxon_id == 9823),
  identical(unique(reference$db_version), parent$database_version))
members <- reference %>%
  mutate(pig_gene = if_else(!is.na(ncbi_gene), paste0("NCBI:", ncbi_gene), ensembl_gene)) %>%
  select(gs_name, pig_gene) %>% distinct()
pathways <- split(members$pig_gene, members$gs_name)
matched_size <- vapply(pathways, function(x) sum(x %in% names(ranks[[1]])), integer(1))
coverage <- data.frame(pathway = names(matched_size), reference_size = lengths(pathways),
  matched_size = unname(matched_size), eligible = matched_size >= 15 & matched_size <= 500,
  row.names = NULL)
write.csv(coverage, file.path(output_dir, "pathway_coverage.csv"), row.names = FALSE)
cat("Reference:", unique(reference$db_version), "pig taxon", unique(reference$ortholog_taxon_id),
  "\nRanked features:", length(ranks[[1]]), "\nGO-matched features:",
  sum(names(ranks[[1]]) %in% unique(members$pig_gene)),
  "\nReference/eligible terms:", length(pathways), sum(coverage$eligible), "\n")

results <- plot_tables <- list()
for (contrast in names(ranks)) {
  set.seed(1234)
  result <- as.data.frame(fgseaMultilevel(pathways = pathways, stats = ranks[[contrast]],
    minSize = 15, maxSize = 500, gseaParam = 1, scoreType = "std", eps = 0, nproc = 1)) %>%
    arrange(padj, pval, pathway)
  groups <- strsplit(contrast, "_vs_", fixed = TRUE)[[1]]
  result$direction <- ifelse(result$NES > 0, groups[1], groups[2])
  result$leading_edge_source_genes <- vapply(result$leadingEdge, function(ids)
    paste(mapping$gene[match(ids, mapping$rank_id)], collapse = ";"), character(1))
  results[[contrast]] <- result
  csv <- result
  csv$leadingEdge <- vapply(csv$leadingEdge, paste, collapse = ";", character(1))
  write.csv(csv, file.path(output_dir, paste0(contrast, "_gsea.csv")), row.names = FALSE)
  cat("\nGO:BP", contrast, ":", nrow(result), "terms; FDR<0.01",
    sum(result$padj < 0.01, na.rm = TRUE), "\n")
  print(head(result[, c("pathway", "NES", "padj", "size", "direction")], 10), row.names = FALSE)

  display <- result %>% filter(is.finite(NES), !is.na(padj)) %>%
    group_by(direction) %>% slice_head(n = 8) %>% ungroup() %>% arrange(NES) %>%
    mutate(label = paste0(gsub("_", " ", sub("^GOBP_", "", pathway)),
      if_else(padj < 0.01, " *", "")), neg_log10_fdr = -log10(padj))
  display$label <- factor(display$label, levels = display$label)
  plot_tables[[contrast]] <- display
  write.csv(select(display, -leadingEdge),
    file.path(output_dir, paste0(contrast, "_displayed_pathways.csv")), row.names = FALSE)
  figure <- ScatterPlot(display, x = "NES", y = "label", size_by = "size",
    size_name = "Matched genes", color_by = "neg_log10_fdr", color_name = "-log10(BH FDR)",
    palette = "viridis", lower_cutoff = 0, upper_cutoff = max(display$neg_log10_fdr),
    shape = 21, alpha = 1, border_size = 0.3, aspect.ratio = NULL,
    theme = "theme_classic", theme_args = list(base_size = 10),
    title = "Endothelial cluster 1: GO Biological Process GSEA",
    subtitle = paste(groups[1], "minus", groups[2], "| Five paired pigs | * term FDR < 0.01"),
    xlab = paste("NES: negative favours", groups[2],
      "; positive favours", groups[1]), ylab = "", combine = FALSE, seed = 1234)[[1]] +
    scale_y_discrete(labels = function(x) vapply(x, function(label)
      paste(strwrap(label, width = 32), collapse = "\n"), character(1))) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.25) +
    deg_theme() +
      theme(legend.position = "bottom", legend.box = "vertical") +
      scale_size_continuous(range = c(2, 5)) +
      labs(caption = deg_wrap(paste0("Up to eight terms per direction by FDR; unstarred terms may not be significant. GO terms overlap.\n",
      "Pig orthologue projection of MSigDB ", unique(reference$db_version),
      ". BH across ", nrow(result), " GO:BP terms in this contrast.\n",
      length(ranks[[contrast]]), " ranked features; unchanged mapping; no gene FDR/logFC filtering.")))
  save_deg_figure(figure, output_dir, paste0(contrast, "_gsea_dotplot"), height = 230)

  top <- result %>% filter(is.finite(NES), !is.na(padj)) %>%
    group_by(direction) %>% slice_head(n = 1) %>% ungroup()
  running <- lapply(seq_len(nrow(top)), function(i) {
    fgsea::plotEnrichment(pathways[[top$pathway[i]]], ranks[[contrast]]) +
      labs(title = paste(strwrap(gsub("_", " ", top$pathway[i]), width = 34), collapse = "\n"),
        subtitle = sprintf("%s minus %s | NES %.2f\nBH FDR %.3g | favours %s",
          groups[1], groups[2], top$NES[i], top$padj[i], top$direction[i]),
        x = paste0("Rank: ", groups[1], " higher (left)\nto ", groups[2], " higher (right)"),
        y = "Running enrichment score") + deg_theme()
  })
  figure <- patchwork::wrap_plots(running, ncol = 2) +
    patchwork::plot_annotation(title = paste("GO:BP GSEA:", groups[1], "minus", groups[2]),
      caption = deg_wrap("Lowest-FDR term in each direction, regardless of significance. Pig orthologue projection; no claim of pathway activation."), theme = deg_theme())
  save_deg_figure(figure, output_dir, paste0(contrast, "_running_scores"), height = 140)
}

summary <- bind_rows(lapply(names(results), function(contrast) {
  result <- results[[contrast]]
  data.frame(collection = "GO:BP", contrast, ranked_features = length(ranks[[contrast]]),
    pathway_matched_features = sum(names(ranks[[contrast]]) %in% unique(members$pig_gene)),
    reference_pathways = length(pathways), tested_pathways = nrow(result),
    positive_FDR_lt_0.01 = sum(result$NES > 0 & result$padj < 0.01, na.rm = TRUE),
    negative_FDR_lt_0.01 = sum(result$NES < 0 & result$padj < 0.01, na.rm = TRUE),
    missing_pvalues = sum(is.na(result$pval)))
}))
write.csv(summary, file.path(output_dir, "coverage_summary.csv"), row.names = FALSE)
saveRDS(list(results = results, ranks = ranks, mapping = mapping, pathways = pathways,
  coverage = coverage, summary = summary, plot_tables = plot_tables,
  database_version = unique(reference$db_version), taxon = 9823L, parameters = parent$parameters),
  file.path(output_dir, "gsea_results.rds"))
info[info == "Status: running"] <- "Status: outputs generated; independent validation pending"
writeLines(c(info, paste("Finished:", Sys.time()),
  capture.output(print(summary, row.names = FALSE))), info_path)
print(summary, row.names = FALSE)
