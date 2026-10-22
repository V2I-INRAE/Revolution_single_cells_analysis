# Run from the project root: Rscript R/DEG/endothelial_cluster01_ora_top100.R
# Add --all to use all qualifying genes and save separately from the top-100 pilot.
# Add --relaxed for all genes at FDR 0.05 / absolute log2FC 0.5 and current pig names.
suppressPackageStartupMessages({
  library(dplyr)
  library(msigdbr)
  library(plotthis)
  library(ggplot2)
})
source("R/utils/save_plots.R")
source("R/DEG/publication_style.R")

input_dir <- "results/job-44542861/clusters_plots/subset/Endothelial/job-44831031/DEG/cluster_01"
relaxed <- "--relaxed" %in% commandArgs(trailingOnly = TRUE)
all_genes <- relaxed || "--all" %in% commandArgs(trailingOnly = TRUE)
fdr_cutoff <- if (relaxed) 0.05 else 0.01
logfc_cutoff <- if (relaxed) 0.5 else 1
selection_tag <- if (relaxed) "all_fdr0.05_logfc0.5" else if (all_genes) "all" else "top100"
selection_label <- if (all_genes) "All qualifying genes" else "Top 100 per direction"
output_dir <- file.path(input_dir, paste0("ora_", selection_tag))
mapping_file <- if (relaxed) "results/misc/Ensembl_ID/gene-mapping/pig_to_human_orthologues_annotated.csv" else
  if (all_genes) file.path(input_dir, "ora_top100", "ora_results.rds") else
  "results/misc/Ensembl_ID/gene-mapping/rename_map.csv"
reference_path <- file.path(input_dir, "ora_top100", "hallmark_pig.rds")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
info_path <- file.path(output_dir, ".INFO")
info <- c(
  paste("# Endothelial cluster 1 directional ORA:", selection_label),
  "Status: running", paste("Started:", Sys.time()),
  paste("Script: R/DEG/endothelial_cluster01_ora_top100.R",
    if (relaxed) "--relaxed" else if (all_genes) "--all" else ""),
  paste("Producer: local Rscript; allocation", Sys.getenv("SLURM_JOB_ID", unset = "none")),
  paste("Read-only inputs:", file.path(input_dir, "P4_vs_N4_all.csv"), ";", mapping_file, ";", reference_path),
  "Parents: paired DE allocation 44922357; Endothelial subclustering 44831031; cleanup 44827588; pipeline 44542861.",
  "Comparison: paired P4-N4 at T4H; positive logFC means higher in P4, negative means higher in N4. Not a test of exclusive induction relative to Control.",
  paste(sprintf("Selection: gene BH FDR<=%g and abs(logFC)>=%g;", fdr_cutoff, logfc_cutoff), selection_label,
    "ordered by FDR, PValue, gene ID, before matching; no replacement of unmatched genes."),
  "Species: Sus scrofa, taxon 9823; msigdbr(db_species='HS', species='Sus scrofa', collection='H'); human Hallmark sets projected to pig orthologues, not native pig pathways.",
  "Reference: reuse hallmark_pig.rds if present, otherwise retrieve and save it; retain original orthology provenance in the snapshot.",
  if (relaxed) "Mapping: current annotated orthologue table, distinct input_gene_id/pig_safe_name/pig_rename_status only. Human targets are not used with pig-projected sets; pig safe-name decisions are unchanged from the prior mapping." else
    if (all_genes) "Mapping: reuse the top100 pilot's saved matching; no new annotation or orthology updates. The only analysis change is removing the top100 cap.",
  "Matching: exact pig Ensembl ID or pig symbol; fallback to existing rename_status=rename safe_name only. No human-symbol matching or source-ID changes.",
  "Unit: unique pig NCBI gene ID, or pig Ensembl ID when NCBI is absent. Count multiple source aliases once; exclude aliases matching more than one pig gene, without choosing or expanding them.",
  "Background: all tested genes with unambiguous Hallmark membership, using the same matching and unique-gene rule as each query.",
  "Method: stats::phyper upper-tail hypergeometric ORA; BH separately across all 50 terms per direction, including zero overlaps. No pathway-size or pathway-significance filter.",
  "Display: ten overlapping terms per direction ranked by pathway FDR/PValue/name, regardless of significance; fold enrichment x, overlap size, -log10(BH FDR) colour. All 50 terms saved.",
  "Seed: no stochastic analysis; plotthis seed 1234. No DE refitting or source-object edits.",
  "Dependencies: msigdbr/babelgene installed locally through renv; renv.lock unchanged.",
  "Validation: one saved-output and visual review only; no repeated numerical audits requested.",
  paste0("Evidence: logs/endothelial-cluster01-ora-", selection_tag, "-20261021.log"),
  paste("R:", getRversion(), "; msigdbr:", packageVersion("msigdbr"),
    "; babelgene:", packageVersion("babelgene"), "; dplyr:", packageVersion("dplyr"),
    "; plotthis:", packageVersion("plotthis"), "; ggplot2:", packageVersion("ggplot2"))
)
writeLines(info, info_path)

if (file.exists(reference_path)) {
  hallmark <- readRDS(reference_path)
} else {
  hallmark <- msigdbr(db_species = "HS", species = "Sus scrofa", collection = "H")
  saveRDS(hallmark, reference_path)
}
stopifnot(all(hallmark$ortholog_taxon_id == 9823))
de <- read.csv(file.path(input_dir, "P4_vs_N4_all.csv"))

# Keep unnamed pig genes with Ensembl IDs, and never expand ambiguous symbols.
reference <- hallmark %>% mutate(pig_gene = if_else(!is.na(ncbi_gene),
  paste0("NCBI:", ncbi_gene), ensembl_gene))
if (all_genes && !relaxed) {
  mapping <- readRDS(mapping_file)$mapping
} else {
  rename_map <- if (relaxed) read.csv(mapping_file) %>%
    distinct(ensembl_id = input_gene_id, safe_name = pig_safe_name, rename_status = pig_rename_status) else
    read.csv(mapping_file)
  aliases <- bind_rows(
    reference %>% transmute(alias = gene_symbol, pig_gene),
    reference %>% transmute(alias = ensembl_gene, pig_gene)
  ) %>% filter(!is.na(alias), !is.na(pig_gene)) %>% distinct()
  ambiguous <- aliases %>% count(alias) %>% filter(n > 1) %>% pull(alias)
  lookup <- aliases %>% filter(!alias %in% ambiguous)
  mapping <- data.frame(gene = de$gene,
    safe_name = rename_map$safe_name[match(de$gene,
      ifelse(rename_map$rename_status == "rename", rename_map$ensembl_id, NA))],
    pig_gene = lookup$pig_gene[match(de$gene, lookup$alias)])
  mapping$match_method <- ifelse(!is.na(mapping$pig_gene), "exact_pig_identifier",
    ifelse(mapping$gene %in% ambiguous, "ambiguous_pig_identifier", "no_Hallmark_match"))
  rescue <- which(is.na(mapping$pig_gene) & !mapping$gene %in% ambiguous &
    !is.na(mapping$safe_name))
  mapping$pig_gene[rescue] <- lookup$pig_gene[match(mapping$safe_name[rescue], lookup$alias)]
  mapping$match_method[rescue[!is.na(mapping$pig_gene[rescue])]] <- "existing_safe_name"
}
background <- unique(na.omit(mapping$pig_gene))
term2gene <- reference %>% select(gs_name, pig_gene) %>%
  filter(pig_gene %in% background) %>% distinct()
terms <- sort(unique(reference$gs_name))
write.csv(mapping, file.path(output_dir, "tested_gene_mapping.csv"), row.names = FALSE)

selected <- results <- list()
summary_rows <- list()
for (direction in c("N4_higher", "P4_higher")) {
  sign <- if (direction == "N4_higher") -1 else 1
  eligible <- de %>% filter(FDR <= fdr_cutoff, logFC * sign >= logfc_cutoff) %>% arrange(FDR, PValue, gene)
  query_table <- if (all_genes) eligible else head(eligible, 100)
  query_table <- left_join(query_table, mapping, by = "gene")
  query_table$display_gene <- ifelse(!is.na(query_table$safe_name) & nzchar(query_table$safe_name),
    query_table$safe_name, query_table$gene)
  query <- unique(na.omit(query_table$pig_gene))
  result <- bind_rows(lapply(terms, function(term) {
    members <- term2gene$pig_gene[term2gene$gs_name == term]
    shared <- intersect(query, members)
    k <- length(shared)
    n <- length(query)
    m <- length(members)
    total <- length(background)
    data.frame(pathway = term, overlap = k, query_size = n, pathway_size = m,
      background_size = total, gene_ratio = k / n, background_ratio = m / total,
      fold_enrichment = (k / n) / (m / total),
      PValue = phyper(k - 1, m, total - m, n, lower.tail = FALSE),
      pig_gene_ids = paste(sort(shared), collapse = ";"),
      source_genes = paste(query_table$gene[query_table$pig_gene %in% shared], collapse = ";"))
  })) %>% mutate(FDR = p.adjust(PValue, method = "BH"), .after = PValue) %>%
    arrange(FDR, PValue, pathway)
  selected[[direction]] <- query_table
  results[[direction]] <- result
  summary_rows[[direction]] <- data.frame(direction, eligible_features = nrow(eligible),
    selected_features = nrow(query_table), matched_features = sum(!is.na(query_table$pig_gene)),
    unique_query_genes = length(query), tested_features = nrow(de),
    matched_background_features = sum(!is.na(mapping$pig_gene)),
    unique_background_genes = length(background), tested_pathways = nrow(result))
  write.csv(query_table, file.path(output_dir, paste0(selection_tag, "_", direction, ".csv")), row.names = FALSE)
  write.csv(result, file.path(output_dir, paste0("hallmark_ora_", direction, ".csv")), row.names = FALSE)
}
coverage <- bind_rows(summary_rows)
write.csv(coverage, file.path(output_dir, "coverage_summary.csv"), row.names = FALSE)
saveRDS(list(selected = selected, results = results, mapping = mapping,
  background = background, term2gene = term2gene, coverage = coverage,
  database_version = unique(hallmark$db_version), taxon = 9823L,
  parameters = list(gene_fdr = fdr_cutoff, abs_logfc = logfc_cutoff, mapping_file = mapping_file)),
  file.path(output_dir, "ora_results.rds"))

plot_data <- bind_rows(results, .id = "direction") %>% filter(overlap > 0) %>%
  group_by(direction) %>% slice_head(n = 10) %>% ungroup() %>%
  mutate(neg_log10_fdr = -log10(FDR),
    label = gsub("_", " ", sub("^HALLMARK_", "", pathway)))
colour_max <- max(plot_data$neg_log10_fdr)
x_max <- max(plot_data$fold_enrichment) * 1.1
plots <- lapply(names(results), function(direction) {
  tab <- plot_data[plot_data$direction == direction, ]
  tab$label <- factor(tab$label, levels = rev(tab$label))
  n <- coverage$unique_query_genes[coverage$direction == direction]
  total_selected <- coverage$selected_features[coverage$direction == direction]
  p <- ScatterPlot(tab, x = "fold_enrichment", y = "label",
    color_by = "neg_log10_fdr", color_name = "-log10(BH FDR)",
    size_by = "overlap", size_name = "Overlap (genes)",
    palette = "viridis", lower_cutoff = 0, upper_cutoff = colour_max,
    shape = 21, alpha = 1, border_size = 0.3,
    theme = "theme_classic", theme_args = list(base_size = 10),
    aspect.ratio = NULL, combine = FALSE, seed = 1234,
    title = if (direction == "N4_higher") "Higher in N4 than P4" else "Higher in P4 than N4",
    subtitle = paste(n, "of", total_selected, "selected genes have usable Hallmark membership"),
    xlab = "Fold enrichment (observed / expected)", ylab = "")[[1]] +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60", linewidth = 0.25) +
    scale_y_discrete(labels = function(x) vapply(x, deg_wrap, character(1), width = 32)) +
    scale_size_area(max_size = 5, limits = c(0, max(plot_data$overlap)),
      breaks = if (all_genes) c(10, 20, 30) else c(5, 10)) +
    coord_cartesian(xlim = c(0, x_max)) + deg_theme() +
    theme(legend.box = "vertical")
  p
})
figure <- patchwork::wrap_plots(plots, ncol = 1, guides = "collect") +
  patchwork::plot_annotation(title = "Endothelial cluster 1: directional Hallmark ORA at T4H",
    subtitle = sprintf("%s from paired P4-N4; gene FDR <= %g and |log2 FC| >= %g",
      selection_label, fdr_cutoff, logfc_cutoff),
    caption = deg_wrap(paste0("Pig orthologue projection of human MSigDB ", unique(hallmark$db_version),
      "; background: ", length(background), " tested, annotated pig genes.\n",
      "Ten terms per direction ranked by pathway FDR; displayed terms are not necessarily significant. BH correction covers all 50 terms per direction.")), theme = deg_theme())
save_deg_figure(figure, output_dir, paste0("hallmark_ora_", selection_tag), height = 230)

info[info == "Status: running"] <- "Status: outputs generated; independent validation pending"
writeLines(c(info, paste("Finished:", Sys.time()),
  paste("MSigDB release:", paste(unique(hallmark$db_version), collapse = ","),
    "; pig reference rows:", nrow(hallmark), "; Hallmark terms:", length(terms)),
  paste("Background:", sum(!is.na(mapping$pig_gene)), "matched features ->", length(background), "unique pig genes."),
  paste("Ambiguous tested identifiers excluded from enrichment:",
    paste(mapping$gene[mapping$match_method == "ambiguous_pig_identifier"], collapse = ", ")),
  "Genes without Hallmark membership are not necessarily missing general gene annotation; the mapping CSV records collection-specific matching.",
  capture.output(print(coverage, row.names = FALSE))), info_path)
print(coverage, row.names = FALSE)
for (direction in names(results)) {
  cat("\n", direction, "\n")
  print(head(results[[direction]][, c("pathway", "overlap", "fold_enrichment", "FDR", "source_genes")], 10), row.names = FALSE)
}
