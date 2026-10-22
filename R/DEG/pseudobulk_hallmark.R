source("R/DEG/publication_style.R")

# Pig identifiers only. Human orthologue names are not pig feature renames.
prepare_hallmark <- function(common_genes, mapping_file, reference_dir) {
  hallmark <- readRDS(file.path(reference_dir, "hallmark_pig.rds"))
  reference <- dplyr::bind_rows(hallmark, readRDS(file.path(reference_dir, "reactome_pig.rds")))
  stopifnot(all(reference$ortholog_taxon_id == 9823))
  reference$pig_gene <- ifelse(!is.na(reference$ncbi_gene),
    paste0("NCBI:", reference$ncbi_gene), reference$ensembl_gene)
  aliases <- unique(rbind(
    data.frame(alias = reference$gene_symbol, pig_gene = reference$pig_gene),
    data.frame(alias = reference$ensembl_gene, pig_gene = reference$pig_gene)))
  aliases <- aliases[!is.na(aliases$alias) & nzchar(aliases$alias) & !is.na(aliases$pig_gene), ]
  ambiguous <- names(which(table(aliases$alias) > 1))
  lookup <- aliases[!aliases$alias %in% ambiguous, ]
  pig_names <- read.csv(mapping_file)
  pig_names <- unique(pig_names[pig_names$pig_rename_status == "rename",
    c("input_gene_id", "pig_safe_name")])
  mapping <- data.frame(gene = common_genes,
    safe_name = pig_names$pig_safe_name[match(common_genes, pig_names$input_gene_id)],
    pig_gene = lookup$pig_gene[match(common_genes, lookup$alias)])
  mapping$match_method <- ifelse(!is.na(mapping$pig_gene), "exact_pig_identifier",
    ifelse(mapping$gene %in% ambiguous, "ambiguous_pig_identifier", "no_collection_match"))
  rescue <- which(is.na(mapping$pig_gene) & !mapping$gene %in% ambiguous &
    !is.na(mapping$safe_name) & nzchar(mapping$safe_name))
  mapping$pig_gene[rescue] <- lookup$pig_gene[match(mapping$safe_name[rescue], lookup$alias)]
  mapping$match_method[rescue[!is.na(mapping$pig_gene[rescue])]] <- "existing_safe_name"
  mapping$display_gene <- ifelse(!is.na(mapping$safe_name) & nzchar(mapping$safe_name),
    mapping$safe_name, mapping$gene)
  mapping$excluded <- mapping$gene %in% c("B2M", "ENSSSCG00000004687", "LDHB",
    "ENSSSCG00000000576", "ARFIP2", "TRIM3", "GNA12", "ENSSSCG00000032549")
  mapping$rank_id <- ifelse(!is.na(mapping$pig_gene), mapping$pig_gene,
    paste0("SOURCE:", mapping$gene))
  hallmark$pig_gene <- ifelse(!is.na(hallmark$ncbi_gene),
    paste0("NCBI:", hallmark$ncbi_gene), hallmark$ensembl_gene)
  members <- unique(hallmark[, c("gs_name", "pig_gene")])
  list(mapping = mapping, pathways = split(members$pig_gene, members$gs_name),
    database_version = unique(hallmark$db_version))
}

run_pseudobulk_hallmark <- function(pb, reference, output_dir, label,
    pathway_fdr = 0.01, seed = 1234, excluded_features = character()) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  mapping <- reference$mapping[match(rownames(pb$dge), reference$mapping$gene), ]
  mapping$excluded <- mapping$excluded | mapping$gene %in% excluded_features
  write.csv(mapping, file.path(output_dir, "tested_gene_mapping.csv"), row.names = FALSE)
  write.csv(mapping[mapping$excluded, ], file.path(output_dir, "excluded_features.csv"), row.names = FALSE)
  ranked_mapping <- mapping[!mapping$excluded, ]
  collisions <- duplicated(ranked_mapping$rank_id) | duplicated(ranked_mapping$rank_id, fromLast = TRUE)
  # Additional collision pairs need a scientific decision, never arbitrary rank selection.
  if (any(collisions)) {
    write.csv(ranked_mapping[collisions, ], file.path(output_dir, "unresolved_collisions.csv"), row.names = FALSE)
    return(list(status = "blocked: additional identifier collisions", summary = NULL))
  }
  unlink(file.path(output_dir, "unresolved_collisions.csv"))
  pathways <- reference$pathways
  matched_size <- vapply(pathways, function(x) sum(x %in% ranked_mapping$rank_id), integer(1))
  coverage <- data.frame(pathway = names(pathways), reference_size = lengths(pathways),
    matched_size = unname(matched_size), eligible = matched_size >= 15 & matched_size <= 500,
    row.names = NULL)
  write.csv(coverage, file.path(output_dir, "pathway_coverage.csv"), row.names = FALSE)
  results <- ranks <- rank_tables <- plot_tables <- list()
  for (contrast in names(pb$results)) {
    tab <- merge(pb$results[[contrast]], mapping, by = "gene", sort = FALSE)
    tab <- tab[!tab$excluded, ]
    # Near-zero deviance differences can be negative numerically; edgeR assigns P=1.
    tab$rank_statistic <- sign(tab$logFC) * sqrt(pmax(tab$F, 0))
    tab <- tab[order(-tab$rank_statistic, tab$gene), ]
    rank_tables[[contrast]] <- tab
    ranks[[contrast]] <- setNames(tab$rank_statistic, tab$rank_id)
    write.csv(tab, file.path(output_dir, paste0(contrast, "_ranked_genes.csv")), row.names = FALSE)
    set.seed(seed)
    result <- as.data.frame(fgsea::fgseaMultilevel(pathways, ranks[[contrast]],
      minSize = 15, maxSize = 500, gseaParam = 1, scoreType = "std", eps = 0, nproc = 1))
    result <- result[order(result$padj, result$pval, result$pathway), ]
    groups <- strsplit(contrast, "_vs_", fixed = TRUE)[[1]]
    result$direction <- ifelse(result$NES > 0, groups[1], groups[2])
    result$leading_edge_source_genes <- vapply(result$leadingEdge, function(ids)
      paste(mapping$gene[match(ids, mapping$rank_id)], collapse = ";"), character(1))
    result$leading_edge_display_genes <- vapply(result$leadingEdge, function(ids)
      paste(mapping$display_gene[match(ids, mapping$rank_id)], collapse = ";"), character(1))
    results[[contrast]] <- result
    csv <- result
    csv$leadingEdge <- vapply(csv$leadingEdge, paste, collapse = ";", character(1))
    write.csv(csv, file.path(output_dir, paste0(contrast, "_gsea.csv")), row.names = FALSE)
    display <- result[is.finite(result$NES) & !is.na(result$padj), ]
    display <- dplyr::slice_head(dplyr::group_by(display, direction), n = 8)
    display <- as.data.frame(dplyr::ungroup(display))
    display <- display[order(display$NES), ]
    display$label <- factor(paste0(gsub("_", " ", sub("^HALLMARK_", "", display$pathway)),
      ifelse(display$padj < pathway_fdr, " *", "")))
    display$label <- factor(display$label, levels = as.character(display$label))
    display$neg_log10_fdr <- -log10(display$padj)
    plot_tables[[contrast]] <- display
    write.csv(display[, setdiff(names(display), "leadingEdge")],
      file.path(output_dir, paste0(contrast, "_displayed_pathways.csv")), row.names = FALSE)
    render_pseudobulk_hallmark(result, display, ranks[[contrast]], pathways,
      reference$database_version, nlevels(pb$sample_info$pig), contrast, output_dir,
      label, pathway_fdr, seed)
  }
  summary <- do.call(rbind, lapply(names(results), function(contrast) {
    x <- results[[contrast]]
    data.frame(contrast, input_features = nrow(mapping), excluded_features = sum(mapping$excluded),
      ranked_features = length(ranks[[contrast]]),
      pathway_matched_features = sum(names(ranks[[contrast]]) %in% unique(unlist(pathways))),
      reference_pathways = length(pathways), tested_pathways = nrow(x),
      positive_FDR_lt_0.01 = sum(x$NES > 0 & x$padj < pathway_fdr, na.rm = TRUE),
      negative_FDR_lt_0.01 = sum(x$NES < 0 & x$padj < pathway_fdr, na.rm = TRUE),
      missing_pvalues = sum(is.na(x$pval)))
  }))
  write.csv(summary, file.path(output_dir, "coverage_summary.csv"), row.names = FALSE)
  saveRDS(list(results = results, ranks = ranks, rank_tables = rank_tables, mapping = mapping,
    pathways = pathways, coverage = coverage, summary = summary, plot_tables = plot_tables,
    database_version = reference$database_version, taxon = 9823L,
    parameters = list(minSize = 15L, maxSize = 500L, gseaParam = 1, scoreType = "std",
      eps = 0, nproc = 1L, seed = seed, pathway_fdr = pathway_fdr,
      excluded_features = excluded_features)),
    file.path(output_dir, "gsea_results.rds"))
  list(status = "completed", summary = summary)
}

# Plot existing enrichment results; this function never runs fgseaMultilevel.
render_pseudobulk_hallmark <- function(result, display, ranks, pathways,
    database_version, n_pigs, contrast, output_dir, label, pathway_fdr = 0.01, seed = 1234) {
  groups <- strsplit(contrast, "_vs_", fixed = TRUE)[[1]]
  figure <- plotthis::ScatterPlot(display, x = "NES", y = "label", size_by = "size",
    size_name = "Matched genes", color_by = "neg_log10_fdr", color_name = "-log10(BH FDR)",
    palette = "viridis", lower_cutoff = 0, upper_cutoff = max(display$neg_log10_fdr),
    shape = 21, alpha = 1, border_size = 0.3, aspect.ratio = NULL,
    theme = "theme_classic", theme_args = list(base_size = 10),
    title = paste(label, "Hallmark GSEA"),
    subtitle = paste(groups[1], "minus", groups[2], "|", n_pigs,
      "paired pigs | * pathway FDR <", pathway_fdr),
    xlab = paste("NES: negative favours", groups[2], "; positive favours", groups[1]),
    ylab = "", combine = FALSE, seed = seed)[[1]] +
    ggplot2::scale_y_discrete(labels = function(x) vapply(x, function(z)
      paste(strwrap(z, width = 32), collapse = "\n"), character(1))) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.25) +
    deg_theme() +
    ggplot2::theme(legend.position = "bottom", legend.box = "vertical") +
    ggplot2::scale_size_continuous(range = c(2, 5)) +
    ggplot2::labs(caption = deg_wrap(paste0("Up to eight terms per direction by FDR; unstarred terms may not be significant.\n",
      "Pig orthologue projection of MSigDB ", database_version,
      ". BH across ", nrow(result), " pathways in this cluster/contrast.\n",
      length(ranks), " ranked features; no gene FDR/logFC filtering.")))
  save_deg_figure(figure, output_dir, paste0(contrast, "_gsea_dotplot"), height = 170)
  top <- as.data.frame(dplyr::ungroup(dplyr::slice_head(dplyr::group_by(
    result[is.finite(result$NES) & !is.na(result$padj), ], direction), n = 1)))
  running <- lapply(seq_len(nrow(top)), function(i) {
    fgsea::plotEnrichment(pathways[[top$pathway[i]]], ranks) +
      ggplot2::labs(title = paste(strwrap(gsub("_", " ", top$pathway[i]), width = 34), collapse = "\n"),
        subtitle = sprintf("%s minus %s | NES %.2f\nBH FDR %.3g | favours %s",
          groups[1], groups[2], top$NES[i], top$padj[i], top$direction[i]),
        x = paste0("Rank: ", groups[1], " higher (left)\nto ", groups[2], " higher (right)"),
        y = "Running enrichment score") + deg_theme()
  })
  figure <- patchwork::wrap_plots(running, ncol = 2) +
    patchwork::plot_annotation(title = paste(label, "Hallmark:", groups[1], "minus", groups[2]),
      caption = deg_wrap("Lowest-FDR pathway in each direction, regardless of significance. Pig orthologue projection; no claim of pathway activation."), theme = deg_theme())
  save_deg_figure(figure, output_dir, paste0(contrast, "_running_scores"), height = 110)
}
