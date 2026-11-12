source("R/plots/publications/publication_style.R")

pseudobulk_gsea_labels <- c(
  hallmark = "Hallmark",
  reactome = "Reactome",
  gobp = "GO:BP"
)

# Pig identifiers only; human orthologue names are not pig feature renames.
# Hallmark and Reactome share one alias union; GO:BP has an independent mapping.
prepare_pseudobulk_gsea <- function(common_genes, collection) {
  reference_dir <- deg_config$reference_dir
  hallmark <- readRDS(file.path(reference_dir, "hallmark_pig.rds"))
  reactome <- readRDS(file.path(reference_dir, "reactome_pig.rds"))
  selected <- if (collection == "hallmark") hallmark else reactome
  reference <- dplyr::bind_rows(hallmark, reactome)
  if (collection == "gobp") {
    selected <- readRDS(file.path(reference_dir, "gobp_pig.rds"))
    reference <- selected
  }
  reference$pig_gene <- ifelse(
    !is.na(reference$ncbi_gene),
    paste0("NCBI:", reference$ncbi_gene),
    reference$ensembl_gene
  )
  aliases <- unique(rbind(
    data.frame(alias = reference$gene_symbol, pig_gene = reference$pig_gene),
    data.frame(alias = reference$ensembl_gene, pig_gene = reference$pig_gene)
  ))
  aliases <- aliases[
    !is.na(aliases$alias) & nzchar(aliases$alias) & !is.na(aliases$pig_gene),
  ]
  ambiguous <- names(which(table(aliases$alias) > 1))
  lookup <- aliases[!aliases$alias %in% ambiguous, ]
  pig_names <- read.csv(deg_config$mapping_file)
  pig_names <- unique(pig_names[
    pig_names$pig_rename_status == "rename",
    c("input_gene_id", "pig_safe_name")
  ])
  mapping <- data.frame(
    gene = common_genes,
    safe_name = pig_names$pig_safe_name[match(
      common_genes,
      pig_names$input_gene_id
    )],
    pig_gene = lookup$pig_gene[match(common_genes, lookup$alias)]
  )
  mapping$match_method <- ifelse(
    !is.na(mapping$pig_gene),
    "exact_pig_identifier",
    ifelse(
      mapping$gene %in% ambiguous,
      "ambiguous_pig_identifier",
      "no_collection_match"
    )
  )
  rescue <- which(
    is.na(mapping$pig_gene) &
      !mapping$gene %in% ambiguous &
      !is.na(mapping$safe_name) &
      nzchar(mapping$safe_name)
  )
  mapping$pig_gene[rescue] <- lookup$pig_gene[match(
    mapping$safe_name[rescue],
    lookup$alias
  )]
  mapping$match_method[rescue[
    !is.na(mapping$pig_gene[rescue])
  ]] <- "existing_safe_name"
  mapping$display_gene <- ifelse(
    !is.na(mapping$safe_name) & nzchar(mapping$safe_name),
    mapping$safe_name,
    mapping$gene
  )
  # Approved exclusions: identifier collisions in all collections.
  mapping$excluded <- mapping$gene %in%
    c(
      "B2M",
      "ENSSSCG00000004687",
      "LDHB",
      "ENSSSCG00000000576",
      "ARFIP2",
      "TRIM3",
      "GNA12",
      "ENSSSCG00000032549"
    )
  if (collection == "gobp") {
    # Approved GO-only corrections; the unmodified reference assignment is kept for audit.
    mapping$reference_pig_gene <- mapping$pig_gene
    mapping$mapping_decision <- "reference_alias"
    il18bp <- mapping$gene == "ENSSSCG00000014800"
    mapping$pig_gene[il18bp] <- "NCBI:100525765"
    mapping$match_method[il18bp] <- "approved_Ensembl_NCBI_IL18BP"
    mapping$mapping_decision[il18bp] <- "approved_GO_only_IL18BP_correction"
    mapping$excluded_reason <- ifelse(
      mapping$excluded,
      "historical_approved_exclusion",
      ""
    )
    split_models <- mapping$gene %in%
      c(
        "ENSSSCG00000007481",
        "ZFP64",
        "ENSSSCG00000014098",
        "ENSSSCG00000051837"
      )
    mapping$excluded[split_models] <- TRUE
    mapping$excluded_reason[
      split_models
    ] <- "approved_GO_only_unresolved_split_models"
  }
  mapping$rank_id <- ifelse(
    !is.na(mapping$pig_gene),
    mapping$pig_gene,
    paste0("SOURCE:", mapping$gene)
  )
  if (collection == "gobp") {
    historical <- prepare_pseudobulk_gsea(common_genes, "hallmark")$mapping
    mapping$historical_rank_id <- historical$rank_id
    mapping$mapping_change <- ifelse(
      mapping$rank_id == historical$rank_id,
      "unchanged",
      ifelse(
        is.na(historical$pig_gene) & !is.na(mapping$pig_gene),
        "new_GO_mapping",
        ifelse(
          !is.na(historical$pig_gene) & is.na(mapping$pig_gene),
          "not_mapped_in_GO",
          "conflicting_canonical_ID"
        )
      )
    )
  }
  selected$pig_gene <- ifelse(
    !is.na(selected$ncbi_gene),
    paste0("NCBI:", selected$ncbi_gene),
    selected$ensembl_gene
  )
  members <- unique(selected[, c("gs_name", "pig_gene")])
  terms <- NULL
  if (collection == "gobp") {
    terms <- unique(as.data.frame(selected[, c(
      "gs_name",
      "gs_exact_source",
      "gs_id",
      "gs_description"
    )]))
    names(terms) <- c("pathway", "go_id", "msigdb_id", "description")
  }
  list(
    mapping = mapping,
    pathways = split(members$pig_gene, members$gs_name),
    terms = terms,
    database_version = unique(selected$db_version),
    collection = collection
  )
}

run_pseudobulk_gsea <- function(
  pb,
  reference,
  output_dir,
  label,
  excluded_features = character()
) {
  collection <- reference$collection
  collection_label <- pseudobulk_gsea_labels[[collection]]
  significant_only <- collection == "gobp"
  contrast <- pb$contrast
  seed <- deg_config$seed
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  mapping <- reference$mapping[
    match(rownames(pb$dge), reference$mapping$gene),
  ]
  mapping$excluded <- mapping$excluded | mapping$gene %in% excluded_features
  if (collection == "gobp") {
    mapping$excluded_reason[
      mapping$excluded & mapping$excluded_reason == ""
    ] <- "additional_historical_approved_exclusion"
  }
  write.csv(
    mapping,
    file.path(output_dir, "tested_gene_mapping.csv"),
    row.names = FALSE
  )
  write.csv(
    mapping[mapping$excluded, ],
    file.path(output_dir, "excluded_features.csv"),
    row.names = FALSE
  )
  ranked_mapping <- mapping[!mapping$excluded, ]
  # A GO:BP canonical ID that disagrees with Hallmark/Reactome would rank one gene under two identities.
  if (
    collection == "gobp" &&
      any(ranked_mapping$mapping_change == "conflicting_canonical_ID")
  ) {
    write.csv(
      ranked_mapping[
        ranked_mapping$mapping_change == "conflicting_canonical_ID",
      ],
      file.path(output_dir, "unresolved_mapping_conflicts.csv"),
      row.names = FALSE
    )
    return(list(
      status = "blocked: conflicting canonical identifiers",
      summary = NULL
    ))
  }
  # Two tested features with one canonical ID need an owner decision, never arbitrary rank selection.
  collisions <- duplicated(ranked_mapping$rank_id) |
    duplicated(ranked_mapping$rank_id, fromLast = TRUE)
  if (any(collisions)) {
    write.csv(
      ranked_mapping[collisions, ],
      file.path(output_dir, "unresolved_collisions.csv"),
      row.names = FALSE
    )
    return(list(
      status = "blocked: additional identifier collisions",
      summary = NULL
    ))
  }
  pathways <- reference$pathways
  matched_size <- vapply(
    pathways,
    function(x) sum(x %in% ranked_mapping$rank_id),
    integer(1)
  )
  coverage <- data.frame(
    pathway = names(pathways),
    reference_size = lengths(pathways),
    matched_size = unname(matched_size),
    eligible = matched_size >= deg_config$gsea_min_size &
      matched_size <= deg_config$gsea_max_size,
    row.names = NULL
  )
  if (!is.null(reference$terms)) {
    coverage <- cbind(
      coverage,
      reference$terms[match(coverage$pathway, reference$terms$pathway), -1]
    )
    write.csv(
      reference$terms,
      file.path(output_dir, "pathway_metadata.csv"),
      row.names = FALSE
    )
  }
  write.csv(
    coverage,
    file.path(output_dir, "pathway_coverage.csv"),
    row.names = FALSE
  )

  # Display annotations must never feed back into canonical IDs or rank joins.
  statistics <- pb$results[,
    setdiff(
      names(pb$results),
      c("gene_symbol", "display_gene", "mapping_source", "mapping_status")
    ),
    drop = FALSE
  ]
  tab <- merge(statistics, mapping, by = "gene", sort = FALSE)
  tab <- tab[!tab$excluded, ]
  # Near-zero deviance differences can be negative numerically; edgeR assigns P=1.
  tab$rank_statistic <- sign(tab$logFC) * sqrt(pmax(tab$F, 0))
  tab <- tab[order(-tab$rank_statistic, tab$gene), ]
  ranks <- setNames(tab$rank_statistic, tab$rank_id)
  write.csv(
    tab,
    file.path(output_dir, paste0(contrast, "_ranked_genes.csv")),
    row.names = FALSE
  )

  run_fgsea <- function(n_perm) {
    set.seed(seed)
    as.data.frame(fgsea::fgseaMultilevel(
      pathways,
      ranks,
      minSize = deg_config$gsea_min_size,
      maxSize = deg_config$gsea_max_size,
      gseaParam = 1,
      scoreType = "std",
      eps = 0,
      nproc = 1,
      nPermSimple = n_perm
    ))
  }
  unavailable <- function(x) !is.finite(x$NES) | is.na(x$pval) | is.na(x$padj)
  result <- run_fgsea(deg_config$gsea_nperm)
  initial_unavailable <- unavailable(result)
  retried <- any(initial_unavailable)
  if (retried) {
    write.csv(
      result[initial_unavailable, c("pathway", "NES", "pval", "padj")],
      file.path(output_dir, paste0(contrast, "_initial_unavailable.csv")),
      row.names = FALSE
    )
    message(
      collection_label,
      " retry for ",
      contrast,
      ": nPermSimple=",
      deg_config$gsea_nperm_retry,
      "; initial unavailable pathways=",
      sum(initial_unavailable)
    )
    result <- run_fgsea(deg_config$gsea_nperm_retry)
  }
  write.csv(
    data.frame(
      contrast,
      initial_unavailable = sum(initial_unavailable),
      retried,
      nPermSimple = if (retried) {
        deg_config$gsea_nperm_retry
      } else {
        deg_config$gsea_nperm
      },
      final_unavailable = sum(unavailable(result)),
      seed
    ),
    file.path(output_dir, "retry_audit.csv"),
    row.names = FALSE
  )
  if (any(unavailable(result))) {
    stop(
      collection_label,
      " enrichment still has unavailable estimates in ",
      label,
      " ",
      contrast,
      " after the retry; no term deletion."
    )
  }

  result <- result[order(result$padj, result$pval, result$pathway), ]
  if (!is.null(reference$terms)) {
    result <- cbind(
      result,
      reference$terms[match(result$pathway, reference$terms$pathway), -1]
    )
  }
  groups <- strsplit(contrast, "_vs_", fixed = TRUE)[[1]]
  result$direction <- ifelse(result$NES > 0, groups[1], groups[2])
  result$leading_edge_source_genes <- vapply(
    result$leadingEdge,
    function(ids) {
      paste(mapping$gene[match(ids, mapping$rank_id)], collapse = ";")
    },
    character(1)
  )
  result$leading_edge_display_genes <- vapply(
    result$leadingEdge,
    function(ids) {
      paste(mapping$display_gene[match(ids, mapping$rank_id)], collapse = ";")
    },
    character(1)
  )
  csv <- result
  csv$leadingEdge <- vapply(
    csv$leadingEdge,
    paste,
    collapse = ";",
    character(1)
  )
  write.csv(
    csv,
    file.path(output_dir, paste0(contrast, "_gsea.csv")),
    row.names = FALSE
  )

  pathway_fdr <- deg_config$pathway_fdr
  display <- result
  if (significant_only) {
    display <- display[display$padj < pathway_fdr, ]
  }
  display <- as.data.frame(dplyr::ungroup(dplyr::slice_head(
    dplyr::group_by(display, direction),
    n = 8
  )))
  display <- display[order(display$NES), ]
  display$label <- if (nrow(display)) {
    factor(paste0(
      gsub(
        "_",
        " ",
        sub(paste0("^", toupper(collection), "_"), "", display$pathway)
      ),
      ifelse(display$padj < pathway_fdr, " *", "")
    ))
  } else {
    factor(character())
  }
  display$label <- factor(display$label, levels = as.character(display$label))
  display$neg_log10_fdr <- -log10(display$padj)
  write.csv(
    display[, setdiff(names(display), "leadingEdge")],
    file.path(output_dir, paste0(contrast, "_displayed_pathways.csv")),
    row.names = FALSE
  )
  render_pseudobulk_gsea(
    result,
    display,
    ranks,
    pathways,
    reference$database_version,
    pb$sample_info,
    contrast,
    output_dir,
    label,
    collection_label,
    significant_only
  )

  summary <- data.frame(
    contrast,
    input_features = nrow(mapping),
    excluded_features = sum(mapping$excluded),
    ranked_features = length(ranks),
    pathway_matched_features = sum(names(ranks) %in% unique(unlist(pathways))),
    reference_pathways = length(pathways),
    tested_pathways = nrow(result),
    significant_positive = sum(result$NES > 0 & result$padj < pathway_fdr),
    significant_negative = sum(result$NES < 0 & result$padj < pathway_fdr)
  )
  write.csv(
    summary,
    file.path(output_dir, "coverage_summary.csv"),
    row.names = FALSE
  )
  list(status = "completed", summary = summary)
}

render_pseudobulk_gsea <- function(
  result,
  display,
  ranks,
  pathways,
  database_version,
  sample_info,
  contrast,
  output_dir,
  label,
  collection_label,
  significant_only
) {
  pathway_fdr <- deg_config$pathway_fdr
  groups <- strsplit(contrast, "_vs_", fixed = TRUE)[[1]]
  n_samples <- vapply(
    groups,
    function(group) sum(sample_info$condition == group),
    integer(1)
  )
  n_pigs <- length(unique(sample_info$pig))
  figure <- if (!nrow(display)) {
    ggplot2::ggplot() +
      ggplot2::annotate(
        "text",
        x = 0,
        y = 0,
        label = paste0(
          "No significant ",
          collection_label,
          " terms (FDR < ",
          pathway_fdr,
          ")"
        )
      ) +
      ggplot2::theme_void() +
      ggplot2::labs(
        title = paste(label, collection_label, "GSEA"),
        subtitle = contrast
      )
  } else {
    plotthis::ScatterPlot(
      display,
      x = "NES",
      y = "label",
      size_by = "size",
      size_name = "Matched genes",
      color_by = "neg_log10_fdr",
      color_name = "-log10(BH FDR)",
      palette = "viridis",
      lower_cutoff = 0,
      upper_cutoff = max(display$neg_log10_fdr),
      shape = 21,
      alpha = 1,
      border_size = 0.3,
      aspect.ratio = NULL,
      theme = "theme_classic",
      theme_args = list(base_size = 10),
      title = paste(label, collection_label, "GSEA"),
      subtitle = sprintf(
        "%s (n=%d samples) minus %s (n=%d samples) | %d paired pigs | * pathway FDR < %g",
        groups[1],
        n_samples[1],
        groups[2],
        n_samples[2],
        n_pigs,
        pathway_fdr
      ),
      xlab = paste(
        "NES: negative favours",
        groups[2],
        "; positive favours",
        groups[1]
      ),
      ylab = "",
      combine = FALSE,
      seed = deg_config$seed
    )[[1]] +
      ggplot2::scale_y_discrete(labels = function(x) {
        vapply(
          x,
          function(z) {
            paste(strwrap(z, width = 32), collapse = "\n")
          },
          character(1)
        )
      }) +
      ggplot2::geom_vline(
        xintercept = 0,
        linetype = "dashed",
        colour = "grey50",
        linewidth = 0.25
      ) +
      deg_theme() +
      ggplot2::theme(legend.position = "bottom", legend.box = "vertical") +
      ggplot2::scale_size_continuous(range = c(2, 5)) +
      ggplot2::labs(
        caption = deg_wrap(paste0(
          if (significant_only) {
            "Up to eight significant terms per direction by FDR.\n"
          } else {
            "Up to eight terms per direction by FDR; unstarred terms may not be significant.\n"
          },
          "Pig orthologue projection of MSigDB ",
          database_version,
          ". BH across ",
          nrow(result),
          " pathways in this collection/cluster/contrast.\n",
          length(ranks),
          " ranked features; no gene FDR/logFC filtering."
        ))
      )
  }
  if (!nrow(display)) {
    figure <- figure +
      ggplot2::labs(
        caption = deg_wrap(paste0(
          "MSigDB ",
          database_version,
          "; BH across ",
          nrow(result),
          " terms; ",
          length(ranks),
          " ranked features. All tested terms remain in tables."
        ))
      )
  }
  # Long term names need room; keep the historical Hallmark dimensions.
  dotplot_height <- 170
  if (collection_label != "Hallmark" && nrow(display)) {
    label_lines <- vapply(
      as.character(display$label),
      function(x) length(strwrap(x, width = 32)),
      integer(1)
    )
    dotplot_height <- max(
      dotplot_height,
      90 + 4 * nrow(display) * max(label_lines)
    )
  }
  save_deg_figure(
    figure,
    output_dir,
    paste0(contrast, "_gsea_dotplot"),
    height = dotplot_height
  )
  candidates <- result
  if (significant_only) {
    candidates <- candidates[candidates$padj < pathway_fdr, ]
  }
  top <- as.data.frame(dplyr::ungroup(dplyr::slice_head(
    dplyr::group_by(candidates, direction),
    n = 1
  )))
  running <- lapply(seq_len(nrow(top)), function(i) {
    fgsea::plotEnrichment(pathways[[top$pathway[i]]], ranks) +
      ggplot2::labs(
        title = paste(
          strwrap(
            gsub(
              "_",
              " ",
              if (collection_label == "GO:BP") {
                sub("^GOBP_", "", top$pathway[i])
              } else {
                top$pathway[i]
              }
            ),
            width = if (collection_label == "GO:BP") 24 else 34
          ),
          collapse = "\n"
        ),
        subtitle = sprintf(
          "%s minus %s | NES %.2f\nBH FDR %.3g | favours %s",
          groups[1],
          groups[2],
          top$NES[i],
          top$padj[i],
          top$direction[i]
        ),
        x = paste0(
          "Rank: ",
          groups[1],
          " higher (left)\nto ",
          groups[2],
          " higher (right)"
        ),
        y = "Running enrichment score"
      ) +
      deg_theme()
  })
  if (!length(running)) {
    running <- list(
      ggplot2::ggplot() +
        ggplot2::annotate(
          "text",
          x = 0,
          y = 0,
          label = paste0(
            "No significant ",
            collection_label,
            " terms (FDR < ",
            pathway_fdr,
            ")"
          )
        ) +
        ggplot2::theme_void()
    )
  }
  figure <- patchwork::wrap_plots(running, ncol = 2) +
    patchwork::plot_annotation(
      title = paste(
        label,
        paste0(collection_label, ":"),
        groups[1],
        "minus",
        groups[2]
      ),
      caption = deg_wrap(paste(
        if (significant_only) {
          "Lowest-FDR significant term in each available direction."
        } else {
          "Lowest-FDR pathway in each direction, regardless of significance."
        },
        "Pig orthologue projection; no claim of pathway activation."
      )),
      theme = deg_theme()
    )
  save_deg_figure(
    figure,
    output_dir,
    paste0(contrast, "_running_scores"),
    height = if (collection_label == "GO:BP") 140 else 110
  )
}
