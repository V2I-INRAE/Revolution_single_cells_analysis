# Paired pseudobulk DEG, GSEA and ORA orchestrator.
# CLI: Rscript R/DEG/main.R --help. Inputs and new output root are explicit;
# frozen 2026.1 pig pathway references are shared and read-only.
suppressPackageStartupMessages({
  library(Seurat)
  library(edgeR)
  library(ggplot2)
})
source("R/utils/save_plots.R")
source("R/DEG/pseudobulk_de.R")
source("R/DEG/pseudobulk_diagnostics.R")
source("R/DEG/pseudobulk_plots.R")
source("R/DEG/pseudobulk_gsea.R")
source("R/DEG/pseudobulk_ora.R")
source("R/plots/publications/plot_DEG_for_pub.R")
source("R/plots/publications/collect_figures.R")

run_pseudobulk <- function(input, subset_label, cluster_col, clusters, output_root,
    mapping_file = "results/misc/Ensembl_ID/gene-mapping/pig_to_human_orthologues_annotated.csv",
    reference_dir = "results/misc/msigdb_pig_2026.1", gsea_exclusions = list(),
    pigs = NULL, min_matched_pigs = 3L, log_path = NULL,
    collections = c("hallmark", "reactome", "gobp")) {
  specs <- pseudobulk_gsea_collections()
  if (!length(collections) || anyNA(collections) || anyDuplicated(collections) ||
      !all(collections %in% specs$collection)) stop("Select unique collections from hallmark, reactome, gobp.")
  if (length(min_matched_pigs) != 1L || !is.finite(min_matched_pigs) ||
      min_matched_pigs < 2L || min_matched_pigs %% 1 != 0) {
    stop("min_matched_pigs must be an integer >=2; single-pig inference is not supported.")
  }
  if (!is.null(pigs) && (!length(pigs) || anyNA(pigs) || any(!nzchar(pigs)) || anyDuplicated(pigs))) {
    stop("Supply a nonempty, unique pig list, or NULL for all pigs.")
  }
  if (missing(output_root) || length(output_root) != 1L || is.na(output_root) || !nzchar(output_root)) {
    stop("Supply an explicit new output_root; historical analysis outputs must not be overwritten.")
  }
  output_root <- path.expand(output_root)
  if (file.exists(output_root)) stop("output_root already exists; choose a new analysis directory.")
  object <- readRDS(input)
  metadata <- object[[]]
  required <- c("sample", "pig", "pressure", "time_point", cluster_col)
  if (!all(required %in% names(metadata))) stop("Required sample/cluster metadata columns are missing.")
  if (anyDuplicated(colnames(object)) || anyDuplicated(rownames(metadata)) ||
      !identical(colnames(object), rownames(metadata))) {
    stop("Cell IDs must be unique and exactly aligned with metadata rows.")
  }
  if (anyNA(metadata[, required, drop = FALSE]) ||
      any(vapply(metadata[, required, drop = FALSE], function(x) any(!nzchar(as.character(x))), logical(1)))) {
    stop("Sample/cluster metadata contains missing or empty values.")
  }
  if (!length(clusters) || anyDuplicated(as.character(clusters)) ||
      any(!as.character(clusters) %in% as.character(metadata[[cluster_col]]))) {
    stop("Requested clusters must be unique, observed cluster labels.")
  }
  # Validate original tuples before relabeling, so T8H/T10H collisions cannot be merged.
  metadata$time_point_original <- as.character(metadata$time_point)
  sample_panel <- unique(metadata[, c("sample", "pig", "pressure", "time_point", "time_point_original")])
  if (anyDuplicated(sample_panel$sample)) stop("Each sample must have exactly one original metadata tuple.")
  sample_panel$sample <- as.character(sample_panel$sample)
  sample_panel$time_point <- as.character(sample_panel$time_point)
  sample_panel$time_point[sample_panel$time_point == "T8H"] <- "T10H"
  condition_levels <- c("Control", "P4", "N4", "P10", "N10")
  condition_map <- c(T0H_None = "Control", T4H_positive = "P4", T4H_negative = "N4",
    T10H_positive = "P10", T10H_negative = "N10")
  sample_panel$condition <- factor(unname(condition_map[paste(sample_panel$time_point,
    sample_panel$pressure, sep = "_")]), levels = condition_levels)
  if (anyNA(sample_panel$condition)) stop("Unexpected time/pressure combination in the observed panel.")
  if (anyDuplicated(sample_panel[, c("pig", "condition")])) {
    stop("Multiple original samples map to one pig-condition; no collisions may be merged.")
  }
  sample_panel$pig <- factor(sample_panel$pig)
  rownames(sample_panel) <- sample_panel$sample
  if (anyNA(match(as.character(metadata$sample), sample_panel$sample))) stop("Cell-to-sample metadata join failed.")
  object$time_point_original <- metadata$time_point_original
  # Assign a complete character vector, not a new value into a factor.
  object$time_point <- as.character(metadata$time_point)
  object$time_point <- replace(object$time_point, object$time_point == "T8H", "T10H")
  count_layers <- SeuratObject::Layers(object[["RNA"]], search = "^counts[.]")
  if (!length(count_layers)) stop("Original RNA counts.* layers are required.")
  layer_cells <- lapply(count_layers, function(layer) SeuratObject::Cells(object[["RNA"]], layer = layer))
  raw_cells <- unlist(layer_cells, use.names = FALSE)
  if (anyDuplicated(raw_cells) || !setequal(raw_cells, colnames(object))) {
    stop("Original count-layer cell IDs must cover every metadata cell exactly once.")
  }
  for (i in seq_along(count_layers)) {
    layer_samples <- unique(as.character(metadata[layer_cells[[i]], "sample"]))
    if (length(layer_samples) != 1L || !identical(layer_samples, sub("^counts[.]", "", count_layers[i]))) {
      stop("Original count layer and sample metadata disagree: ", count_layers[i])
    }
  }
  if (!setequal(sub("^counts[.]", "", count_layers), sample_panel$sample)) {
    stop("Original count layers do not match the full observed specimen panel.")
  }
  common_genes <- Reduce(intersect, lapply(count_layers, function(layer)
    SeuratObject::Features(object[["RNA"]], layer = layer)))
  if (!length(common_genes) || anyDuplicated(common_genes)) stop("Invalid original common-layer gene intersection.")
  if (!is.null(pigs) && !all(pigs %in% as.character(sample_panel$pig))) {
    stop("Unknown pigs: ", paste(setdiff(pigs, as.character(sample_panel$pig)), collapse = ", "))
  }
  # Select specimens only AFTER validating all source layers and fixing the gene universe.
  cohort_audit <- sample_panel
  cohort_audit$selected <- if (is.null(pigs)) TRUE else as.character(sample_panel$pig) %in% pigs
  cohort_audit$exclusion_reason <- ifelse(cohort_audit$selected, "", "pig_not_selected")
  sample_panel <- sample_panel[cohort_audit$selected, , drop = FALSE]
  sample_panel$pig <- droplevels(sample_panel$pig)
  if (nlevels(sample_panel$pig) < min_matched_pigs) {
    stop("Selected cohort has fewer pigs than min_matched_pigs; choose the intended pair minimum explicitly.")
  }
  if (!dir.create(output_root, recursive = TRUE, showWarnings = FALSE)) stop("Could not create a new output_root.")
  info_path <- file.path(output_root, ".INFO")
  if (is.null(log_path)) log_path <- paste0("logs/", tolower(subset_label), "-",
    tolower(basename(output_root)), "-", format(Sys.Date(), "%Y%m%d"), ".log")
  info <- c(paste("#", subset_label, "paired pseudobulk pipeline"), "Status: running",
    paste("Started:", Sys.time()), paste("Producer allocation:", Sys.getenv("SLURM_JOB_ID", unset = "local")),
    "Script: R/DEG/main.R; count/DE, diagnostics, DE and publication plotting, GSEA and ORA.",
    paste("Selected GSEA collections:", paste(collections, collapse = ", ")),
    paste("Input (read-only):", input), paste("Source producer:", basename(dirname(input))),
    paste("Population:", subset_label, "; cluster column:", cluster_col, "; clusters:", paste(clusters, collapse = ", ")),
    "Groups (ordered): Control/T0H, P4/T4H positive, N4/T4H negative, P10/T10H positive, N10/T10H negative airway-pressure ventilation.",
    paste("Cohort:", paste(levels(sample_panel$pig), collapse = ", "), "; selected/source specimens:",
      nrow(sample_panel), "/", nrow(cohort_audit), "; Control/P4/N4/P10/N10:",
      paste(as.integer(table(sample_panel$condition)), collapse = "/")),
    "Cohort exclusions: cohort_audit.csv records every source specimen and pig_not_selected reasons.",
    if (min_matched_pigs == 2L) "Interpretation: owner-approved post-hoc exploratory testing with >=2 matched pigs; biological replication remains limited.",
    "T8H becomes T10H character-safely in memory only; time_point_original and original sample IDs (including REVO31-N8) retained for audit.",
    "Exclude only individual observed-panel aggregates with <10 cells (including absent cluster populations); sample_audit.csv retains exact original IDs, n_cells and exclusion reasons. No zero-expression libraries.",
    "Retain ALL passing specimens in one joint model, including late groups with only 1 or 2 replicates. No imputation, additional deletion or downsampling.",
    paste0("Test each requested contrast only with >=", min_matched_pigs,
      " DIRECTLY MATCHED PIGS after sample filtering; missing groups or fewer pairs are explicitly skipped. All six requests/counts/pairs/reasons appear in contrast_status.csv and run summaries."),
    "No eligible contrasts => no fit. Rank-deficient design or residual df<=0 => blocked, no fit; no model simplification or further sample removal.",
    "Candidates: intersection of ALL original RNA counts.* layers, before sample/time/cluster selection.",
    "Filtering per cluster: filterByExpr(group=condition), default min.count=10/min.total.count=15/large.n=10/min.prop=0.7. Observed retained groups only; a retained singleton reduces minimum group size to 1. This consequence is retained, not silently changed.",
    "Normalization: recalculate filtered library sizes, then TMM. No integrated or normalized cell counts used.",
    "Model: ~ pig + condition; edgeR glmQLFit(legacy=FALSE), dispersion=NULL, robust=FALSE; glmQLFTest.",
    "Requested contrasts (ordered, first-minus-second; intercept/pig coefficients zero): P4_vs_Control, N4_vs_Control, P4_vs_N4, P4_vs_P10, N4_vs_N10, P10_vs_N10. One joint fit on observed condition levels per eligible cluster; full rank and positive residual df required. BH separately per cluster/tested contrast.",
    "DE exports: full statistics and FDR<=0.05 tables, without an effect-size cutoff.",
    "Volcano/heatmap selection: FDR<=0.05 AND abs(log2FC)>=0.5 within the same tested contrast; heatmap union across tested members of the six requested contrasts.",
    "MDS: same filtered TMM DGEList as DE; top500 pairwise leading-logFC distances, prior.count=2; no pig-effect removal.",
    "Heatmap: TMM log2 CPM prior.count=2, row z-scores over all samples, Euclidean/complete row clustering, fixed Control/P4/N4/P10/N10 then pig order; RdBu blue-low/white-zero/red-high saturated at +/-2.5; pig effects retained.",
    "Volcano: exact logFC and -log10(BH FDR), four classes, +/-0.5 and FDR0.05 dashed guides; five labels per direction.",
    paste("Pig display-name mapping:", mapping_file, "; pig_safe_name only when pig_rename_status=rename. Original feature IDs retained; no human rename."),
    paste("References:", reference_dir, "; frozen MSigDB pig projections, taxon9823. H/Reactome remain the historical DE-display alias source; GO:BP (C5:GO:BP) uses its own mapping when selected. Versions and reference hashes saved per collection."),
    "GSEA exclusions: eight previously approved features B2M/ENSSSCG00000004687, LDHB/ENSSSCG00000000576, ARFIP2/TRIM3, GNA12/ENSSSCG00000032549; absent/filtered members need no exclusion.",
    paste("Additional approved GSEA exclusions by cluster:", paste(names(gsea_exclusions),
      vapply(gsea_exclusions, paste, character(1), collapse = ","), sep = "=", collapse = "; ")),
    "GSEA mapping: exact unambiguous pig symbol/Ensembl aliases, then existing safe pig name; NCBI canonical ID or Ensembl fallback; unmatched/ambiguous features remain SOURCE:<feature> nonhits. Additional tested-ID collisions block GSEA for that cluster pending owner decision.",
    "GSEA ranks: all tested features minus approved exclusions; sign(logFC)*sqrt(pmax(F,0)), decreasing with gene-ID tie-breaker. Numerical negative F values receive zero rank, matching edgeR P=1; original DE statistics unchanged. No gene FDR/fold-change filtering.",
    "GSEA: fgseaMultilevel(minSize=15,maxSize=500,gseaParam=1,scoreType=std,eps=0,nproc=1); seed1234 reset per contrast. BH separately per collection/cluster/contrast, not across collections/clusters/contrasts; pathway reporting FDR<0.01.",
    "GSEA direction: positive NES favours first phenotype; negative favours second. Hallmark/Reactome display up to eight terms/direction including nonsignificant terms, stars mark FDR<0.01. GO:BP dotplots/running scores show significant terms only; explicit no-significant-term output.",
    if ("gobp" %in% collections) "GO-only approved workaround (Nov4): ENSSSCG00000014800->NCBI:100525765/IL18BP; RNF121 retained separately. Additionally exclude ENSSSCG00000007481,ZFP64,ENSSSCG00000014098,ENSSSCG00000051837; original BD provenance unresolved. GO ranks differ only by these four extra exclusions and canonical-ID annotation, not source statistics. No GO depth/redundancy filtering.",
    "Interpretation: pressure is confounded with left/right lung, and Control is azygo lung at T0H. Baseline contrasts also combine time/perfusion effects; no claim of isolated pressure effects or pathway activation.",
    "Hallmark/Reactome retain historical ranks/exclusions; GO:BP has separately audited mapping and approved GO-only exclusions. Source RDS/references read-only; no automatic reference download or historical backfill.",
    "ORA: all genes with DEG BH FDR<=0.05 combined, no fold-change cutoff; tested unambiguous genes including unannotated nonhits as background. Duplicate canonical IDs are collapsed once per gene; any significant source feature selects that gene. 15–500 genes per pathway; hypergeometric upper tail; BH per collection/cluster/contrast; FDR<=0.05. GSEA-only exclusions are not automatically applied to ORA; ambiguous mappings and collapsed features are audited.",
    "Verification: one basic reopening of saved results and figure inspection after execution; no repeated audits or reruns.",
    paste("Evidence:", log_path),
    paste("R:", getRversion(), "; Seurat:", packageVersion("Seurat"), "; edgeR:", packageVersion("edgeR"),
      "; fgsea:", packageVersion("fgsea"), "; plotthis:", packageVersion("plotthis")))
  writeLines(info, info_path)
  finished <- FALSE
  on.exit(if (!finished) writeLines(c(info, "Status: failed; see execution log."), info_path))
  write.csv(sample_panel, file.path(output_root, "sample_panel.csv"), row.names = FALSE)
  write.csv(cohort_audit, file.path(output_root, "cohort_audit.csv"), row.names = FALSE)
  reference <- references <- NULL
  summaries <- eligibility <- list()
  for (cluster in as.character(clusters)) {
    directory <- file.path(output_root, sprintf("cluster_%02d", as.integer(cluster)))
    dir.create(file.path(directory, "gsea"), recursive = TRUE, showWarnings = FALSE)
    saveRDS(collections, file.path(directory, "gsea", "selected_collections.rds"))
    label <- paste(subset_label, "cluster", cluster)
    cluster_info <- c(paste("#", label), "Status: running", paste("Started:", Sys.time()),
      paste("Run settings:", info_path), paste("Input:", input), paste("Cluster column:", cluster_col),
      paste("Producer allocation:", Sys.getenv("SLURM_JOB_ID", unset = "local")),
      "Seed: 1234 for GSEA and label placement; DE/MDS/heatmap deterministic.", paste("Evidence:", log_path))
    writeLines(cluster_info, file.path(directory, ".INFO"))
    cat("\n", label, "started", format(Sys.time()), "\n")
    pb <- prepare_pseudobulk(object, common_genes, cluster_col, cluster, sample_panel,
      min_matched_pigs = min_matched_pigs)
    write.csv(pb$sample_info, file.path(directory, "sample_metadata.csv"), row.names = FALSE)
    write.csv(pb$sample_audit, file.path(directory, "sample_audit.csv"), row.names = FALSE)
    write.csv(pb$group_counts, file.path(directory, "group_counts.csv"), row.names = FALSE)
    write.csv(pb$contrast_status, file.path(directory, "contrast_status.csv"), row.names = FALSE)
    s <- data.frame(cluster, status = pb$status, reason = pb$reason,
      cells = sum(pb$sample_info$n_cells), observed_cells = sum(pb$sample_audit$n_cells),
      panel_samples = nrow(pb$sample_audit), samples = nrow(pb$sample_info),
      pigs = nlevels(pb$sample_info$pig),
      min_cells = if (nrow(pb$sample_info)) min(pb$sample_info$n_cells) else NA_integer_,
      excluded_samples = sum(!pb$sample_audit$retained),
      absent_samples = sum(pb$sample_audit$n_cells == 0L),
      below_10_samples = sum(pb$sample_audit$n_cells > 0L & pb$sample_audit$n_cells < 10L),
      candidate_genes = length(common_genes), design_rank = pb$design_rank, residual_df = pb$residual_df,
      filter_min_group_size = pb$filter_min_group_size, requested_contrasts = nrow(pb$contrast_status),
      eligible_contrasts = sum(pb$contrast_status$eligible), tested_contrasts = 0L)
    eligibility[[cluster]] <- s
    if (pb$status != "eligible") {
      summaries[[cluster]] <- data.frame(s, pb$contrast_status, tested = NA_integer_,
        up = NA_integer_, down = NA_integer_, heatmap_genes = NA_integer_,
        hallmark_status = if ("hallmark" %in% collections) paste("not run:", pb$status) else "not selected",
        reactome_status = if ("reactome" %in% collections) paste("not run:", pb$status) else "not selected",
        gobp_status = if ("gobp" %in% collections) paste("not run:", pb$status) else "not selected")
      write.csv(pb$contrast_status, file.path(directory, "DEG_summary.csv"), row.names = FALSE)
      cluster_info[cluster_info == "Status: running"] <- paste("Status:", pb$status, "; no fit, diagnostics, DE plots or GSEA.")
      writeLines(c(cluster_info, paste("Reason:", pb$reason), paste("Finished:", Sys.time()),
        paste("Retained/excluded specimens:", nrow(pb$sample_info), "/", sum(!pb$sample_audit$retained)),
        paste("Joint design rank:", pb$design_rank, "; residual df:", pb$residual_df),
        "Audits: sample_audit.csv, group_counts.csv and contrast_status.csv; no zero-expression libraries created."),
        file.path(directory, ".INFO"))
      cat(label, pb$status, ":", pb$reason, "\n")
      next
    }
    if (is.null(reference)) {
      reference <- prepare_hallmark(common_genes, mapping_file, reference_dir)
      references <- setNames(lapply(collections, function(collection) {
        if (collection == "hallmark") return(reference)
        selected <- prepare_pseudobulk_gsea(common_genes, mapping_file, reference_dir, collection)
        if (collection == "reactome") stopifnot(identical(reference$mapping, selected$mapping))
        selected
      }), collections)
      write.csv(reference$mapping, file.path(output_root, "candidate_gene_mapping.csv"), row.names = FALSE)
    }
    plot_pseudobulk_diagnostics(pb, file.path(directory, "diagnostics"), label)
    pb <- fit_pseudobulk(pb, directory)
    plots <- plot_pseudobulk_de(pb, file.path(directory, "plots"), label, reference$mapping)
    gsea <- setNames(lapply(collections, function(collection)
      run_pseudobulk_gsea(pb, references[[collection]], file.path(directory, "gsea", collection), label,
        excluded_features = gsea_exclusions[[cluster]])), collections)
    ora <- setNames(lapply(collections, function(collection)
      run_pseudobulk_ora(pb, references[[collection]], file.path(directory, "ora", collection), label)), collections)
    s$status <- pb$status
    s$tested_contrasts <- length(pb$results)
    eligibility[[cluster]] <- s
    stats <- data.frame(pb$summary,
      up_fdr_logfc = vapply(pb$results, function(x) sum(x$FDR <= 0.05 & x$logFC >= 0.5), integer(1)),
      down_fdr_logfc = vapply(pb$results, function(x) sum(x$FDR <= 0.05 & x$logFC <= -0.5), integer(1)))
    contrast_summary <- merge(pb$contrast_status, stats, by = "contrast", all.x = TRUE, sort = FALSE)
    for (collection in specs$collection) {
      contrast_summary[[paste0(collection, "_status")]] <- if (!collection %in% collections) "not selected" else
        ifelse(contrast_summary$contrast_status == "tested", gsea[[collection]]$status, "not run: skipped contrast")
      contrast_summary[[paste0(collection, "_ora_status")]] <- if (!collection %in% collections) "not selected" else
        ifelse(contrast_summary$contrast_status == "tested", ora[[collection]]$status, "not run: skipped contrast")
    }
    s <- data.frame(s, heatmap_genes = plots$n_heatmap_genes, contrast_summary, row.names = NULL)
    for (collection in collections) {
      coverage <- gsea[[collection]]$summary
      if (is.null(coverage)) next
      # Historical unprefixed coverage columns remain Hallmark-only.
      if (collection != "hallmark") names(coverage)[names(coverage) != "contrast"] <-
        paste0(collection, "_", names(coverage)[names(coverage) != "contrast"])
      s <- merge(s, coverage, by = "contrast", all.x = TRUE, sort = FALSE)
    }
    s <- s[match(pb$contrast_status$contrast, s$contrast), , drop = FALSE]
    summaries[[cluster]] <- s
    collection_status <- paste(collections, "GSEA", vapply(gsea, function(x) x$status, character(1)),
      "ORA", vapply(ora, function(x) x$status, character(1)), collapse = "; ")
    cluster_info[cluster_info == "Status: running"] <- paste("Status: DE/diagnostics/plots completed;", collection_status)
    writeLines(c(cluster_info, paste("Finished:", Sys.time()),
      paste("Cells:", sum(pb$sample_info$n_cells), "; samples:", ncol(pb$counts),
        "; pigs:", nlevels(pb$sample_info$pig), "; min cells per sample:", min(pb$sample_info$n_cells)),
      paste("Candidate genes:", length(common_genes), "; low-expression genes removed:", sum(!pb$keep),
        "; tested:", sum(pb$keep), "; genes outside original common-layer intersection:", nrow(object) - length(common_genes)),
      paste("Excluded specimens:", sum(!pb$sample_audit$retained), "; requested/tested contrasts:",
        nrow(pb$contrast_status), "/", length(pb$results)),
      paste("Joint design rank:", pb$design_rank, "; residual df:", pb$residual_df,
        "; filterByExpr minimum retained group size:", pb$filter_min_group_size),
      paste("Heatmap genes:", plots$n_heatmap_genes),
      "Validation: one saved-output review and figure inspection pending."), file.path(directory, ".INFO"))
    cat(label, "completed;", collection_status, "\n")
    print(s, row.names = FALSE)
  }
  # Some clusters may require an additional mapping decision before GSEA.
  summary <- dplyr::bind_rows(summaries)
  write.csv(summary, file.path(output_root, "cluster_summary.csv"), row.names = FALSE)
  write.csv(dplyr::bind_rows(eligibility), file.path(output_root, "cluster_eligibility.csv"), row.names = FALSE)
  completed <- summary$contrast_status == "tested"
  info[info == "Status: running"] <- if (!any(completed)) {
    "Status: no contrasts tested; see cluster statuses and requested-contrast skip reasons. No models or figures run."
  } else if (all(vapply(collections, function(collection)
    all(summary[[paste0(collection, "_status")]][completed] == "completed") &&
      all(summary[[paste0(collection, "_ora_status")]][completed] == "completed"), logical(1)))) {
    "Status: supported contrasts completed; saved-output and figure review pending. See skipped-contrast/blocked-cluster audits."
  } else {
    "Status: supported contrasts DE/plots completed; some enrichment analyses blocked by identifier mapping. See contrast/cluster audits."
  }
  writeLines(c(info, paste("Finished:", Sys.time()), paste("Original count layers:", length(count_layers),
    "; candidate genes:", length(common_genes)), "Summary: cluster_summary.csv"), info_path)
  for (contrast in unique(summary$contrast[completed])) {
    plot_DEG_for_pub(output_root, contrast = contrast)
  }
  publication_figures <- collect_deg_publication_figures(output_root)
  info <- c(info, paste("Publication figure copies:", length(publication_figures),
    "; directory:", file.path(dirname(output_root), "publications_figures", basename(output_root)),
    "; copies verified against source figures by MD5."))
  writeLines(c(info, paste("Finished:", Sys.time()), "Summary: cluster_summary.csv"), info_path)
  finished <- TRUE
  invisible(summary)
}

if (sys.nframe() == 0L) {
  source("R/DEG/cli.R")
  run_deg_cli(commandArgs(trailingOnly = TRUE))
}
