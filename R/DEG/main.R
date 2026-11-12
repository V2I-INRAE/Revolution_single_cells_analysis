# Paired pseudobulk DEG and GSEA orchestrator.
# CLI: Rscript R/DEG/main.R --help. Inputs and new output root are explicit;
# frozen 2026.1 pig pathway references are shared and read-only.
suppressPackageStartupMessages({
  library(Seurat)
  library(edgeR)
  library(ggplot2)
})
source("R/utils/save_plots.R")
source("R/DEG/pseudobulk_de.R")
source("R/utils/pig_gene_mapping.R")
source("R/DEG/pseudobulk_plots.R")
source("R/DEG/pseudobulk_gsea.R")
source("R/plots/publications/plot_DEG_for_pub.R")
source("R/plots/publications/publication_paths.R")

run_pseudobulk <- function(input, subset_label, cluster_col, clusters, output_root,
    mapping_file = "results/misc/Ensembl_ID/gene-mapping/pig_to_human_orthologues_annotated.csv",
    reference_dir = "results/misc/msigdb_pig_2026.1", gsea_exclusions = list(),
    pigs = NULL, min_matched_pigs = 3L, log_path = NULL,
    collections = c("hallmark", "reactome", "gobp")) {
  specs <- pseudobulk_gsea_collections()
  if (!length(collections) || anyNA(collections) || anyDuplicated(collections) ||
      !all(collections %in% specs$collection)) stop("Select unique collections from hallmark, reactome, gobp.")
  if (length(min_matched_pigs) != 1L || !is.finite(min_matched_pigs) ||
      min_matched_pigs < 3L || min_matched_pigs %% 1 != 0) {
    stop("min_matched_pigs must be an integer >=3; the complete-pair rule cannot be relaxed.")
  }
  if (!is.null(pigs) && (!length(pigs) || anyNA(pigs) || any(!nzchar(pigs)) || anyDuplicated(pigs))) {
    stop("Supply a nonempty, unique pig list, or NULL for all pigs.")
  }
  if (missing(output_root) || length(output_root) != 1L || is.na(output_root) || !nzchar(output_root)) {
    stop("Supply the resolution-scoped DEG output_root; existing analysis outputs must not be overwritten.")
  }
  output_root <- path.expand(output_root)
  if (file.exists(output_root)) stop("DEG results already exist; refusing to overwrite or create a duplicate run directory.")
  resolution <- sub("^harmony_snn_res[.]", "", cluster_col)
  resolution_dir <- subset_resolution_dir(input, resolution)
  expected_output <- file.path(resolution_dir, "DEG")
  if (fs::path_abs(output_root) != fs::path_abs(expected_output)) {
    stop("Use the single resolution-scoped DEG directory: ", expected_output)
  }
  publication_dir <- subset_publication_dir(input, resolution)
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
  if (identical(clusters, "all")) clusters <- as.character(sort(unique(as.integer(as.character(metadata[[cluster_col]])))))
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
    "Script: R/DEG/main.R; count/DE, annotated display tables, publication DEG plots, GSEA. No redundant diagnostic/volcano/heatmap image exports.",
    paste("Selected GSEA collections:", paste(collections, collapse = ", ")),
    paste("Input (read-only):", input), paste("Source producer:", basename(dirname(input))),
    paste("Population:", subset_label, "; cluster column:", cluster_col, "; clusters:", paste(clusters, collapse = ", ")),
    "Groups (ordered): Control/T0H, P4/T4H positive, N4/T4H negative, P10/T10H positive, N10/T10H negative airway-pressure ventilation.",
    paste("Cohort:", paste(levels(sample_panel$pig), collapse = ", "), "; selected/source specimens:",
      nrow(sample_panel), "/", nrow(cohort_audit), "; Control/P4/N4/P10/N10:",
      paste(as.integer(table(sample_panel$condition)), collapse = "/")),
    "Cohort exclusions: cohort_audit.csv records every source specimen and pig_not_selected reasons.",
    "Model scope: separate comparison-specific complete pig pairs; not a joint all-condition model.",
    "T8H becomes T10H character-safely in memory only; time_point_original and original sample IDs (including REVO31-N8) retained for audit.",
    "Exclude individual specimen/population aggregates with <25 cells (including absent populations). Then exclude specimens outside the comparison and pigs without both compared conditions. Audit all original specimens; no zero-expression libraries or imputation.",
    paste0("Test each requested contrast only with >=", min_matched_pigs,
      " COMPLETE MATCHED PIG PAIRS after the 25-cell filter; each retained pig contributes exactly one specimen in each compared condition. All six requests/counts/pairs/skip reasons are recorded."),
    "No eligible contrasts => no fit. Rank-deficient design or residual df<=0 => blocked, no fit; no model simplification or further sample removal.",
    "Candidates: intersection of ALL original RNA counts.* layers, before sample/time/cluster selection.",
    "Filtering independently per cluster/contrast: filterByExpr(group=condition), min.count=10/min.total.count=15/large.n=10/min.prop=0.7, using only complete pairs. No 20%-of-cells prevalence filter; minimum retained group size equals paired pigs (>=3).",
    "Normalization: recalculate filtered library sizes, then TMM. No integrated or normalized cell counts used.",
    "Model: ~ pig + condition; edgeR glmQLFit(legacy=FALSE), dispersion=NULL, robust=FALSE; glmQLFTest.",
    "Requested contrasts (ordered, first-minus-second; intercept/pig coefficients zero): P4_vs_Control, N4_vs_Control, P4_vs_N4, P4_vs_P10, N4_vs_N10, P10_vs_N10. One separate paired fit per eligible cluster/contrast; full rank and positive residual df required. BH separately per cluster/contrast.",
    "DE exports: full statistics and FDR<=0.05 tables, without an effect-size cutoff.",
    "Publication DEG display selection: FDR<=0.05 AND abs(log2FC)>=0.5 within the same tested contrast.",
    "Diagnostic figures: not automatically exported; saved fits and sample/model audits retained.",
    "No pooled heatmap matrix: sample sets and tested genes are comparison-specific; no heatmap/diagnostic image exports.",
    "Volcano: exact logFC and -log10(BH FDR), four classes, +/-0.5 and FDR0.05 dashed guides; five labels per direction.",
    "Pig display-name mapping: R/utils/pig_gene_mapping.R; unflagged pig Ensembl annotation and unambiguous pig-side frozen MSigDB rescue. Feature keys retained; duplicate symbols ID-qualified; review-required symbols unresolved. GSEA canonical mapping remains independent.",
    paste("References:", reference_dir, "; frozen MSigDB pig projections, taxon9823. Historical H/Reactome canonical mapping and independent GO:BP mapping retained. Display annotation is separate. Versions and reference hashes saved per collection."),
    "GSEA exclusions: eight previously approved features B2M/ENSSSCG00000004687, LDHB/ENSSSCG00000000576, ARFIP2/TRIM3, GNA12/ENSSSCG00000032549; absent/filtered members need no exclusion.",
    paste("Additional approved GSEA exclusions by cluster or cluster/contrast:", paste(names(gsea_exclusions),
      vapply(gsea_exclusions, paste, character(1), collapse = ","), sep = "=", collapse = "; ")),
    "GSEA mapping: exact unambiguous pig symbol/Ensembl aliases, then existing safe pig name; NCBI canonical ID or Ensembl fallback; unmatched/ambiguous features remain SOURCE:<feature> nonhits. Additional tested-ID collisions block GSEA for that cluster pending owner decision.",
    "GSEA ranks: all tested features minus approved exclusions; sign(logFC)*sqrt(pmax(F,0)), decreasing with gene-ID tie-breaker. Numerical negative F values receive zero rank, matching edgeR P=1; original DE statistics unchanged. No gene FDR/fold-change filtering.",
    "GSEA: fgseaMultilevel(minSize=15,maxSize=500,gseaParam=1,scoreType=std,eps=0,nproc=1); seed1234 reset per contrast. BH separately per collection/cluster/contrast, not across collections/clusters/contrasts; pathway reporting FDR<0.01.",
    "Uniform GSEA retry rule for every collection/contrast: initial nPermSimple=1000; if any estimates unavailable, retry the whole contrast once at10000 with seed1234; preserve initial unavailable terms and retry_audit.csv; stop if still unavailable. No collection-specific overrides.",
    "GSEA direction: positive NES favours first phenotype; negative favours second. Hallmark/Reactome display up to eight terms/direction including nonsignificant terms, stars mark FDR<0.01. GO:BP dotplots/running scores show significant terms only; explicit no-significant-term output.",
    if ("gobp" %in% collections) "GO-only approved workaround (Nov4): ENSSSCG00000014800->NCBI:100525765/IL18BP; RNF121 retained separately. Additionally exclude ENSSSCG00000007481,ZFP64,ENSSSCG00000014098,ENSSSCG00000051837; original BD provenance unresolved. GO ranks differ only by these four extra exclusions and canonical-ID annotation, not source statistics. No GO depth/redundancy filtering.",
    "Interpretation: pressure is confounded with left/right lung, and Control is azygo lung at T0H. Baseline contrasts also combine time/perfusion effects; no claim of isolated pressure effects or pathway activation.",
    "Canonical mapping and approved collection/cluster exclusions retained; ranks are recalculated from the new paired fits. Source objects/references read-only; no automatic reference download or historical backfill.",
    "Verification: one basic reopening of saved results and figure inspection after execution; no repeated audits or reruns.",
    paste("Evidence:", log_path),
    paste("R:", getRversion(), "; Seurat:", packageVersion("Seurat"), "; edgeR:", packageVersion("edgeR"),
      "; fgsea:", packageVersion("fgsea"), "; plotthis:", packageVersion("plotthis")))
  writeLines(info, info_path)
  finished <- FALSE
  on.exit(if (!finished) writeLines(c(info, "Status: failed; see execution log."), info_path))
  write.csv(sample_panel, file.path(output_root, "sample_panel.csv"), row.names = FALSE)
  write.csv(cohort_audit, file.path(output_root, "cohort_audit.csv"), row.names = FALSE)
  gene_mapping <- prepare_pig_gene_mapping(rownames(object[["RNA"]]), reference_dir = reference_dir)
  write.csv(gene_mapping, file.path(output_root, "gene_display_mapping.csv"), row.names = FALSE)
  reference <- references <- NULL
  summaries <- eligibility <- list()
  for (cluster in as.character(clusters)) {
    directory <- file.path(output_root, sprintf("cluster_%02d", as.integer(cluster)))
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    label <- paste(subset_label, "cluster", cluster)
    cluster_info <- c(paste("#", label), "Status: running", paste("Started:", Sys.time()),
      paste("Run settings:", info_path), paste("Input:", input), paste("Cluster column:", cluster_col),
      "Each comparison below has its own complete-pair model, filtered gene universe and normalization.",
      "population_sample_audit.csv records the 25-cell gate only; comparison sample_audit.csv records actual pair retention.",
      paste("Evidence:", log_path))
    writeLines(cluster_info, file.path(directory, ".INFO"))
    cat("\n", label, "started", format(Sys.time()), "\n")
    population <- prepare_pseudobulk(object, common_genes, cluster_col, cluster, sample_panel,
      min_matched_pigs = min_matched_pigs)
    write.csv(population$sample_audit, file.path(directory, "population_sample_audit.csv"), row.names = FALSE)
    rows <- statuses <- plot_results <- list()
    for (contrast in pseudobulk_comparisons()$contrast) {
      comparison_dir <- file.path(directory, contrast)
      dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)
      comparison_info <- c(paste("#", label, contrast), "Status: running", paste("Started:", Sys.time()),
        paste("Run settings:", info_path), paste("Input:", input), paste("Cluster column:", cluster_col),
        "Method: >=25 cells/specimen, >=3 complete pig pairs; separate filterByExpr/TMM/~pig+condition edgeR fit.",
        "Seed: 1234 for GSEA/labels; deterministic DE.", paste("Evidence:", log_path))
      writeLines(comparison_info, file.path(comparison_dir, ".INFO"))
      pb <- prepare_paired_comparison(population, contrast)
      write.csv(pb$sample_audit, file.path(comparison_dir, "sample_audit.csv"), row.names = FALSE)
      write.csv(pb$sample_info, file.path(comparison_dir, "sample_metadata.csv"), row.names = FALSE)
      write.csv(pb$contrast_status, file.path(comparison_dir, "contrast_status.csv"), row.names = FALSE)
      s <- data.frame(cluster, status = pb$status, reason = pb$reason,
        cells = sum(pb$sample_info$n_cells), samples = nrow(pb$sample_info), pigs = nlevels(pb$sample_info$pig),
        min_retained_cells = if (nrow(pb$sample_info)) min(pb$sample_info$n_cells) else NA_integer_,
        below_25_samples = sum(pb$sample_audit$in_comparison & !pb$sample_audit$passes_cell_filter),
        unmatched_samples = sum(pb$sample_audit$exclusion_reason == "unmatched_pig"),
        candidate_genes = length(common_genes), design_rank = pb$design_rank, residual_df = pb$residual_df,
        filter_min_group_size = pb$filter_min_group_size, pb$contrast_status,
        tested = NA_integer_, up = NA_integer_, down = NA_integer_,
        up_fdr_logfc = NA_integer_, down_fdr_logfc = NA_integer_, row.names = NULL)
      for (collection in specs$collection) s[[paste0(collection, "_status")]] <-
        if (collection %in% collections) "not run: skipped contrast" else "not selected"
      if (pb$status == "eligible") {
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
        pb <- fit_pseudobulk(pb, comparison_dir, gene_mapping)
        plot_results[[contrast]] <- pb$results[[contrast]]
        gsea_dir <- file.path(comparison_dir, "gsea")
        dir.create(gsea_dir, showWarnings = FALSE)
        saveRDS(collections, file.path(gsea_dir, "selected_collections.rds"))
        gsea <- setNames(lapply(collections, function(collection)
          run_pseudobulk_gsea(pb, references[[collection]], file.path(gsea_dir, collection), label,
            excluded_features = unique(c(gsea_exclusions[[cluster]],
              gsea_exclusions[[paste(cluster, contrast, sep = "/")]])))), collections)
        s$status <- "completed"
        s$contrast_status <- "tested"
        s$tested <- pb$summary$tested
        s$up <- pb$summary$up
        s$down <- pb$summary$down
        tab <- pb$results[[contrast]]
        s$up_fdr_logfc <- sum(tab$FDR <= 0.05 & tab$logFC >= 0.5)
        s$down_fdr_logfc <- sum(tab$FDR <= 0.05 & tab$logFC <= -0.5)
        for (collection in collections) {
          s[[paste0(collection, "_status")]] <- gsea[[collection]]$status
          coverage <- gsea[[collection]]$summary
          if (is.null(coverage)) next
          coverage$contrast <- NULL
          names(coverage) <- paste0(collection, "_", names(coverage))
          s <- cbind(s, coverage)
        }
        comparison_info[comparison_info == "Status: running"] <- paste("Status: DE completed;",
          paste(collections, vapply(gsea, function(x) x$status, character(1)), collapse = "; "))
        comparison_info <- c(comparison_info, paste("Tested genes:", sum(pb$keep), "; filtered:", sum(!pb$keep)),
          paste("Design rank:", pb$design_rank, "; residual df:", pb$residual_df))
      } else {
        write.csv(pb$contrast_status, file.path(comparison_dir, "DEG_summary.csv"), row.names = FALSE)
        comparison_info[comparison_info == "Status: running"] <- paste("Status: skipped;", pb$reason, "; no model or GSEA run.")
      }
      writeLines(c(comparison_info, paste("Matched pigs:", pb$contrast_status$matched_pig_ids),
        paste("Retained specimens:", nrow(pb$sample_info), "; minimum cells:", s$min_retained_cells),
        paste("Finished:", Sys.time()), "Validation: saved-output and figure review pending."),
        file.path(comparison_dir, ".INFO"))
      rows[[contrast]] <- s
      statuses[[contrast]] <- pb$contrast_status
      cat(label, contrast, s$status, "; pairs", s$matched_pigs, "; genes", s$tested, "\n")
    }
    if (length(plot_results)) prepare_pseudobulk_plot_data(list(results = plot_results),
      file.path(directory, "plots"), gene_mapping)
    cluster_summary <- dplyr::bind_rows(rows)
    write.csv(dplyr::bind_rows(statuses), file.path(directory, "contrast_status.csv"), row.names = FALSE)
    write.csv(cluster_summary, file.path(directory, "DEG_summary.csv"), row.names = FALSE)
    summaries[[cluster]] <- cluster_summary
    eligibility[[cluster]] <- data.frame(cluster, requested = nrow(cluster_summary),
      tested = sum(cluster_summary$contrast_status == "tested"),
      skipped = sum(cluster_summary$contrast_status == "skipped"),
      observed_cells = sum(population$sample_audit$n_cells),
      samples_passing_25 = sum(population$sample_audit$passes_cell_filter))
    cluster_info[cluster_info == "Status: running"] <- paste("Status: comparisons processed; tested",
      length(plot_results), "/", nrow(cluster_summary), "; see per-comparison enrichment statuses.")
    writeLines(c(cluster_info, paste("Finished:", Sys.time()),
      "Summary: DEG_summary.csv; comparison-specific models and audits in <contrast>/; no shared joint fit.",
      "Validation: saved-output and figure review pending."), file.path(directory, ".INFO"))
  }
  summary <- dplyr::bind_rows(summaries)
  write.csv(summary, file.path(output_root, "cluster_summary.csv"), row.names = FALSE)
  write.csv(dplyr::bind_rows(eligibility), file.path(output_root, "cluster_eligibility.csv"), row.names = FALSE)
  completed <- summary$contrast_status == "tested"
  info[info == "Status: running"] <- if (!any(completed)) {
    "Status: no comparisons eligible; explicit skip audits saved; no models or figures run."
  } else if (all(vapply(collections, function(collection)
    all(summary[[paste0(collection, "_status")]][completed] == "completed"), logical(1)))) {
    "Status: supported comparisons completed; saved-output and figure review pending."
  } else "Status: DE completed; some enrichment blocked pending explicit identifier decisions."
  writeLines(c(info, paste("Finished:", Sys.time()), paste("Original count layers:", length(count_layers),
    "; candidate genes:", length(common_genes)), "Summary: cluster_summary.csv"), info_path)
  for (contrast in pseudobulk_comparisons()$contrast) {
    plot_DEG_for_pub(output_root, output_dir = publication_dir, contrast = contrast)
  }
  finished <- TRUE
  invisible(summary)
}

if (sys.nframe() == 0L) {
  source("R/DEG/cli.R")
  run_deg_cli(commandArgs(trailingOnly = TRUE))
}
