# From the project root, source this file and call run_pseudobulk() explicitly.
# output_root must be a NEW directory, never an existing historical DEG directory.
# Supply reference_dir pointing to the existing frozen Hallmark/Reactome inputs;
# the default references location is not populated or downloaded by this driver.
# Example (choose an unused output path and supply the true frozen references):
# run_pseudobulk(
#   input = "data/job-44542861/subsets/Endothelial/job-44831031/subclustered.rds",
#   subset_label = "Endothelial", cluster_col = "harmony_snn_res.0.3", clusters = 1:5,
#   output_root = "results/T10H_Endothelial_joint5group_new",
#   reference_dir = "/path/to/existing/frozen/references",
#   gsea_exclusions = list("2" = c("TNFSF15", "ENSSSCG00000034087"),
#     "5" = c("TNFSF15", "ENSSSCG00000034087")))
suppressPackageStartupMessages({
  library(Seurat)
  library(edgeR)
  library(ggplot2)
})
source("R/utils/save_plots.R")
source("R/DEG/pseudobulk_de.R")
source("R/DEG/pseudobulk_diagnostics.R")
source("R/DEG/pseudobulk_plots.R")
source("R/DEG/pseudobulk_hallmark.R")

run_pseudobulk <- function(input, subset_label, cluster_col, clusters, output_root,
    mapping_file = "results/misc/Ensembl_ID/gene-mapping/pig_to_human_orthologues_annotated.csv",
    reference_dir = file.path(output_root, "references"), gsea_exclusions = list()) {
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
  if (nrow(sample_panel) != 23L || !identical(as.integer(table(sample_panel$condition)), c(5L, 5L, 5L, 4L, 4L)) ||
      !identical(levels(sample_panel$pig), c("REVO26", "REVO27", "REVO29", "REVO30", "REVO31")) ||
      any(sample_panel$pig == "REVO29" & sample_panel$condition %in% c("P10", "N10")) ||
      !identical(sample_panel$sample[sample_panel$pig == "REVO31" & sample_panel$condition == "N10"], "REVO31-N8")) {
    stop("Observed panel must contain the original 23 specimens (5/5/5/4/4), no REVO29 late specimens, and REVO31-N8 unchanged.")
  }
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
  if (!dir.create(output_root, recursive = TRUE, showWarnings = FALSE)) stop("Could not create a new output_root.")
  info_path <- file.path(output_root, ".INFO")
  log_path <- paste0("logs/", tolower(subset_label), "-all-clusters-pseudobulk-",
    format(Sys.Date(), "%Y%m%d"), ".log")
  info <- c(paste("#", subset_label, "paired pseudobulk pipeline"), "Status: running",
    paste("Started:", Sys.time()), paste("Producer allocation:", Sys.getenv("SLURM_JOB_ID", unset = "local")),
    "Script: R/DEG/run_pseudobulk.R; count/DE, diagnostics, DE plotting and Hallmark modules.",
    paste("Input (read-only):", input), paste("Source producer:", basename(dirname(input))),
    paste("Population:", subset_label, "; cluster column:", cluster_col, "; clusters:", paste(clusters, collapse = ", ")),
    "Groups (ordered): Control/T0H, P4/T4H positive, N4/T4H negative, P10/T10H positive, N10/T10H negative airway-pressure ventilation.",
    "Observed panel before cluster selection: 23 original specimens; Control/P4/N4=5 each, P10/N10=4 each; REVO29 has no late specimens.",
    "T8H becomes T10H character-safely in memory only; time_point_original and original sample IDs (including REVO31-N8) retained for audit.",
    "Exclude only individual observed-panel aggregates with <10 cells (including absent cluster populations); sample_audit.csv retains exact original IDs, n_cells and exclusion reasons. No zero-expression libraries.",
    "Retain ALL passing specimens in one joint model, including late groups with only 1 or 2 replicates. No imputation, additional deletion or downsampling.",
    "Test each requested contrast only with >=3 DIRECTLY MATCHED PIGS after sample filtering; missing groups or fewer pairs are explicitly skipped. All six requests/counts/pairs/reasons appear in contrast_status.csv and run summaries.",
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
    paste("References:", reference_dir, "; MSigDB2026.1.Hs Hallmark pig projection, taxon9823. Reactome used only for consistent pig aliases, not enrichment tests."),
    "GSEA exclusions: eight previously approved features B2M/ENSSSCG00000004687, LDHB/ENSSSCG00000000576, ARFIP2/TRIM3, GNA12/ENSSSCG00000032549; absent/filtered members need no exclusion.",
    paste("Additional approved GSEA exclusions by cluster:", paste(names(gsea_exclusions),
      vapply(gsea_exclusions, paste, character(1), collapse = ","), sep = "=", collapse = "; ")),
    "GSEA mapping: exact unambiguous pig symbol/Ensembl aliases, then existing safe pig name; NCBI canonical ID or Ensembl fallback; unmatched/ambiguous features remain SOURCE:<feature> nonhits. Additional tested-ID collisions block GSEA for that cluster pending owner decision.",
    "GSEA ranks: all tested features minus approved exclusions; sign(logFC)*sqrt(pmax(F,0)), decreasing with gene-ID tie-breaker. Numerical negative F values receive zero rank, matching edgeR P=1; original DE statistics unchanged. No gene FDR/fold-change filtering.",
    "GSEA: fgseaMultilevel(minSize=15,maxSize=500,gseaParam=1,scoreType=std,eps=0,nproc=1); seed1234 reset per contrast. BH separately per cluster/contrast; pathway reporting FDR<0.01.",
    "GSEA direction: positive NES favours first phenotype; negative favours second. Display up to eight terms/direction, including nonsignificant terms; stars mark FDR<0.01. Lowest-FDR running score in each direction.",
    "Interpretation: baseline contrasts combine time/perfusion and pressure effects; no claim of exclusive induction or pathway activation.",
    "No ORA/GO/Reactome analysis in this run. Source RDS and input references are read-only.",
    "Verification: one basic reopening of saved results and figure inspection after execution; no repeated audits or reruns.",
    paste("Evidence:", log_path),
    paste("R:", getRversion(), "; Seurat:", packageVersion("Seurat"), "; edgeR:", packageVersion("edgeR"),
      "; fgsea:", packageVersion("fgsea"), "; plotthis:", packageVersion("plotthis")))
  writeLines(info, info_path)
  finished <- FALSE
  on.exit(if (!finished) writeLines(c(info, "Status: failed; see execution log."), info_path))
  write.csv(sample_panel, file.path(output_root, "sample_panel.csv"), row.names = FALSE)
  reference <- NULL
  summaries <- eligibility <- list()
  for (cluster in as.character(clusters)) {
    directory <- file.path(output_root, sprintf("cluster_%02d", as.integer(cluster)))
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    label <- paste(subset_label, "cluster", cluster)
    cluster_info <- c(paste("#", label), "Status: running", paste("Started:", Sys.time()),
      paste("Run settings:", info_path), paste("Input:", input), paste("Cluster column:", cluster_col),
      paste("Producer allocation:", Sys.getenv("SLURM_JOB_ID", unset = "local")),
      "Seed: 1234 for GSEA and label placement; DE/MDS/heatmap deterministic.", paste("Evidence:", log_path))
    writeLines(cluster_info, file.path(directory, ".INFO"))
    cat("\n", label, "started", format(Sys.time()), "\n")
    pb <- prepare_pseudobulk(object, common_genes, cluster_col, cluster, sample_panel)
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
        hallmark_status = paste("not run:", pb$status))
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
      write.csv(reference$mapping, file.path(output_root, "candidate_gene_mapping.csv"), row.names = FALSE)
    }
    plot_pseudobulk_diagnostics(pb, file.path(directory, "diagnostics"), label)
    pb <- fit_pseudobulk(pb, directory)
    plots <- plot_pseudobulk_de(pb, file.path(directory, "plots"), label, reference$mapping)
    gsea <- run_pseudobulk_hallmark(pb, reference, file.path(directory, "gsea", "hallmark"), label,
      excluded_features = gsea_exclusions[[cluster]])
    s$status <- pb$status
    s$tested_contrasts <- length(pb$results)
    eligibility[[cluster]] <- s
    stats <- data.frame(pb$summary,
      up_fdr_logfc = vapply(pb$results, function(x) sum(x$FDR <= 0.05 & x$logFC >= 0.5), integer(1)),
      down_fdr_logfc = vapply(pb$results, function(x) sum(x$FDR <= 0.05 & x$logFC <= -0.5), integer(1)))
    contrast_summary <- merge(pb$contrast_status, stats, by = "contrast", all.x = TRUE, sort = FALSE)
    contrast_summary$hallmark_status <- ifelse(contrast_summary$contrast_status == "tested",
      gsea$status, "not run: skipped contrast")
    s <- data.frame(s, heatmap_genes = plots$n_heatmap_genes, contrast_summary, row.names = NULL)
    if (!is.null(gsea$summary)) s <- merge(s, gsea$summary, by = "contrast", all.x = TRUE, sort = FALSE)
    s <- s[match(pb$contrast_status$contrast, s$contrast), , drop = FALSE]
    summaries[[cluster]] <- s
    cluster_info[cluster_info == "Status: running"] <- paste("Status: DE/diagnostics/plots completed; Hallmark", gsea$status)
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
    cat(label, "completed; Hallmark", gsea$status, "\n")
    print(s, row.names = FALSE)
  }
  # Some clusters may require an additional mapping decision before GSEA.
  summary <- dplyr::bind_rows(summaries)
  write.csv(summary, file.path(output_root, "cluster_summary.csv"), row.names = FALSE)
  write.csv(dplyr::bind_rows(eligibility), file.path(output_root, "cluster_eligibility.csv"), row.names = FALSE)
  completed <- summary$contrast_status == "tested"
  info[info == "Status: running"] <- if (!any(completed)) {
    "Status: no contrasts tested; see cluster statuses and requested-contrast skip reasons. No models or figures run."
  } else if (all(summary$hallmark_status[completed] == "completed")) {
    "Status: supported contrasts completed; saved-output and figure review pending. See skipped-contrast/blocked-cluster audits."
  } else {
    "Status: supported contrasts DE/plots completed; some Hallmark analyses blocked by identifier collisions. See contrast/cluster audits."
  }
  writeLines(c(info, paste("Finished:", Sys.time()), paste("Original count layers:", length(count_layers),
    "; candidate genes:", length(common_genes)), "Summary: cluster_summary.csv"), info_path)
  finished <- TRUE
  invisible(summary)
}

if (sys.nframe() == 0L) stop(paste(
  "Source R/DEG/run_pseudobulk.R and call run_pseudobulk() with an explicit NEW output_root,",
  "the original input/cluster settings, existing reference_dir, and approved gsea_exclusions.",
  "No default production destination or automatic analysis is permitted."))
