# Run from the project root: Rscript R/DEG/run_pseudobulk.R
# Or source this file and call run_pseudobulk() for another subset or cluster list.
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
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  info_path <- file.path(output_root, ".INFO")
  log_path <- paste0("logs/", tolower(subset_label), "-all-clusters-pseudobulk-",
    format(Sys.Date(), "%Y%m%d"), ".log")
  info <- c(paste("#", subset_label, "paired pseudobulk pipeline"), "Status: running",
    paste("Started:", Sys.time()), paste("Producer allocation:", Sys.getenv("SLURM_JOB_ID", unset = "local")),
    "Script: R/DEG/run_pseudobulk.R; count/DE, diagnostics, DE plotting and Hallmark modules.",
    paste("Input (read-only):", input), paste("Source producer:", basename(dirname(input))),
    paste("Population:", subset_label, "; cluster column:", cluster_col, "; clusters:", paste(clusters, collapse = ", ")),
    "Groups: paired Control/T0H, P4/T4H positive airway-pressure ventilation, N4/T4H negative-pressure ventilation.",
    "T8H becomes T10H in memory only. Late samples are outside this analysis panel.",
    "All panel samples and cells retained within each requested cluster; no downsampling. Per-sample cell counts are recorded in sample_metadata.csv.",
    "Candidates: intersection of ALL original RNA counts.* layers, before sample/time/cluster selection.",
    "Filtering per cluster: filterByExpr(group=condition), default min.count=10/min.total.count=15/large.n=10/min.prop=0.7.",
    "Normalization: recalculate filtered library sizes, then TMM. No integrated or normalized cell counts used.",
    "Model: ~ pig + condition; edgeR glmQLFit(legacy=FALSE), dispersion=NULL, robust=FALSE; glmQLFTest.",
    "Contrasts: P4-Control, N4-Control, P4-N4. BH FDR separately per cluster/contrast; both gene directions retained.",
    "DE exports: full statistics and FDR<=0.05 tables, without an effect-size cutoff.",
    "Volcano/heatmap selection: FDR<=0.05 AND abs(log2FC)>=0.5 within the same contrast; heatmap union across contrasts.",
    "MDS: same filtered TMM DGEList as DE; top500 pairwise leading-logFC distances, prior.count=2; no pig-effect removal.",
    "Heatmap: TMM log2 CPM prior.count=2, row z-scores over all samples, Euclidean/complete row clustering, fixed Control/P4/N4 then pig order; RdBu blue-low/red-high saturated at +/-2.5; pig effects retained.",
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
  object <- readRDS(input)
  object$time_point[object$time_point == "T8H"] <- "T10H"
  count_layers <- SeuratObject::Layers(object[["RNA"]], search = "^counts\\.")
  common_genes <- Reduce(intersect, lapply(count_layers, function(layer)
    SeuratObject::Features(object[["RNA"]], layer = layer)))
  reference <- prepare_hallmark(common_genes, mapping_file, reference_dir)
  write.csv(reference$mapping, file.path(output_root, "candidate_gene_mapping.csv"), row.names = FALSE)
  summaries <- list()
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
    pb <- prepare_pseudobulk(object, common_genes, cluster_col, cluster)
    plot_pseudobulk_diagnostics(pb, file.path(directory, "diagnostics"), label)
    pb <- fit_pseudobulk(pb, directory)
    plots <- plot_pseudobulk_de(pb, file.path(directory, "plots"), label, reference$mapping)
    gsea <- run_pseudobulk_hallmark(pb, reference, file.path(directory, "gsea", "hallmark"), label,
      excluded_features = gsea_exclusions[[cluster]])
    s <- data.frame(cluster, cells = sum(pb$sample_info$n_cells), samples = ncol(pb$counts),
      pigs = nlevels(pb$sample_info$pig), min_cells = min(pb$sample_info$n_cells),
      candidate_genes = nrow(pb$counts), heatmap_genes = plots$n_heatmap_genes, pb$summary,
      up_fdr_logfc = vapply(pb$results, function(x) sum(x$FDR <= 0.05 & x$logFC >= 0.5), integer(1)),
      down_fdr_logfc = vapply(pb$results, function(x) sum(x$FDR <= 0.05 & x$logFC <= -0.5), integer(1)),
      hallmark_status = gsea$status, row.names = NULL)
    if (!is.null(gsea$summary)) s <- merge(s, gsea$summary, by = "contrast", sort = FALSE)
    summaries[[cluster]] <- s
    cluster_info[cluster_info == "Status: running"] <- paste("Status: DE/diagnostics/plots completed; Hallmark", gsea$status)
    writeLines(c(cluster_info, paste("Finished:", Sys.time()),
      paste("Cells:", sum(pb$sample_info$n_cells), "; samples:", ncol(pb$counts),
        "; pigs:", nlevels(pb$sample_info$pig), "; min cells per sample:", min(pb$sample_info$n_cells)),
      paste("Candidate genes:", length(common_genes), "; low-expression genes removed:", sum(!pb$keep),
        "; tested:", sum(pb$keep), "; genes outside original common-layer intersection:", nrow(object) - length(common_genes)),
      paste("Heatmap genes:", plots$n_heatmap_genes),
      "Validation: one saved-output review and figure inspection pending."), file.path(directory, ".INFO"))
    cat(label, "completed; Hallmark", gsea$status, "\n")
    print(s, row.names = FALSE)
  }
  # Some clusters may require an additional mapping decision before GSEA.
  summary <- dplyr::bind_rows(summaries)
  write.csv(summary, file.path(output_root, "cluster_summary.csv"), row.names = FALSE)
  info[info == "Status: running"] <- if (all(summary$hallmark_status == "completed"))
    "Status: execution completed; one saved-output and figure review pending." else
    "Status: DE and plots completed; some Hallmark analyses blocked by identifier collisions."
  writeLines(c(info, paste("Finished:", Sys.time()), paste("Original count layers:", length(count_layers),
    "; candidate genes:", length(common_genes)), "Summary: cluster_summary.csv"), info_path)
  finished <- TRUE
  invisible(summary)
}

if (sys.nframe() == 0L) run_pseudobulk(
  input = "data/job-44542861/subsets/Endothelial/job-44831031/subclustered.rds",
  subset_label = "Endothelial", cluster_col = "harmony_snn_res.0.3", clusters = 1:5,
  output_root = "results/job-44542861/clusters_plots/subset/Endothelial/job-44831031/DEG",
  gsea_exclusions = list("2" = c("TNFSF15", "ENSSSCG00000034087"),
    "5" = c("TNFSF15", "ENSSSCG00000034087")))
