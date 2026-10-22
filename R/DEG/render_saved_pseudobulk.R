# Rendering from frozen analysis outputs only; never fit models or run enrichment.
source("R/utils/save_plots.R")
source("R/DEG/pseudobulk_diagnostics.R")
source("R/DEG/pseudobulk_plots.R")
source("R/DEG/pseudobulk_hallmark.R")

render_saved_pseudobulk <- function(input_dir, output_dir, label, fdr_cutoff, logfc_cutoff, seed, diagnostics_only = FALSE) {
  pb <- readRDS(file.path(input_dir, "pseudobulk_fit.rds"))
  mds <- readRDS(file.path(input_dir, "diagnostics", "mds.rds"))
  mds_data <- read.csv(file.path(input_dir, "diagnostics", "mds_coordinates.csv"))
  stopifnot(!anyDuplicated(mds_data$sample), !anyNA(mds_data),
    identical(mds_data$sample, as.character(mds$sample_info$sample)),
    identical(mds_data$sample, as.character(pb$sample_info$sample)),
    isTRUE(all.equal(mds_data$MDS1, unname(mds$mds$x), check.attributes = FALSE)),
    isTRUE(all.equal(mds_data$MDS2, unname(mds$mds$y), check.attributes = FALSE)))
  mds_data$condition <- factor(mds_data$condition, levels = levels(mds$sample_info$condition))
  mds_data$pig <- factor(mds_data$pig, levels = levels(mds$sample_info$pig))


  if (diagnostics_only) {
    render_pseudobulk_diagnostics(mds_data, mds$mds, nrow(pb$dge),
      file.path(output_dir, "diagnostics"), label, seed)
    return(invisible(c(samples = nrow(mds_data), genes = nrow(pb$dge))))
  }
  sample_info <- pb$sample_info
  sample_info$condition <- factor(sample_info$condition, levels = c("Control", "P4", "N4"))
  sample_info <- sample_info[order(sample_info$condition, sample_info$pig), ]
  log2_cpm <- edgeR::cpm(pb$dge, log = TRUE, prior.count = 2)
  passes <- vapply(pb$results, function(tab) {
    rownames(log2_cpm) %in% tab$gene[tab$FDR <= fdr_cutoff & abs(tab$logFC) >= logfc_cutoff]
  }, logical(nrow(log2_cpm)))
  selected <- rownames(log2_cpm)[rowSums(passes) > 0]
  z_scores <- t(scale(t(log2_cpm[selected, sample_info$sample, drop = FALSE])))
  row_hc <- hclust(dist(z_scores, method = "euclidean"), method = "complete")
  saved <- read.csv(file.path(input_dir, "plots", "heatmap_zscores.csv"), check.names = FALSE)
  stopifnot(identical(saved$gene, rownames(z_scores)[row_hc$order]),
    identical(names(saved)[-1], sample_info$sample),
    isTRUE(all.equal(unname(as.matrix(saved[-1])),
      unname(z_scores[saved$gene, , drop = FALSE]), tolerance = 1e-12)))

  volcano <- read.csv(file.path(input_dir, "plots", "volcano_data.csv"))
  stopifnot(!anyNA(volcano), all(table(volcano$contrast) == nrow(pb$dge)))
  for (contrast in names(pb$results)) {
    tab <- volcano[volcano$contrast == contrast, ]
    original <- pb$results[[contrast]]
    stopifnot(identical(tab$gene, original$gene),
      isTRUE(all.equal(tab$logFC, original$logFC)),
      isTRUE(all.equal(tab$FDR, original$FDR)))
  }


  gsea <- readRDS(file.path(input_dir, "gsea", "hallmark", "gsea_results.rds"))
  stopifnot(identical(names(gsea$results), names(pb$results)),
    identical(names(gsea$plot_tables), names(pb$results)),
    identical(gsea$parameters$seed, seed))
  for (contrast in names(gsea$results)) {
    display <- read.csv(file.path(input_dir, "gsea", "hallmark",
      paste0(contrast, "_displayed_pathways.csv")))
    original <- gsea$plot_tables[[contrast]]
    stopifnot(identical(display$pathway, original$pathway))
    for (column in c("NES", "padj", "size", "neg_log10_fdr")) {
      stopifnot(isTRUE(all.equal(display[[column]], original[[column]], tolerance = 1e-12)))
    }
  }
  cat(label, ": inputs verified;", nrow(pb$dge), "genes;", nrow(z_scores), "heatmap genes\n")
  render_pseudobulk_diagnostics(mds_data, mds$mds, nrow(pb$dge),
    file.path(output_dir, "diagnostics"), label, seed)
  render_pseudobulk_heatmap(z_scores, row_hc, sample_info, file.path(output_dir, "plots"),
    label, fdr_cutoff, logfc_cutoff)
  render_pseudobulk_volcano(volcano, file.path(output_dir, "plots"), label,
    fdr_cutoff, logfc_cutoff, seed)
  for (contrast in names(gsea$results)) {
    render_pseudobulk_hallmark(gsea$results[[contrast]], gsea$plot_tables[[contrast]],
      gsea$ranks[[contrast]], gsea$pathways, gsea$database_version,
      nlevels(pb$sample_info$pig), contrast, file.path(output_dir, "gsea", "hallmark"),
      label, gsea$parameters$pathway_fdr, seed)
  }
  invisible(c(samples = nrow(sample_info), genes = nrow(pb$dge), heatmap_genes = nrow(z_scores)))
}
