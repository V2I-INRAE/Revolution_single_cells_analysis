# Rendering from frozen analysis outputs only; never fit models or run enrichment.
source("R/utils/save_plots.R")
source("R/DEG/pseudobulk_diagnostics.R")
source("R/DEG/pseudobulk_plots.R")
source("R/DEG/pseudobulk_hallmark.R")

# Construct only production targets: previews are never discovered recursively.
saved_pseudobulk_figure_paths <- function(directory, contrasts, hallmark_contrasts = character(),
                                          diagnostics_only = FALSE) {
  stems <- file.path("diagnostics", c("sample_mds", "sample_qc"))
  if (!diagnostics_only) {
    stems <- c(stems, file.path("plots", c("heatmap_groups", paste0("volcano_", contrasts))))
    if (length(hallmark_contrasts)) {
      stems <- c(stems,
        file.path("gsea", "hallmark", paste0(hallmark_contrasts, "_gsea_dotplot")),
        file.path("gsea", "hallmark", paste0(hallmark_contrasts, "_running_scores")))
    }
  }
  file.path(directory, paste0(stems, ".pdf"))
}

read_saved_pseudobulk_gsea <- function(input_dir, contrasts) {
  directory <- file.path(input_dir, "gsea", "hallmark")
  path <- file.path(directory, "gsea_results.rds")
  if (file.exists(file.path(directory, "unresolved_collisions.csv")) || !file.exists(path)) {
    status <- if (file.exists(file.path(directory, "unresolved_collisions.csv"))) "blocked by identifier collisions" else "missing saved results"
    stop("Hallmark GSEA ", status, " in ", input_dir,
      "; full regeneration unavailable. Use --diagnostics-only for an explicitly partial replacement.")
  }
  gsea <- readRDS(path)
  if (!is.null(gsea$status) && !identical(gsea$status, "completed")) {
    stop("Hallmark GSEA status: ", gsea$status, "; use --diagnostics-only for a partial replacement.")
  }
  if (!identical(names(gsea$results), contrasts) ||
      !identical(names(gsea$plot_tables), contrasts) || !identical(names(gsea$ranks), contrasts)) {
    stop("Hallmark saved results are incomplete or contrasts differ from the frozen fit in ", input_dir,
      "; full regeneration unavailable. Use --diagnostics-only for a partial replacement.")
  }
  gsea
}

render_saved_pseudobulk <- function(input_dir, output_dir, label, fdr_cutoff, logfc_cutoff, seed, diagnostics_only = FALSE) {
  pb <- readRDS(file.path(input_dir, "pseudobulk_fit.rds"))
  mds <- readRDS(file.path(input_dir, "diagnostics", "mds.rds"))
  mds_data <- read.csv(file.path(input_dir, "diagnostics", "mds_coordinates.csv"))
  required <- c("sample", "pig", "condition", "n_cells", "library_size", "tmm_factor",
    "effective_library_size", "MDS1", "MDS2")
  stopifnot(all(required %in% names(mds_data)))
  stopifnot(!anyDuplicated(mds_data$sample), !anyNA(mds_data[, required, drop = FALSE]),
    identical(mds_data$sample, as.character(mds$sample_info$sample)),
    identical(mds_data$sample, as.character(pb$sample_info$sample)),
    isTRUE(all.equal(mds_data$MDS1, unname(mds$mds$x), check.attributes = FALSE)),
    isTRUE(all.equal(mds_data$MDS2, unname(mds$mds$y), check.attributes = FALSE)))
  condition_levels <- c("Control", "P4", "N4", "P10", "N10")
  stopifnot(!anyDuplicated(pb$sample_info$sample),
    !anyNA(pb$sample_info[, c("sample", "condition", "pig")]),
    all(pb$sample_info$condition %in% condition_levels),
    identical(as.character(mds_data$condition), as.character(mds$sample_info$condition)),
    identical(as.character(mds_data$condition), as.character(pb$sample_info$condition)),
    identical(as.character(mds_data$pig), as.character(mds$sample_info$pig)),
    identical(as.character(mds_data$pig), as.character(pb$sample_info$pig)))
  mds_data$condition <- factor(mds_data$condition,
    levels = condition_levels[condition_levels %in% mds_data$condition])
  mds_data$pig <- factor(mds_data$pig, levels = levels(factor(pb$sample_info$pig)))

  if (diagnostics_only) {
    render_pseudobulk_diagnostics(mds_data, mds$mds, nrow(pb$dge),
      file.path(output_dir, "diagnostics"), label, seed)
    return(invisible(c(samples = nrow(mds_data), genes = nrow(pb$dge))))
  }
  sample_info <- pb$sample_info
  sample_info$condition <- factor(sample_info$condition,
    levels = condition_levels[condition_levels %in% sample_info$condition])
  sample_info$pig <- factor(sample_info$pig)
  sample_info <- sample_info[order(sample_info$condition, sample_info$pig), ]
  log2_cpm <- edgeR::cpm(pb$dge, log = TRUE, prior.count = 2)
  # Preserve the saved run's display selection, including historical DEG-only heatmaps.
  saved <- read.csv(file.path(input_dir, "plots", "heatmap_zscores.csv"), check.names = FALSE)
  saved$gene <- as.character(saved$gene)
  stopifnot(!anyNA(saved$gene), !anyDuplicated(saved$gene),
    all(saved$gene %in% rownames(log2_cpm)))
  selected <- rownames(log2_cpm)[rownames(log2_cpm) %in% saved$gene]
  z_scores <- log2_cpm[selected, sample_info$sample, drop = FALSE]
  if (length(selected)) z_scores <- t(scale(t(z_scores)))
  stopifnot(all(is.finite(z_scores)))
  row_hc <- if (length(selected) > 1) {
    hclust(dist(z_scores, method = "euclidean"), method = "complete")
  } else FALSE
  gene_order <- if (length(selected) > 1) rownames(z_scores)[row_hc$order] else selected
  saved_values <- as.matrix(saved[-1])
  storage.mode(saved_values) <- "double"
  stopifnot(identical(saved$gene, gene_order),
    identical(names(saved)[-1], as.character(sample_info$sample)),
    isTRUE(all.equal(unname(saved_values),
      unname(z_scores[saved$gene, , drop = FALSE]), tolerance = 1e-12)))

  volcano <- read.csv(file.path(input_dir, "plots", "volcano_data.csv"))
  stopifnot(!anyNA(volcano), setequal(unique(volcano$contrast), names(pb$results)),
    all(table(volcano$contrast) == nrow(pb$dge)))
  for (contrast in names(pb$results)) {
    tab <- volcano[volcano$contrast == contrast, ]
    original <- pb$results[[contrast]]
    stopifnot(identical(tab$gene, original$gene),
      isTRUE(all.equal(tab$logFC, original$logFC)),
      isTRUE(all.equal(tab$FDR, original$FDR)))
  }


  gsea <- read_saved_pseudobulk_gsea(input_dir, names(pb$results))
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
      pb$sample_info, contrast, file.path(output_dir, "gsea", "hallmark"),
      label, gsea$parameters$pathway_fdr, seed)
  }
  invisible(c(samples = nrow(sample_info), genes = nrow(pb$dge), heatmap_genes = nrow(z_scores)))
}
