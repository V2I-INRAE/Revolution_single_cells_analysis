args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || length(args) > 2L) stop(
  "Usage: Rscript scripts/soupx/review_soupx.R <soupx_data_dir> [validated_cellbender_run_dir]")
source("scripts/soupx/io.R")
source("scripts/soupx/plots.R")
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(scplotter)
})
future::plan("sequential")
input <- normalizePath(args[1])
source_run <- basename(input)
sample_id <- basename(dirname(input))
stopifnot("Marker-to-cluster mappings are specific to the REVO30-P4 pilot" =
  sample_id == "REVO30-P4", source_run == "job-44404963")
source_results <- file.path("results", "soupx", source_run, sample_id)
source_metadata <- jsonlite::read_json(file.path(source_results, "run.json"),
  simplifyVector = TRUE)
stopifnot(source_metadata$status == "validated_pending_qc")
cb_dir <- if (length(args) == 2L) normalizePath(args[2]) else NULL
if (!is.null(cb_dir)) {
  cb_metadata <- jsonlite::read_json(file.path(cb_dir, "run.json"), simplifyVector = TRUE)
  cb_review <- jsonlite::read_json(file.path(cb_dir, "qc_review.json"), simplifyVector = TRUE)
  cb_validation <- jsonlite::read_json(file.path(cb_dir, "validation.json"), simplifyVector = TRUE)
  stopifnot(cb_metadata$sample == sample_id, cb_metadata$returncode == 0,
    cb_metadata$status %in% c("validated_pending_qc", "accepted"),
    cb_validation$numerical_validation == "passed")
  cb_accepted <- isTRUE(cb_review$elbo_converged) && cb_review$qc_status == "accepted"
}
job_id <- Sys.getenv("SLURM_JOB_ID")
run_id <- if (nzchar(job_id)) paste0("job-", job_id) else
  paste(format(Sys.time(), "%Y%m%d-%H%M%S"), Sys.getpid(), sep = "-")
output <- file.path("results", "soupx", run_id, sample_id)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
if (!dir.create(output)) stop("Existing review directory: ", output)

channel <- readRDS(file.path(input, "soupx_channel.rds"))
matrices <- list(Raw = channel$toc,
  SoupX = readRDS(file.path(input, "corrected_counts.rds")))
validate_soupx_counts(matrices$SoupX, matrices$Raw)
clustering <- readRDS(file.path(input, "preliminary_clustering.rds"))
clusters <- clustering$clusters[colnames(matrices$Raw)]
stopifnot(!anyNA(clusters), identical(names(clusters), colnames(matrices$Raw)))
feature_map <- read.delim(file.path(input, "feature_map.tsv"))
if (!is.null(cb_dir)) {
  # Explicit FPR, full matrix and original Rhapsody calls; extra CB calls are not added.
  h5 <- file.path(cb_dir, paste0(sample_id, "_cellbender_FPR_0.01.h5"))
  # Seurat's reader mistakes CellBender latent groups for expression matrices.
  # Preserve the venv symlink: normalizePath resolves it to the base interpreter.
  reticulate::use_python(file.path(getwd(), ".venv-scvi", "bin", "python"), required = TRUE)
  cb_input <- reticulate::import("cellbender.remove_background.data.io")$load_data(h5)
  cb <- as(t(cb_input$matrix), "CsparseMatrix")
  dimnames(cb) <- list(as.character(cb_input$gene_ids), as.character(cb_input$barcodes))
  rm(cb_input)
  rows <- match(feature_map$feature_id, rownames(cb))
  columns <- match(colnames(matrices$Raw), colnames(cb))
  stopifnot(!anyNA(rows), !anyNA(columns), !anyDuplicated(rownames(cb)),
    nrow(cb) == nrow(matrices$Raw))
  cb <- cb[rows, columns, drop = FALSE]
  dimnames(cb) <- dimnames(matrices$Raw)
  validate_soupx_counts(cb, matrices$Raw)
  matrices$CellBender <- cb
  comparison <- data.frame(cell = colnames(cb), cluster = clusters,
    raw_molecules = colSums(matrices$Raw),
    soupx_molecules = colSums(matrices$SoupX), cellbender_molecules = colSums(cb))
  cb_calls <- readLines(file.path(cb_dir, paste0(sample_id, "_cellbender_cell_barcodes.csv")))
  stopifnot(!anyDuplicated(cb_calls))
  # Compare the fixed original cell set, independently of CellBender's cell calls.
  comparison$cellbender_called <- comparison$cell %in% cb_calls
  connection <- gzfile(file.path(output, "matched_cell_counts.csv.gz"), "wt")
  write.csv(comparison, connection, row.names = FALSE)
  close(connection)
  agreement <- list(cells = ncol(cb), features = nrow(cb),
    cellbender_source = cb_dir, cellbender_fpr = 0.01,
    cellbender_qc_accepted = cb_accepted, cellbender_qc_status = cb_review$qc_status,
    cellbender_elbo_converged = cb_review$elbo_converged,
    original_cells_not_called_by_cellbender = comparison$cell[!comparison$cellbender_called],
    additional_cellbender_calls_not_in_comparison = sum(!cb_calls %in% colnames(cb)),
    cellbender_h5_md5 = unname(tools::md5sum(h5)),
    soupx_fraction_removed = 1 - sum(matrices$SoupX) / sum(matrices$Raw),
    cellbender_fraction_removed = 1 - sum(cb) / sum(matrices$Raw),
    corrected_cell_total_pearson = cor(comparison$soupx_molecules,
      comparison$cellbender_molecules),
    removed_cell_total_pearson = cor(comparison$raw_molecules - comparison$soupx_molecules,
      comparison$raw_molecules - comparison$cellbender_molecules))
  jsonlite::write_json(agreement, file.path(output, "matched_agreement.json"),
    pretty = TRUE, auto_unbox = TRUE, digits = NA, null = "null")
  write.csv(data.frame(gene = rownames(cb), raw_molecules = rowSums(matrices$Raw),
    soupx_molecules = rowSums(matrices$SoupX), cellbender_molecules = rowSums(cb)),
    file.path(output, "matched_gene_counts.csv"), row.names = FALSE)
}

# These qualitative, post-hoc marker patterns guide QC, not final cell annotation.
# Human-PBMC-trained BD predictions are not a pig-lung validation reference.
supported <- c("0", "1", "2", "3", "4", "5", "7", "8", "9", "10", "11",
  "12", "17", "18", "19")
epithelial <- c("10", "12", "17", "19")
panels <- list(
  T_like = list(genes = c("CD3D", "CD3E"), positive = c("0", "1"), exclude = character()),
  cytotoxic_lymphoid_like = list(genes = c("GZMA", "PRF1"), positive = "1", exclude = "0"),
  B_like = list(genes = c("MS4A1", "CD79B"), positive = "5", exclude = "3"),
  plasma_like = list(genes = c("JCHAIN", "MZB1"), positive = "3", exclude = "5"),
  myeloid_like = list(genes = c("CD68", "C1QA", "CSF3R", "TREM1"),
    positive = c("2", "7", "11"), exclude = character()),
  endothelial_like = list(genes = c("PECAM1", "CLDN5", "VWF"),
    positive = c("8", "9", "18"), exclude = character()),
  fibroblast_like = list(genes = c("COL1A1", "DCN", "LUM"), positive = "4", exclude = character()),
  ciliated_like = list(genes = c("FOXJ1", "RSPH1"), positive = "10", exclude = epithelial),
  secretory_airway_like = list(genes = c("SCGB1A1", "SCGB3A2"),
    positive = c("10", "12"), exclude = epithelial),
  alveolar_epithelial_like = list(genes = c("SFTPC", "SFTPB", "AGER"),
    positive = "17", exclude = epithelial),
  POU2F3_epithelial_like = list(genes = "POU2F3", positive = "19", exclude = epithelial))
jsonlite::write_json(list(panels = panels,
  unresolved_clusters = setdiff(unique(clusters), supported),
  note = "Post-hoc raw-marker patterns; unresolved and potentially overlapping lineages excluded from negative references. No cells excluded from corrected data."),
  file.path(output, "provisional_marker_panels.json"), pretty = TRUE, auto_unbox = TRUE)
genes <- unique(unlist(lapply(panels, `[[`, "genes")))
missing <- setdiff(genes, rownames(matrices$Raw))
writeLines(missing, file.path(output, "missing_panel_genes.txt"))
genes <- intersect(genes, rownames(matrices$Raw))

stats <- list()
for (panel_name in names(panels)) {
  panel <- panels[[panel_name]]
  pos <- clusters %in% panel$positive
  neg <- clusters %in% setdiff(supported, c(panel$positive, panel$exclude))
  for (gene in intersect(panel$genes, genes)) {
    for (method in names(matrices)) {
      values <- as.numeric(matrices[[method]][gene, ])
      for (reference in c("positive", "negative")) {
        mask <- if (reference == "positive") pos else neg
        stats[[length(stats) + 1L]] <- data.frame(panel = panel_name, gene = gene,
          method = method, reference = reference, cells = sum(mask),
          molecules = sum(values[mask]), mean_molecules = mean(values[mask]),
          expressing_fraction = mean(values[mask] > 0))
      }
    }
  }
}
stats <- do.call(rbind, stats)
write.csv(stats, file.path(output, "marker_reference_counts.csv"), row.names = FALSE)
raw_stats <- stats[stats$method == "Raw", ]
key <- function(data) paste(data$panel, data$gene, data$reference, sep = ":")
index <- match(key(stats), key(raw_stats))
stats$retained_fraction <- stats$molecules / raw_stats$molecules[index]
stats$detection_change_percentage_points <- 100 *
  (stats$expressing_fraction - raw_stats$expressing_fraction[index])
write.csv(stats, file.path(output, "marker_preservation.csv"), row.names = FALSE)

plot_labels <- names(matrices)
if (!is.null(cb_dir) && !cb_accepted) {
  plot_labels[plot_labels == "CellBender"] <- "CellBender FPR 0.01 (exploratory; QC not accepted)"
}
plot_metadata <- plot_soupx_markers(matrices, clusters, sample_id, output,
  genes = genes, labels = plot_labels)

label_path <- file.path("data", "raw_data", sample_id,
  paste0(sample_id, "_cell_type_experimental.csv"))
labels <- read.csv(label_path, colClasses = "character")
stopifnot(!anyDuplicated(labels$Cell_Index))
label_index <- match(names(clusters), labels$Cell_Index)
stopifnot(!anyNA(label_index))
write.csv(as.data.frame(table(cluster = clusters,
  bd_prediction = labels$Cell_Type_Experimental[label_index])),
  file.path(output, "bd_prediction_cluster_contingency.csv"), row.names = FALSE)
jsonlite::write_json(list(status = "computed_pending_visual_review",
  source_data = input, source_metadata = normalizePath(file.path(source_results, "run.json")),
  cellbender_source = cb_dir, cellbender_comparison = if (is.null(cb_dir))
    "not_requested" else if (cb_accepted) "accepted_comparator_on_original_rhapsody_cells"
    else "exploratory_unaccepted_comparator_on_original_rhapsody_cells",
  fixed_clusters = TRUE, unchanged_count_matrices = TRUE,
  marker_panels_post_hoc = TRUE, unresolved_clusters_excluded_from_negative_references = TRUE,
  bd_labels_not_ground_truth = TRUE, missing_genes = missing,
  plot_library = "scplotter::FeatureStatPlot", plot_layer = "data",
  plot_normalization = "Seurat LogNormalize, scale factor 10000; common color scale",
  marker_plots = plot_metadata,
  versions = list(Seurat = as.character(packageVersion("Seurat")),
    scplotter = as.character(packageVersion("scplotter"))),
  plot_script_md5 = unname(tools::md5sum("scripts/soupx/plots.R")),
  script_md5 = unname(tools::md5sum("scripts/soupx/review_soupx.R"))),
  file.path(output, "review.json"), pretty = TRUE, auto_unbox = TRUE, null = "null")
message("Review tables and figures saved to ", output)
