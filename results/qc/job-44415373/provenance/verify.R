# Independent saved-output checks against the pre-submission parameter capture.
suppressPackageStartupMessages(library(Seurat))
run_id <- "job-44415373"
capture_dir <- ".amp/in/qc-mad3-run"
params <- readRDS(file.path(capture_dir, "params.rds"))
samples <- readLines(file.path(capture_dir, "samples.txt"))
labeled <- readRDS(file.path("data/qc_labeled_data", run_id, "labeled_concatenated.rds"))
clean <- readRDS(file.path("data/clean_concatenated_data", run_id, "clean_concatenated_scrublet.rds"))
metadata <- labeled[[]]
stopifnot(
  identical(labeled@misc$qc$run_id, run_id),
  identical(labeled@misc$qc$params, params),
  identical(clean@misc$qc, labeled@misc$qc),
  setequal(labeled@misc$qc$sample_order, samples), length(samples) == 23L,
  !anyDuplicated(colnames(labeled)), !anyDuplicated(rownames(labeled)),
  identical(rownames(metadata), colnames(labeled)), !anyNA(metadata),
  identical(rownames(clean[[]]), colnames(clean)), !anyNA(clean[[]]),
  setequal(colnames(clean), rownames(metadata)[metadata$keep_scrublet]),
  !any(c("doublet_class", "doublet_score", "keep_doubletfinder", "job_id", "run_id") %in% colnames(metadata)),
  !any(grepl("^(pANN_|DF[.])", colnames(metadata))),
  all(clean$keep_cell), all(clean$doublet_filter_method == "Scrublet"),
  !any(c("qc_umap_1", "qc_umap_2") %in% colnames(clean[[]])),
  length(clean@reductions) == 0L, length(labeled@reductions) == 0L,
  length(VariableFeatures(clean)) == 0L, length(VariableFeatures(labeled)) == 0L,
  setequal(Layers(clean[["RNA"]]), paste0("counts.", samples)),
  setequal(Layers(labeled[["RNA"]]), paste0("counts.", samples)),
  !("DoubletFinder" %in% loadedNamespaces())
)
summary <- list()
f <- params$filtering
for (sample in labeled@misc$qc$sample_order) {
  cat("Checking", sample, "\n")
  sample_qc <- labeled@misc$qc$samples[[sample]]
  input <- readRDS(file.path("data/soupx", sample, "job-44409807", "corrected_counts.rds"))
  stopifnot(!anyDuplicated(colnames(input)), !anyDuplicated(rownames(input)),
    all(is.finite(input@x)), all(input@x >= 0), all(input@x == round(input@x)))
  colnames(input) <- paste(sample, colnames(input), sep = "_")
  idx <- match(params$metrics$mitochondrial_genes, rownames(input))
  stopifnot(!anyNA(idx))
  rownames(input)[idx] <- paste0("MT-", params$metrics$mitochondrial_genes)
  m <- metadata[colnames(input), , drop = FALSE]
  counts <- LayerData(labeled, assay = "RNA", layer = paste0("counts.", sample))
  stopifnot(setequal(rownames(counts), rownames(input)),
    setequal(colnames(counts), colnames(input)),
    Matrix::nnzero(counts[rownames(input), colnames(input)] - input) == 0L)
  totals <- Matrix::colSums(input)
  genes <- Matrix::colSums(input > 0)
  complexity <- log10(genes) / log10(totals)
  ribo <- 100 * Matrix::colSums(input[grepl(params$metrics$ribosomal_pattern, rownames(input)), , drop = FALSE]) / totals
  mito <- 100 * Matrix::colSums(input[grepl(params$metrics$mitochondrial_pattern, rownames(input)), , drop = FALSE]) / totals
  feature_max <- min(f$max_features_cap, median(genes) + f$feature_mad_multiplier * mad(genes, constant = 1.4826))
  initial_bad <- genes <= f$min_features | genes > feature_max | complexity <= f$min_log10_genes_per_umi
  mt_max <- min(f$max_mito_percent_cap, median(mito[!initial_bad]) + f$mito_mad_multiplier * mad(mito[!initial_bad], constant = 1.4826))
  ribo_max <- median(ribo) + f$ribo_mad_multiplier * mad(ribo, constant = 1.4826)
  bad <- initial_bad | mito > mt_max | ribo > ribo_max
  bd <- params$doublets$bd_multiplet_table
  rate <- approx(bd$cells, bd$rate, xout = ncol(input), rule = 2)$y / 100
  condition <- strsplit(sample, "-", fixed = TRUE)[[1]][[2]]
  stopifnot(
    all(m$sample == sample), all(m$orig.ident == sample),
    all(m$pig == strsplit(sample, "-", fixed = TRUE)[[1]][[1]]),
    all(m$pressure == switch(substr(condition, 1, 1), C = "None", P = "positive", N = "negative")),
    all(m$time_point == if (condition == "C") "T0H" else paste0("T", substring(condition, 2), "H")),
    isTRUE(all.equal(unname(m$nCount_RNA), unname(totals))),
    isTRUE(all.equal(unname(m$nFeature_RNA), unname(genes))),
    isTRUE(all.equal(unname(m$log10GenesPerUMI), unname(complexity))),
    isTRUE(all.equal(unname(m$percent.ribo), unname(ribo))),
    isTRUE(all.equal(unname(m$percent.mt), unname(mito))),
    all(m$flag_low_features == (genes <= f$min_features)),
    all(m$flag_high_features == (genes > feature_max)),
    all(m$flag_low_complexity == (complexity <= f$min_log10_genes_per_umi)),
    all(m$flag_high_mt == (mito > mt_max)),
    all(m$flag_high_ribo == (ribo > ribo_max)),
    all(m$ribo_status == ifelse(ribo > ribo_max, "High ribosomal", "Within ribosomal threshold")),
    all(m$qc_outlier == bad), all(m$keep_qc == !bad),
    all(is.finite(m$doublet_scores)),
    all(m$predicted_doublets == (m$doublet_scores > 0.15)),
    all(m$keep_scrublet == (!bad & m$doublet_scores <= 0.15)),
    isTRUE(all.equal(sample_qc$feature_max, unname(feature_max))),
    isTRUE(all.equal(sample_qc$mt_max, unname(mt_max))),
    isTRUE(all.equal(sample_qc$ribo_max, unname(ribo_max))),
    identical(sample_qc$scrublet, list(expected_doublet_rate = rate,
      min_counts = 3L, n_prin_comps = 30L, threshold = 0.15, seed = 0L)),
    all(is.finite(as.matrix(m[, c("qc_umap_1", "qc_umap_2")]))))
  kept <- colnames(input)[!bad & m$doublet_scores <= 0.15]
  detected <- Matrix::rowSums(input[, kept, drop = FALSE] > 0)
  expected <- input[detected >= f$min_cells_per_feature, kept, drop = FALSE]
  observed <- LayerData(clean, assay = "RNA", layer = paste0("counts.", sample))
  meta_cols <- colnames(m)[!grepl("^qc_umap_", colnames(m))]
  stopifnot(
    any(detected == f$min_cells_per_feature - 1), any(detected == f$min_cells_per_feature),
    setequal(rownames(observed), rownames(expected)), setequal(colnames(observed), kept),
    Matrix::nnzero(observed - expected[rownames(observed), colnames(observed)]) == 0L,
    identical(clean[[]][kept, meta_cols, drop = FALSE], m[kept, meta_cols, drop = FALSE]))
  summary[[sample]] <- data.frame(sample = sample, input_cells = ncol(input),
    qc_only = sum(!bad), scrublet_doublets = sum(m$doublet_scores > 0.15),
    retained = length(kept), retained_percent = 100 * length(kept) / ncol(input),
    retained_genes = nrow(expected), feature_max = feature_max,
    mt_max = mt_max, ribo_max = ribo_max, ribo_flags = sum(ribo > ribo_max),
    expected_doublet_rate = rate)
  rm(input, counts, expected, observed)
  gc()
}
summary <- do.call(rbind, summary)
print(summary, row.names = FALSE)
output_dir <- file.path("results/qc", run_id)
expected_figures <- c(paste0("qc_before_after_", c("genes", "umi", "mitochondrial", "complexity", "ribosomal"), ".png"),
  paste0("qc_umap_", c("scrublet", "ribosomal", "outliers"), ".png"))
figures <- list.files(output_dir, pattern = "[.](png|svg)$")
stopifnot(setequal(figures, expected_figures),
  all(file.info(file.path(output_dir, figures))$size > 0),
  sum(summary$input_cells) == 299393L)
for (figure in figures) {
  image <- png::readPNG(file.path(output_dir, figure))
  stopifnot(length(dim(image)) == 3L, all(dim(image)[1:2] > 1000L))
  rm(image)
  gc()
}
write.csv(summary, file.path(output_dir, "validation_summary.csv"), row.names = FALSE)
saveRDS(list(params = params, summary = summary, figures = figures,
  validated_at = Sys.time(), session = sessionInfo()), file.path(output_dir, "validation.rds"))
cat("PASS: all 23 saved input counts, metadata/IDs, independent sample-specific QC flags, Scrublet-only calls/retention, >=10-cell per-sample gene filtering (9/10 boundary), counts-only checkpoints and eight readable PNGs.\n")
