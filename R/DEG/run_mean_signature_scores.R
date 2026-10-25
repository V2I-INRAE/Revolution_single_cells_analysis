# Arithmetic mean of EXISTING RNA LogNormalize values, including expression zeros.
# Not AUCell, AddModuleScore, mean raw counts, or log(mean expression).
# Source from project root; sourcing does not execute a run.
# run_mean_signature_scores(
#   input = "data/job-44542861/subsets/Endothelial/job-44831031/subclustered.rds",
#   mapping_file = "results/misc/Ensembl_ID/signature-scoring-v2/pig_mapping.csv",
#   universe_file = "results/misc/Ensembl_ID/signature-scoring/common_count_universe.rds",
#   output_root = "results/Endothelial_cluster1_mean_scores_new",
#   cluster_col = "harmony_snn_res.0.3", clusters = "1")
# Fixed mapping/universe are verified across ALL original sample layers before
# selecting cells. For a recorded real-cell smoke only, supply cells+selection_note
# instead of cluster_col/clusters. No stochastic operation is used for scoring.
source("R/DEG/signatures.R")
source("R/DEG/signature_gene_sets.R")

run_mean_signature_scores <- function(input, mapping_file, universe_file, output_root,
    cluster_col = NULL, clusters = NULL, batch_size = 1000L,
    cells = NULL, selection_note = NULL) {
  if (file.exists(output_root)) stop("output_root exists; use a new directory.")
  for (p in c("SeuratObject", "Matrix")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
  if (batch_size < 1L || batch_size != as.integer(batch_size)) stop("Use a positive integer batch_size.")
  if (xor(is.null(cluster_col), is.null(clusters))) stop("Supply both cluster_col and clusters.")
  if (!is.null(cluster_col) && !is.null(cells)) stop("Cluster selection scores ALL selected cells; do not also supply cells.")
  if (!is.null(cluster_col)) {
    clusters <- as.character(clusters)
    if (length(cluster_col) != 1L || is.na(cluster_col) || !nzchar(cluster_col) || !length(clusters) ||
        anyNA(clusters) || anyDuplicated(clusters) || any(!nzchar(clusters))) stop("Invalid cluster column/labels.")
    selection_note <- paste0("All cells in ", cluster_col, " clusters ", paste(clusters, collapse = ","),
      "; verify full original sample-layer universe first; no cell/sample subsampling.")
  }
  if (!is.null(cells) && (is.null(selection_note) || !nzchar(selection_note))) stop("Pilot cell selection requires selection_note.")
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  info_path <- file.path(output_root, ".INFO")
  writeLines(c("# Mean log-expression signature scoring", "Status: RUNNING", paste("Started:", Sys.time()),
    paste("Read-only input:", input), paste("Input producer:", basename(dirname(input))),
    paste("Mapping:", mapping_file), paste("Frozen universe:", universe_file),
    paste("Selection:", if (is.null(selection_note)) "all source cells" else selection_note),
    paste("Computational batch_size:", batch_size, "; no stochastic score operation/seed."),
    "Arithmetic per-cell mean of existing RNA data.* LogNormalize values INCLUDING zeros.",
    "Fixed denominator per arm: number of unique retained pig features, not detected genes per cell.",
    "No transformation, ranking, cutoff, controls, subtraction, cell filtering or source-object write-back.",
    "Original human lists unchanged; reviewed mapping exclusions retained; Rooney arms separate."), info_path)
  completed <- FALSE
  on.exit(if (!completed) cat("\nStatus: FAILED; partial outputs not validated.\n", file = info_path, append = TRUE))
  source_before <- file.info(input)[, c("size", "mtime")]
  signature_md5 <- unname(tools::md5sum("R/DEG/signatures.R"))
  universe <- readRDS(universe_file)
  features <- universe$features
  if (!length(features) || anyDuplicated(features)) stop("Invalid frozen feature universe.")
  prepared <- prepare_signature_gene_sets(mapping_file, features, signatures)
  mapping <- prepared$mapping; coverage <- prepared$coverage
  gene_sets <- prepared$gene_sets[lengths(prepared$gene_sets) > 0L]
  write.csv(mapping, file.path(output_root, "mapping_audit.csv"), row.names = FALSE)
  write.csv(coverage, file.path(output_root, "coverage.csv"), row.names = FALSE)
  if (!length(gene_sets)) stop("No nonempty approved gene sets in frozen universe.")
  if (any(coverage$human_coverage < 1)) warning("Incomplete source-human coverage; see coverage.csv.")
  if (any(coverage$small_set_warning)) warning("Small effective sets (<10 genes); warning only, no eligibility filter.")
  object <- readRDS(input)
  metadata <- object[[]]
  required <- c("sample", "pig", "pressure", "time_point")
  if (!all(required %in% names(metadata)) || !identical(rownames(metadata), colnames(object)) ||
      anyDuplicated(rownames(metadata)) || anyNA(metadata[, required, drop = FALSE]) ||
      any(vapply(metadata[, required, drop = FALSE], function(x) any(!nzchar(as.character(x))), logical(1)))) stop("Source metadata/unique cell IDs invalid.")
  panel <- unique(metadata[, required, drop = FALSE])
  if (anyDuplicated(panel$sample) || !setequal(as.character(panel$sample), universe$samples)) stop("Source must retain one original tuple per sample and the full frozen sample panel before selection.")
  normalization <- object@commands[["NormalizeData.RNA"]]
  if (is.null(normalization) || !identical(normalization@assay.used, "RNA") ||
      !identical(normalization@params$normalization.method, "LogNormalize") ||
      length(normalization@params$scale.factor) != 1L || !is.finite(normalization@params$scale.factor) ||
      normalization@params$scale.factor <= 0) stop("Saved RNA LogNormalize command evidence required; do not infer normalization from a data-layer name.")
  normalization_evidence <- list(command = "NormalizeData.RNA", assay = normalization@assay.used,
    parameters = normalization@params, time = normalization@time.stamp, call = normalization@call.string)
  data_layers <- sort(SeuratObject::Layers(object[["RNA"]], search = "^data[.]"), method = "radix")
  count_layers <- sort(SeuratObject::Layers(object[["RNA"]], search = "^counts[.]"), method = "radix")
  if (!length(data_layers) || !setequal(sub("^data[.]", "", data_layers), universe$samples) ||
      !setequal(sub("^counts[.]", "", count_layers), universe$samples)) stop("Original sample-specific data/count layers required for every frozen sample.")
  layer_cells <- lapply(data_layers, function(l) SeuratObject::Cells(object[["RNA"]], layer = l))
  data_features <- lapply(data_layers, function(l) SeuratObject::Features(object[["RNA"]], layer = l))
  count_features <- lapply(count_layers, function(l) SeuratObject::Features(object[["RNA"]], layer = l))
  all_cells <- unlist(layer_cells, use.names = FALSE)
  if (anyDuplicated(all_cells) || !setequal(all_cells, rownames(metadata))) stop("RNA data layers must partition original cells exactly.")
  for (i in seq_along(data_layers)) {
    sample <- sub("^data[.]", "", data_layers[i])
    if (!identical(unique(as.character(metadata[layer_cells[[i]], "sample"])), sample) ||
        !setequal(layer_cells[[i]], SeuratObject::Cells(object[["RNA"]], layer = paste0("counts.", sample)))) stop("RNA data/count layer and sample identity disagree.")
  }
  data_common <- Reduce(intersect, data_features)
  if (!setequal(Reduce(intersect, count_features), features) || !all(features %in% data_common)) stop("Frozen count universe must match source counts and be present in ALL data layers; no zero padding or automatic universe change.")
  # All layer checks above precede cluster or pilot-cell selection.
  if (!is.null(cluster_col)) {
    if (!cluster_col %in% names(metadata) || anyNA(metadata[[cluster_col]]) ||
        any(!nzchar(as.character(metadata[[cluster_col]]))) || !all(clusters %in% as.character(metadata[[cluster_col]]))) stop("Requested source cluster column/labels missing or incomplete.")
    cells <- rownames(metadata)[as.character(metadata[[cluster_col]]) %in% clusters]
  }
  if (is.null(cells)) cells <- rownames(metadata)
  if (!length(cells) || anyDuplicated(cells) || !all(cells %in% rownames(metadata))) stop("Selected cells must be unique original cell IDs.")
  selected <- lapply(layer_cells, function(z) sort(z[z %in% cells], method = "radix"))
  ordered_cells <- unlist(selected, use.names = FALSE)
  metadata <- metadata[ordered_cells, , drop = FALSE]
  scores <- matrix(NA_real_, nrow = length(gene_sets), ncol = length(ordered_cells),
    dimnames = list(names(gene_sets), ordered_cells))
  batches <- list(); batch_id <- 0L
  for (i in seq_along(data_layers)) {
    cc <- selected[[i]]
    if (!length(cc)) next
    for (start in seq.int(1L, length(cc), by = batch_size)) {
      batch_id <- batch_id + 1L
      bc <- cc[start:min(start + batch_size - 1L, length(cc))]
      mat <- SeuratObject::LayerData(object[["RNA"]], layer = data_layers[i], features = features, cells = bc)
      mat <- mat[features, bc, drop = FALSE]
      if (!inherits(mat, "dgCMatrix") || any(!is.finite(mat@x)) || any(mat@x < 0)) stop("Expected finite nonnegative existing sparse LogNormalize RNA data.")
      for (sn in names(gene_sets)) {
        # Matrix::colMeans includes structural zeros; denominator never varies by cell.
        scores[sn, bc] <- Matrix::colMeans(mat[gene_sets[[sn]], , drop = FALSE])
      }
      batches[[batch_id]] <- data.frame(batch = batch_id, layer = data_layers[i], n_cells = length(bc),
        first_cell = bc[1], last_cell = tail(bc, 1), stringsAsFactors = FALSE)
      message("Mean log-expression batch ", batch_id, ": ", length(bc), " cells")
    }
  }
  if (any(!is.finite(scores)) || any(scores < 0)) stop("Mean score output is nonfinite or negative.")
  settings <- list(method = "mean_log_expression", input = input, source_size = source_before$size,
    source_mtime = source_before$mtime, mapping_file = mapping_file, universe_file = universe_file,
    signature_md5 = signature_md5, input_md5 = tools::md5sum(c(mapping_file, universe_file,
      if (!is.null(prepared$adjudication)) prepared$adjudication_file)), adjudication = prepared$adjudication,
    universe = universe, data_common_features = length(data_common), normalization = normalization_evidence,
    formula = "sum(existing RNA LogNormalize values for fixed unique retained pig features) / number of those features; zeros included",
    denominator = lengths(gene_sets), cluster_col = cluster_col, clusters = clusters,
    selection = if (is.null(selection_note)) "all source cells" else selection_note,
    batch_size = batch_size, batch_schedule = do.call(rbind, batches),
    feature_order = "frozen radix order", cell_order = "radix data layer then radix barcode",
    source_n_cells = ncol(object), scored_n_cells = ncol(scores),
    versions = vapply(c("SeuratObject", "Matrix"), function(p) as.character(utils::packageVersion(p)), character(1)),
    session = utils::sessionInfo())
  result <- list(scores = scores, metadata = metadata, gene_sets = gene_sets, coverage = coverage,
    mapping = mapping, settings = settings)
  saveRDS(result, file.path(output_root, "signature_scores.rds"))
  write.csv(settings$batch_schedule, file.path(output_root, "batch_schedule.csv"), row.names = FALSE)
  write.csv(data.frame(cell_id = rownames(metadata), metadata, check.names = FALSE), file.path(output_root, "cell_metadata.csv"), row.names = FALSE)
  reopened <- readRDS(file.path(output_root, "signature_scores.rds"))
  stopifnot(identical(reopened, result), identical(colnames(scores), rownames(metadata)),
    identical(source_before, file.info(input)[, c("size", "mtime")]),
    identical(signature_md5, unname(tools::md5sum("R/DEG/signatures.R"))))
  cat(paste0("\nStatus: COMPLETED; exact saved-RDS reopen PASS; source size/mtime/signature hash unchanged.\n",
    "Cells: ", ncol(scores), "; arms: ", nrow(scores), "; fixed universe: ", length(features),
    "; data-layer intersection: ", length(data_common), "\n",
    "Saved NormalizeData.RNA: LogNormalize; scale.factor=", normalization@params$scale.factor,
    "; timestamp=", normalization@time.stamp, ". Existing values used unchanged.\n",
    "Mean score range: ", paste(range(scores), collapse = " to "), "; not constrained to [0,1].\n",
    "No stochastic scoring, ranks, cutoff, reference/control-gene subtraction or cell filtering.\n"), file = info_path, append = TRUE)
  completed <- TRUE
  invisible(result)
}
