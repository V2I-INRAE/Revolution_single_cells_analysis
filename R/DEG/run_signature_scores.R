# Source from the project root. This workflow is separate from pseudobulk DE.
# Example (explicit NEW directory; full scoring is not run by sourcing):
# run_signature_scores(
#   input = "data/job-44542861/subsets/Endothelial/job-44831031/subclustered.rds",
#   mapping_file = "results/misc/Ensembl_ID/signature-scoring-v2/pig_mapping.csv",
#   universe_file = "results/misc/Ensembl_ID/signature-scoring/common_count_universe.rds",
#   output_root = "results/Endothelial_cluster1_signature_scores_new",
#   cluster_col = "harmony_snn_res.0.3", clusters = "1")
# The full original subset count layers are validated BEFORE selecting all cells
# of the requested clusters. Omit cluster_col/clusters to score all source cells.
# For a recorded pilot only, use cells + selection_note instead of cluster selection.
# Changing selected cells/batches can change randomized ties. Replot saved scores;
# do not rerank for each condition. Cross-run comparisons require the same protocol.
source("R/DEG/signatures.R")
source("R/DEG/signature_gene_sets.R")

run_signature_scores <- function(input, mapping_file, universe_file, output_root,
    seed = 1234L, batch_size = 1000L, cells = NULL, selection_note = NULL,
    cluster_col = NULL, clusters = NULL) {
  if (file.exists(output_root)) stop("output_root exists; use a new directory.")
  for (p in c("SeuratObject", "Matrix", "AUCell", "BiocParallel", "DelayedArray")) {
    if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
  }
  if (batch_size < 1L || batch_size != as.integer(batch_size) ||
      seed < 0 || seed != as.integer(seed)) stop("Use a positive integer batch_size and nonnegative integer seed.")
  if (xor(is.null(cluster_col), is.null(clusters))) stop("Supply both cluster_col and clusters.")
  if (!is.null(cluster_col) && !is.null(cells)) stop("Cluster selection scores ALL selected cluster cells; do not also supply cells.")
  if (!is.null(cluster_col)) {
    clusters <- as.character(clusters)
    if (length(cluster_col) != 1L || !nzchar(cluster_col) || !length(clusters) ||
        anyNA(clusters) || anyDuplicated(clusters) || any(!nzchar(clusters))) stop("Use one cluster column and unique nonempty cluster labels.")
    selection_note <- paste0("All cells in ", cluster_col, " clusters ", paste(clusters, collapse = ","),
      "; frozen universe verified from all original subset count layers before selection; no cell/sample subsampling.")
  }
  if (!is.null(cells) && (is.null(selection_note) || !nzchar(selection_note))) {
    stop("Explicit cell selection requires selection_note; never silently downsample.")
  }
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  info_path <- file.path(output_root, ".INFO")
  writeLines(c("# AUCell signature scoring", "Status: RUNNING", paste("Started:", Sys.time()),
    paste("Input (read-only):", input), paste("Input producer:", basename(dirname(input))),
    paste("Mapping:", mapping_file), paste("Fixed count universe:", universe_file),
    paste("Seed:", seed, "; batch_size:", batch_size),
    paste("Selection:", if (is.null(selection_note)) "all source cells" else selection_note),
    "Original RNA counts; no normalization, integration, scaling, cell filtering or object write-back.",
    "Fixed common-layer universe; top5%; separate signed arms; all nonempty approved sets.",
    "Random ties, including zeros: installed AUCell default keepZeroesAsNA=FALSE; serial backend.",
    "Source lists unchanged; derived mapping excludes unresolved cases; unique pig features per arm."), info_path)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  completed <- FALSE
  on.exit(if (!completed) cat("\nStatus: FAILED; partial outputs are not validated scores.\n", file = info_path, append = TRUE))
  source_before <- file.info(input)[, c("size", "mtime")]
  signature_md5 <- unname(tools::md5sum("R/DEG/signatures.R"))
  universe <- readRDS(universe_file)
  features <- universe$features
  if (!length(features) || anyDuplicated(features)) stop("Invalid frozen feature universe.")
  prepared <- prepare_signature_gene_sets(mapping_file, features, signatures)
  mapping <- prepared$mapping; gene_sets <- prepared$gene_sets; coverage <- prepared$coverage
  adjudication <- prepared$adjudication; adjudication_file <- prepared$adjudication_file
  # <10 is a descriptive warning, never an eligibility threshold.
  write.csv(mapping, file.path(output_root, "mapping_audit.csv"), row.names = FALSE)
  write.csv(coverage, file.path(output_root, "coverage.csv"), row.names = FALSE)
  gene_sets <- gene_sets[lengths(gene_sets) > 0L]
  if (!length(gene_sets)) stop("No nonempty approved gene sets in frozen universe.")
  if (any(coverage$human_coverage < 1)) warning("Incomplete source-human coverage; see coverage.csv.")
  if (any(coverage$small_set_warning)) warning("Small effective sets (<10 genes); interpret cautiously; see coverage.csv.")
  object <- readRDS(input)
  metadata <- object[[]]
  required <- c("sample", "pig", "pressure", "time_point")
  if (!all(required %in% names(metadata)) || !identical(rownames(metadata), colnames(object)) ||
      anyDuplicated(rownames(metadata)) || anyNA(metadata[, required, drop = FALSE]) ||
      any(vapply(metadata[, required, drop = FALSE], function(x) any(!nzchar(as.character(x))), logical(1)))) {
    stop("Required metadata missing or cell IDs/metadata not uniquely aligned.")
  }
  panel <- unique(metadata[, required, drop = FALSE])
  if (anyDuplicated(panel$sample)) stop("Each sample must have one pig/pressure/time tuple.")
  if (!setequal(as.character(panel$sample), universe$samples)) stop("Input must retain the frozen source sample panel; select cells only after validation.")
  layers <- SeuratObject::Layers(object[["RNA"]], search = "^counts[.]")
  if (!length(layers)) stop("Original sample-specific RNA counts.* layers required.")
  layers <- sort(layers, method = "radix")
  layer_cells <- lapply(layers, function(l) SeuratObject::Cells(object[["RNA"]], layer = l))
  layer_features <- lapply(layers, function(l) SeuratObject::Features(object[["RNA"]], layer = l))
  all_layer_cells <- unlist(layer_cells, use.names = FALSE)
  if (anyDuplicated(all_layer_cells) || !setequal(all_layer_cells, rownames(metadata))) stop("Count layers must partition source cells exactly.")
  for (i in seq_along(layers)) {
    if (!identical(unique(as.character(metadata[layer_cells[[i]], "sample"])), sub("^counts[.]", "", layers[i]))) {
      stop("Sample/count-layer identity mismatch: ", layers[i])
    }
  }
  if (!setequal(sub("^counts[.]", "", layers), universe$samples) ||
      !setequal(Reduce(intersect, layer_features), features)) {
    stop("Input common count-layer universe differs from frozen universe; no zero padding or automatic universe change.")
  }
  # Selection occurs only after full source panel/layer/universe validation above.
  if (!is.null(cluster_col)) {
    if (!cluster_col %in% names(metadata) || anyNA(metadata[[cluster_col]]) ||
        any(!nzchar(as.character(metadata[[cluster_col]]))) ||
        !all(clusters %in% as.character(metadata[[cluster_col]]))) stop("Requested cluster labels/column are absent or incomplete in source metadata.")
    cells <- rownames(metadata)[as.character(metadata[[cluster_col]]) %in% clusters]
  }
  if (is.null(cells)) cells <- rownames(metadata)
  if (!length(cells) || anyDuplicated(cells) || !all(cells %in% rownames(metadata))) stop("Selected cells must be unique source cell IDs.")
  # Fixed sample layer order, then radix-sorted barcodes; batch boundaries never cross layers.
  selected <- lapply(layer_cells, function(z) sort(z[z %in% cells], method = "radix"))
  ordered_cells <- unlist(selected, use.names = FALSE)
  metadata <- metadata[ordered_cells, , drop = FALSE]
  auc_max_rank <- ceiling(0.05 * length(features))
  scores <- matrix(NA_real_, nrow = length(gene_sets), ncol = length(ordered_cells),
    dimnames = list(names(gene_sets), ordered_cells))
  depth <- setNames(integer(length(ordered_cells)), ordered_cells)
  batch_rows <- list()
  batch_id <- 0L
  for (i in seq_along(layers)) {
    cc <- selected[[i]]
    if (!length(cc)) next
    for (start in seq.int(1L, length(cc), by = batch_size)) {
      batch_id <- batch_id + 1L
      bc <- cc[start:min(start + batch_size - 1L, length(cc))]
      batch_seed <- seed + batch_id - 1L
      mat <- SeuratObject::LayerData(object[["RNA"]], layer = layers[i], features = features, cells = bc)
      mat <- mat[features, bc, drop = FALSE]
      if (!inherits(mat, "dgCMatrix") || any(!is.finite(mat@x)) || any(mat@x < 0)) stop("Expected finite nonnegative original sparse RNA counts.")
      depth[bc] <- Matrix::colSums(mat > 0)
      set.seed(batch_seed)
      rankings <- AUCell::AUCell_buildRankings(mat, plotStats = FALSE, splitByBlocks = TRUE,
        BPPARAM = BiocParallel::SerialParam(RNGseed = batch_seed), keepZeroesAsNA = FALSE, verbose = FALSE)
      auc <- AUCell::getAUC(AUCell::AUCell_calcAUC(gene_sets, rankings,
        nCores = 1L, normAUC = TRUE, aucMaxRank = auc_max_rank, verbose = FALSE))
      if (!identical(rownames(auc), names(gene_sets)) || !identical(colnames(auc), bc) ||
          any(!is.finite(auc)) || any(auc < 0 | auc > 1)) stop("AUCell output IDs or normalized values invalid.")
      scores[, bc] <- auc
      batch_rows[[batch_id]] <- data.frame(batch = batch_id, layer = layers[i], seed = batch_seed,
        n_cells = length(bc), first_cell = bc[1], last_cell = tail(bc, 1), stringsAsFactors = FALSE)
      message("AUCell batch ", batch_id, ": ", length(bc), " cells")
      rm(mat, rankings, auc)
    }
  }
  metadata$signature_detected_features <- unname(depth[rownames(metadata)])
  metadata$signature_below_auc_cutoff <- metadata$signature_detected_features < auc_max_rank
  versions <- vapply(c("AUCell", "SeuratObject", "Matrix", "BiocParallel", "DelayedArray", "DelayedMatrixStats"),
    function(p) as.character(utils::packageVersion(p)), character(1))
  settings <- list(method = "AUCell", input = input, source_size = source_before$size, source_mtime = source_before$mtime,
    mapping_file = mapping_file, universe_file = universe_file, signature_md5 = signature_md5,
    input_md5 = tools::md5sum(c(mapping_file, universe_file,
      if (!is.null(adjudication)) adjudication_file)), adjudication = adjudication, universe = universe,
    auc_max_rank = auc_max_rank, normAUC = TRUE, keepZeroesAsNA = FALSE,
    seed = seed, batch_size = batch_size, backend = "BiocParallel::SerialParam(RNGseed=batch_seed)",
    feature_order = "frozen radix order", cell_order = "radix layer then radix barcode",
    rng_kind = RNGkind(), delayed_array_block_size = DelayedArray::getAutoBlockSize(),
    versions = versions, selection = if (is.null(selection_note)) "all source cells" else selection_note,
    cluster_col = cluster_col, clusters = clusters,
    source_n_cells = ncol(object), scored_n_cells = ncol(scores),
    batch_schedule = do.call(rbind, batch_rows), session = utils::sessionInfo())
  result <- list(scores = scores, metadata = metadata, gene_sets = gene_sets,
    coverage = coverage, mapping = mapping, settings = settings)
  saveRDS(result, file.path(output_root, "signature_scores.rds"))
  write.csv(settings$batch_schedule, file.path(output_root, "batch_schedule.csv"), row.names = FALSE)
  write.csv(data.frame(cell_id = rownames(metadata), metadata, check.names = FALSE),
    file.path(output_root, "cell_metadata.csv"), row.names = FALSE)
  reopened <- readRDS(file.path(output_root, "signature_scores.rds"))
  stopifnot(identical(reopened, result), identical(colnames(scores), rownames(metadata)),
    identical(source_before, file.info(input)[, c("size", "mtime")]),
    identical(signature_md5, unname(tools::md5sum("R/DEG/signatures.R"))))
  cat(paste0("\nStatus: COMPLETED; saved RDS reopened and exact equality PASS.\n",
    "Source size/mtime and signature hash unchanged.\n",
    "Scored cells: ", ncol(scores), "; gene-set arms: ", nrow(scores),
    "; features: ", length(features), "; aucMaxRank: ", auc_max_rank, "\n",
    "Cells below cutoff (retained): ", sum(metadata$signature_below_auc_cutoff), "\n",
    "AUCell: ", versions[["AUCell"]], "; recorded seed/order/backend/batch plan required for reproduction.\n",
    "No activity thresholds, statistical tests, composite signed score or cell filtering.\n"),
    file = info_path, append = TRUE)
  completed <- TRUE
  invisible(result)
}
