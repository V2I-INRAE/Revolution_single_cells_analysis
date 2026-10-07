args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript scripts/soupx/run_soupx.R <sample>")
sample_id <- args[1]
source("scripts/soupx/io.R")
source("scripts/soupx/plots.R")
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(SoupX)
})
future::plan("sequential")

job_id <- Sys.getenv("SLURM_JOB_ID")
array_id <- Sys.getenv("SLURM_ARRAY_JOB_ID")
run_id <- if (nzchar(array_id)) paste0("job-", array_id) else if (nzchar(job_id)) paste0("job-", job_id) else
  paste(format(Sys.time(), "%Y%m%d-%H%M%S"), Sys.getpid(), sep = "-")
data_dir <- file.path("data", "soupx", sample_id, run_id)
results_dir <- file.path("results", "soupx", run_id, sample_id)
for (directory in c(data_dir, results_dir)) {
  dir.create(dirname(directory), recursive = TRUE, showWarnings = FALSE)
  if (!dir.create(directory, showWarnings = FALSE)) {
    stop("Cannot create new run directory: ", directory)
  }
}

settings <- list(seed = 1234L, normalization = "LogNormalize", scale_factor = 10000,
  hvg_method = "vst", n_hvgs = 2000L, computed_pcs = 30L, neighbor_pcs = 1:20,
  k = 20L, algorithm = "Louvain", resolution = 0.6, regression = "none",
  integration = "none", soup_range = c(0, 100), knee_min_exclusive = 1000,
  min_noncalled_background = 2000L, force_accept = FALSE,
  adjustment = "subtraction", round_to_int = TRUE)
input_dir <- file.path("data", "raw_data", sample_id)
inputs <- c(
  unfiltered = file.path(input_dir, paste0(sample_id, "_RSEC_MolsPerCell_Unfiltered_MEX.zip")),
  filtered = file.path(input_dir, paste0(sample_id, "_RSEC_MolsPerCell_MEX.zip")),
  metrics = file.path(input_dir, paste0(sample_id, "_Metrics_Summary.csv")))
started <- Sys.time()
provenance <- list(sample = sample_id, run_id = run_id, status = "running",
  slurm_job_id = job_id, slurm_array_job_id = array_id,
  slurm_array_task_id = Sys.getenv("SLURM_ARRAY_TASK_ID"),
  started = format(started, tz = "UTC", usetz = TRUE), settings = settings,
  inputs = as.list(inputs), warnings = character(),
  versions = setNames(lapply(c("Seurat", "SoupX", "Matrix", "jsonlite", "scplotter"),
    function(package) as.character(packageVersion(package))),
    c("Seurat", "SoupX", "Matrix", "jsonlite", "scplotter")), rng_kind = RNGkind(),
  scripts_md5 = as.list(tools::md5sum(c("scripts/soupx/io.R",
    "scripts/soupx/plots.R", "scripts/soupx/run_soupx.R", "scripts/sbatch_soupx.sh"))))
write_provenance <- function() {
  jsonlite::write_json(provenance, file.path(results_dir, "run.json"),
    pretty = TRUE, auto_unbox = TRUE, null = "null")
}
write_provenance()
writeLines(capture.output(sessionInfo()), file.path(results_dir, "sessionInfo.txt"))

run <- function() {
  stopifnot(all(file.exists(inputs)))
  provenance$inputs <<- lapply(inputs, function(path) list(
    path = normalizePath(path), bytes = file.info(path)$size,
    md5 = unname(tools::md5sum(path))))
  manifest <- Sys.getenv("SOUPX_SAMPLE_MANIFEST")
  if (nzchar(manifest)) {
    provenance$sample_manifest <<- list(path = manifest,
      md5 = unname(tools::md5sum(manifest)))
  }
  message("Reading and aligning raw RSEC counts for ", sample_id)
  extraction <- tempfile(paste0("soupx-", sample_id, "-"))
  dir.create(extraction)
  on.exit(unlink(extraction, recursive = TRUE), add = TRUE)
  raw <- read_soupx_mex(inputs[["unfiltered"]], file.path(extraction, "raw"))
  filtered <- read_soupx_mex(inputs[["filtered"]], file.path(extraction, "filtered"))
  aligned <- align_soupx_counts(raw, filtered)
  expected <- read_putative_cell_count(inputs[["metrics"]])
  stopifnot(ncol(aligned$toc) == expected)
  provenance$input_validation <<- list(features = nrow(aligned$tod),
    barcodes = ncol(aligned$tod), called_cells = ncol(aligned$toc),
    zero_padded_features = aligned$zero_padded_features,
    repaired_symbols = sum(aligned$feature_map$symbol != aligned$feature_map$soupx_name),
    raw_called_counts_equal_filtered = TRUE)
  write.table(aligned$feature_map, file.path(data_dir, "feature_map.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
  background <- soupx_background_check(aligned$tod, aligned$toc, expected)
  provenance$background <<- background
  write_provenance()
  jsonlite::write_json(background, file.path(results_dir, "background_check.json"),
    pretty = TRUE, auto_unbox = TRUE)
  if (!background$passed) stop(
    "Background gate failed: inspect background_check.json; no automatic range change.")
  message("Background gate passed: knee proxy = ", background$molecules_at_called_cell_rank,
    "; non-called background barcodes = ", background$noncalled_barcodes_in_range)
  rm(raw, filtered)
  gc()

  message("Clustering all ", ncol(aligned$toc), " called cells before correction")
  set.seed(settings$seed)
  obj <- CreateSeuratObject(aligned$toc, project = sample_id,
    min.cells = 0, min.features = 0)
  stopifnot(identical(colnames(obj), colnames(aligned$toc)))
  obj <- NormalizeData(obj, normalization.method = settings$normalization,
    scale.factor = settings$scale_factor)
  obj <- FindVariableFeatures(obj, selection.method = settings$hvg_method,
    nfeatures = settings$n_hvgs)
  obj <- ScaleData(obj, features = VariableFeatures(obj))
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = settings$computed_pcs,
    seed.use = settings$seed)
  obj <- FindNeighbors(obj, dims = settings$neighbor_pcs, k.param = settings$k)
  obj <- FindClusters(obj, resolution = settings$resolution, algorithm = 1,
    random.seed = settings$seed)
  clusters <- setNames(as.character(Idents(obj)), colnames(obj))
  stopifnot(identical(names(clusters), colnames(aligned$toc)), !anyNA(clusters))
  saveRDS(list(clusters = clusters, pca = Embeddings(obj, "pca"),
    variable_features = VariableFeatures(obj), settings = settings,
    commands = obj@commands), file.path(data_dir, "preliminary_clustering.rds"))
  write.csv(data.frame(cell = names(clusters), cluster = clusters),
    file.path(results_dir, "cluster_assignments.csv"), row.names = FALSE)
  provenance$cluster_sizes <<- as.list(table(clusters))
  rm(obj)
  gc()

  # SoupChannel estimates the soup profile once using its default range.
  message("Estimating contamination with official SoupX ", packageVersion("SoupX"))
  sc <- SoupChannel(tod = aligned$tod, toc = aligned$toc)
  sc <- setClusters(sc, clusters)
  saveRDS(sc, file.path(data_dir, "soupx_before_estimation.rds"), compress = FALSE)
  provenance$auto_estimation_defaults <<- lapply(formals(autoEstCont)[-1], deparse)
  png(file.path(results_dir, "contamination_estimation.png"),
    width = 2400, height = 1800, res = 300)
  sc <- tryCatch(autoEstCont(sc, forceAccept = FALSE), finally = dev.off())
  stopifnot(length(sc$fit$rhoEst) == 1L, is.finite(sc$fit$rhoEst),
    sum(sc$fit$dd$useEst) > 0, all(is.finite(sc$fit$posterior)))
  write.csv(sc$soupProfile, file.path(results_dir, "soup_profile.csv"))
  write.csv(sc$fit$dd, file.path(results_dir, "contamination_estimates.csv"),
    row.names = FALSE)
  write.csv(sc$fit$markersUsed, file.path(results_dir, "estimation_markers.csv"),
    row.names = FALSE)
  provenance$estimation <<- list(rho = sc$fit$rhoEst, rho_fwhm = sc$fit$rhoFWHM,
    independent_estimates = sum(sc$fit$dd$useEst),
    candidate_marker_count = nrow(sc$fit$markersUsed),
    marker_count = length(unique(sc$fit$dd$gene)))
  aligned$tod <- NULL
  gc()
  set.seed(settings$seed)
  corrected <- adjustCounts(sc, method = settings$adjustment,
    roundToInt = settings$round_to_int)
  validate_soupx_counts(corrected, sc$toc)
  saveRDS(corrected, file.path(data_dir, "corrected_counts.rds"), compress = FALSE)
  saveRDS(sc, file.path(data_dir, "soupx_channel.rds"), compress = FALSE)
  # Reopen both count sources: a successful write is not serialization validation.
  reopened <- readRDS(file.path(data_dir, "corrected_counts.rds"))
  channel <- readRDS(file.path(data_dir, "soupx_channel.rds"))
  validate_soupx_counts(reopened, channel$toc)
  stopifnot(identical(reopened, corrected), identical(channel$toc, sc$toc))
  rm(reopened, channel)

  raw_total <- Matrix::colSums(sc$toc)
  corrected_total <- Matrix::colSums(corrected)
  per_cell <- data.frame(cell = colnames(corrected), cluster = clusters,
    raw_molecules = raw_total, corrected_molecules = corrected_total,
    removed_fraction = 1 - corrected_total / raw_total,
    raw_genes = Matrix::colSums(sc$toc > 0),
    corrected_genes = Matrix::colSums(corrected > 0))
  connection <- gzfile(file.path(results_dir, "cell_counts.csv.gz"), "wt")
  write.csv(per_cell, connection, row.names = FALSE)
  close(connection)
  # Group-level expression is exported for review, without declaring cluster cell types.
  group <- Matrix::sparseMatrix(i = seq_along(clusters),
    j = match(clusters, unique(clusters)), x = 1,
    dims = c(length(clusters), length(unique(clusters))))
  raw_group <- sc$toc %*% group
  corrected_group <- corrected %*% group
  cluster_genes <- data.frame(
    gene = rep(rownames(sc$toc), times = ncol(group)),
    cluster = rep(unique(clusters), each = nrow(sc$toc)),
    cells = rep(as.numeric(table(factor(clusters, levels = unique(clusters)))),
      each = nrow(sc$toc)),
    raw_molecules = as.vector(as.matrix(raw_group)),
    corrected_molecules = as.vector(as.matrix(corrected_group)))
  connection <- gzfile(file.path(results_dir, "gene_counts_by_cluster.csv.gz"), "wt")
  write.csv(cluster_genes, connection, row.names = FALSE)
  close(connection)
  provenance$validation <<- list(status = "passed", saved_counts_reopened = TRUE,
    called_cells = ncol(corrected), features = nrow(corrected),
    raw_molecules = sum(sc$toc), corrected_molecules = sum(corrected),
    removed_fraction = 1 - sum(corrected) / sum(sc$toc),
    nonnegative_integer_counts = TRUE, no_added_counts = TRUE)
  provenance$marker_plots <<- plot_soupx_markers(
    list(Raw = sc$toc, SoupX = corrected), clusters, sample_id, results_dir,
    labels = c("Raw", paste0("SoupX (automatic rho = ",
      format(100 * sc$fit$rhoEst, digits = 3), "%)")))
  provenance$status <<- "validated_pending_qc"
  provenance$qc_note <<- paste(
    "Numerical validation passed; diagnostic and marker preservation review required.",
    "No downstream checkpoint replaced; no CellBender comparator accepted here.")
  # This provisional checkpoint is redundant once the fitted channel is verified.
  unlink(file.path(data_dir, "soupx_before_estimation.rds"))
}

tryCatch(withCallingHandlers(run(), warning = function(w) {
  provenance$warnings <<- c(provenance$warnings, conditionMessage(w))
}), error = function(e) {
  provenance$status <<- "failed_requires_review"
  provenance$error <<- conditionMessage(e)
  provenance$finished <<- format(Sys.time(), tz = "UTC", usetz = TRUE)
  provenance$elapsed_seconds <<- as.numeric(difftime(Sys.time(), started, units = "secs"))
  write_provenance()
  stop(e)
})
provenance$finished <- format(Sys.time(), tz = "UTC", usetz = TRUE)
provenance$elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
write_provenance()
message("SoupX numerical validation complete; status: ", provenance$status)
