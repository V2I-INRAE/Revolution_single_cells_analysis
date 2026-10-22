# Refit marker entry point: cluster markers for the refitted cleaned object,
# then a top-5 dotplot and a top-25 heatmap via the shared cluster plots.
# Usage: Rscript R/clusters/refit_markers.R <data/job-XXX/refit_clean_annotated.rds> <resolution>
# Discovery mirrors R/clusters/subset_markers.R (same settings, ranking, CSV
# schemas); object preparation mirrors the parent path with scale.data kept
# for display-scaling provenance. No per-gene feature panels.
suppressPackageStartupMessages({
  library(Seurat)
  library(presto)
})
source("R/config.R")
source("R/cell_cycle.R") # for rename_cc_features() used below
source("R/utils/save_plots.R")
source("R/plots/marker_cluster_plots.R")

run_refit_markers <- function(input, resolution) {
  if (basename(input) != "refit_clean_annotated.rds") {
    stop("Input must be a refit_clean_annotated.rds file, got: ", input)
  }
  job_id <- regmatches(normalizePath(input, mustWork = FALSE),
    regexpr("job-[0-9]+", normalizePath(input, mustWork = FALSE)))
  if (!length(job_id)) {
    stop("Input path must contain a job-<id> directory, got: ", input)
  }
  cfg <- config$markers
  sel.clust <- paste0(cfg$cluster_prefix, resolution)
  directory <- file.path("results", job_id, "clusters_plots",
    paste0("markers_res", resolution))
  if (dir.exists(directory)) {
    stop("Output directory already exists; refusing to overwrite it: ", directory)
  }
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)

  info_path <- file.path(directory, ".INFO")
  info <- c(
    "# Refit markers and plots",
    "Status: running",
    paste("Started:", Sys.time()),
    paste("Marker producer job:", Sys.getenv("SLURM_JOB_ID", unset = "manual")),
    paste("Input (read-only):", input),
    paste("Output:", directory),
    paste("Cluster column:", sel.clust),
    paste("Seed:", config$run$seed),
    paste("FindAllMarkers: RNA data; test.use=", cfg$test_use,
      "; only.pos=", cfg$only_pos, "; min.pct=", cfg$min_pct,
      "; logfc.threshold=", cfg$logfc_threshold, sep = ""),
    "min.diff.pct=-Inf; max.cells.per.ident=Inf; return.thresh=Inf.",
    paste("Significance: adjusted p <", cfg$adjusted_p_cutoff,
      "; ranking: descending avg_log2FC, adjusted p, gene."),
    "Reuse normalized RNA; join data/scale.data layers and translate existing feature IDs in memory only.",
    "Dotplot: RNA data slot; scale=FALSE (no gene-wise scaling); red colour scale, limits 0 to plot maximum.",
    "No normalization, QC, reclustering or source RDS writes. Exploratory cluster markers, not condition DE.",
    paste("Heatmap: display scaling recomputed over exactly the plotted top25 genes;",
      "at most", cluster_heatmap_max_cells_per_ident,
      "cells per cluster; RdBu palette, fixed limits",
      paste(cluster_heatmap_disp_limits, collapse = " to "), ";",
      "dotplot top5 and heatmap top25; PNG and PDF exports."),
    paste("Seurat:", packageVersion("Seurat"), "; presto:", packageVersion("presto"))
  )
  completed <- FALSE
  validation <- "Validation: pending"
  writeLines(info, info_path)
  on.exit({
    info[info == "Status: running"] <- paste("Status:",
      if (completed) "completed" else "failed/incomplete; see execution log")
    writeLines(c(info, paste("Finished:", Sys.time()), validation), info_path)
  }, add = TRUE)
  input_info <- file.info(input)[, c("size", "mtime")]

  message("Loading refit marker input: ", input, " (resolution ", resolution, ")")
  alldata <- readRDS(input)
  cells <- colnames(alldata)
  metadata <- alldata[[]]
  stopifnot(!anyDuplicated(cells), identical(rownames(metadata), cells),
    !anyNA(metadata$sample), !anyNA(metadata[[sel.clust]]))
  DefaultAssay(alldata) <- "RNA"
  # Keep normalized RNA for discovery plus existing scaled values; drop
  # counts, graphs and reductions from the working copy (no feature panels).
  layers <- Layers(alldata[["RNA"]], search = "^(data|scale\\.data)($|\\.)")
  alldata <- DietSeurat(alldata, assays = "RNA", layers = layers,
    dimreducs = NULL, graphs = NULL, misc = TRUE)
  alldata <- JoinLayers(alldata, assay = "RNA", layers = "^data($|\\.)", new = "data")
  scale_layers <- Layers(alldata[["RNA"]], search = "^scale\\.data($|\\.)")
  if (length(scale_layers) > 1L) {
    alldata <- JoinLayers(alldata, assay = "RNA",
      layers = "^scale\\.data($|\\.)", new = "scale.data")
  }
  message("RNA layers kept: ",
    paste(Layers(alldata[["RNA"]], search = "^(data|scale\\.data)($|\\.)"),
      collapse = ", "))

  renames_applied <- intersect(names(cfg$feature_renames), rownames(alldata[["RNA"]]))
  alldata <- rename_cc_features(alldata, cfg$feature_renames)
  message("Feature renames applied: ", length(renames_applied), "/", length(cfg$feature_renames))

  # Refit objects restore active identities to the cleaned annotation; select
  # the requested resolution explicitly.
  Idents(alldata) <- alldata[[sel.clust, drop = TRUE]]
  stopifnot(identical(colnames(alldata), cells),
    identical(as.character(alldata$subset), as.character(metadata$subset)))
  info <- c(info, paste("Actual:", ncol(alldata), "cells;",
    nlevels(Idents(alldata)), "clusters;", length(unique(alldata$sample)), "samples."),
    paste("Feature renames applied:", length(renames_applied), "/", length(cfg$feature_renames)),
    capture.output(print(cfg$feature_renames[renames_applied])))
  writeLines(info, info_path)

  set.seed(config$run$seed)
  markers <- FindAllMarkers(alldata, assay = "RNA", slot = "data",
    test.use = cfg$test_use, only.pos = cfg$only_pos,
    min.pct = cfg$min_pct, logfc.threshold = cfg$logfc_threshold,
    min.diff.pct = -Inf, max.cells.per.ident = Inf, return.thresh = Inf,
    verbose = FALSE)
  markers <- markers |>
    dplyr::mutate(cluster = factor(cluster, levels = levels(Idents(alldata))),
      pct_difference = pct.1 - pct.2) |>
    dplyr::arrange(cluster, dplyr::desc(avg_log2FC), p_val_adj, gene)
  significant <- dplyr::filter(markers, p_val_adj < cfg$adjusted_p_cutoff)
  top5 <- significant |> dplyr::group_by(cluster) |> dplyr::slice_head(n = 5) |> dplyr::ungroup()
  top25 <- significant |> dplyr::group_by(cluster) |> dplyr::slice_head(n = 25) |> dplyr::ungroup()

  write.csv(markers, file.path(directory, "markers_prefiltered.csv"), row.names = FALSE)
  write.csv(significant, file.path(directory, "markers_significant.csv"), row.names = FALSE)
  saved_tables <- list(markers_prefiltered.csv = markers,
    markers_significant.csv = significant)
  for (n in c(5, 20, 25)) {
    top <- significant |> dplyr::group_by(cluster) |> dplyr::slice_head(n = n) |> dplyr::ungroup()
    write.csv(top, file.path(directory, paste0("markers_top", n, ".csv")), row.names = FALSE)
    saved_tables[[paste0("markers_top", n, ".csv")]] <- top
  }
  counts <- table(Idents(alldata))
  summary <- data.frame(cluster = names(counts), n_cells = as.integer(counts),
    n_prefiltered = as.integer(table(markers$cluster)),
    n_significant = as.integer(table(significant$cluster)),
    n_adjusted_p_zero = as.integer(table(markers$cluster[markers$p_val_adj == 0])))
  composition <- as.data.frame(table(sample = alldata$sample,
    cluster = alldata[[sel.clust, drop = TRUE]]), responseName = "n_cells")
  write.csv(summary, file.path(directory, "cluster_summary.csv"), row.names = FALSE)
  write.csv(composition, file.path(directory, "sample_cluster_counts.csv"), row.names = FALSE)
  saved_tables$cluster_summary.csv <- summary
  saved_tables$sample_cluster_counts.csv <- composition

  pipeline_csvs <- c("markers_prefiltered.csv", "markers_significant.csv",
    "markers_top5.csv", "markers_top20.csv", "markers_top25.csv",
    "cluster_summary.csv", "sample_cluster_counts.csv")
  for (csv in pipeline_csvs) {
    reopened <- read.csv(file.path(directory, csv), stringsAsFactors = FALSE, check.names = FALSE)
    expected <- as.data.frame(saved_tables[[csv]])
    for (column in names(expected)) {
      if (is.factor(expected[[column]]) || is.character(expected[[column]])) {
        expected[[column]] <- as.character(expected[[column]])
        reopened[[column]] <- as.character(reopened[[column]])
      }
    }
    stopifnot(identical(names(reopened), names(expected)),
      isTRUE(all.equal(as.list(reopened), as.list(expected),
        check.attributes = FALSE, tolerance = 1e-12)))
  }
  stopifnot(sum(summary$n_cells) == ncol(alldata), sum(composition$n_cells) == ncol(alldata))
  rm(saved_tables, reopened, expected)

  heatmap_genes <- unique(as.character(top25$gene))
  hvg_overlap <- intersect(heatmap_genes,
    rownames(GetAssayData(alldata, assay = "RNA", layer = "scale.data")))
  alldata <- Seurat::ScaleData(alldata, assay = "RNA",
    features = heatmap_genes,
    vars.to.regress = config$norm$vars_to_regress,
    verbose = FALSE)
  info <- c(info, paste(
    "Heatmap display scaling: ScaleData over exactly the plotted top25 genes",
    paste0("(", length(heatmap_genes), " genes),"),
    "same regression covariates as the pipeline; stored HVG-only scale.data untouched.",
    paste0("Top25 genes already scaled: ", length(hvg_overlap),
      "/", length(heatmap_genes), ".")))
  writeLines(info, info_path)

  plot_cluster_dotplot(alldata, top5, sel.clust, directory,
    level_label = "Cluster")
  plot_cluster_heatmap(alldata, heatmap_genes, sel.clust, directory,
    seed = config$run$seed, level_label = "Cluster")
  figures <- file.path(directory, c("dotplot_top5.png", "dotplot_top5.pdf",
    "heatmap_top25.png", "heatmap_top25.pdf"))
  stopifnot(all(file.exists(figures)), all(file.size(figures) > 0))
  for (path in figures[grepl("\\.png$", figures)]) {
    image <- png::readPNG(path, native = TRUE)
    stopifnot(all(dim(image) > 0))
    rm(image)
  }
  for (path in figures[grepl("\\.pdf$", figures)]) {
    stopifnot(identical(rawToChar(readBin(path, what = "raw", n = 5L)), "%PDF-"))
  }
  stopifnot(identical(file.info(input)[, c("size", "mtime")], input_info))
  validation <- c(
    "Validation: all seven CSVs reopened and values matched; cluster/sample counts conserve cells; PNGs decoded; nonempty PDF signatures checked; source RDS size/mtime unchanged.",
    "PDF page rendering and visual/biological review remain pending; no PDF parser is installed.",
    paste("Marker rows:", nrow(markers), "; significant rows:", nrow(significant)),
    capture.output(print(counts))
  )
  completed <- TRUE
  message("Refit marker tables and plots completed: ", directory)
  invisible(directory)
}

# Sourcing defines the runner without executing marker discovery.
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) {
    stop("Usage: Rscript R/clusters/refit_markers.R <data/job-XXX/refit_clean_annotated.rds> <resolution>")
  }
  run_refit_markers(args[[1L]], as.numeric(args[[2L]]))
}
