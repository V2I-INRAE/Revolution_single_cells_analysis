# Subset marker entry point: cluster markers for one subclustered subset,
# with pig gene-name annotation. Publication plots are generated explicitly.
# Usage: Rscript R/clusters/subset_markers.R <data/job-XXX/subsets/<Folder>/<subrun>/subclustered.rds> <resolution>
suppressPackageStartupMessages({
  library(Seurat)
  library(presto)
})
source("R/config.R")
source("R/cell_cycle.R") # for rename_cc_features() used below
source("R/utils/pig_gene_mapping.R")
source("R/plots/publications/publication_paths.R")

run_subset_markers <- function(input, resolution) {
  if (basename(input) != "subclustered.rds") {
    stop("Input must be a subclustered.rds file, got: ", input)
  }
  parts <- strsplit(
    normalizePath(input, mustWork = FALSE),
    .Platform$file.sep
  )[[1L]]
  subsets_idx <- match("subsets", parts)
  if (is.na(subsets_idx) || length(parts) < subsets_idx + 3L) {
    stop(
      "Input path must look like data/job-XXX/subsets/<Folder>/<subrun>/subclustered.rds, got: ",
      input
    )
  }
  parent_job <- parts[subsets_idx - 1L]
  folder <- parts[subsets_idx + 1L]
  subrun <- parts[subsets_idx + 2L]
  cfg <- config$markers
  sel.clust <- paste0(cfg$cluster_prefix, resolution)
  directory <- file.path(
    subset_resolution_dir(input, as.character(resolution)),
    "Markers"
  )
  if (dir.exists(directory)) {
    stop("Marker results already exist; refusing to overwrite: ", directory)
  }
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)

  info_path <- file.path(directory, ".INFO")
  info <- c(
    "# Subset marker tables",
    "Status: running",
    paste("Started:", Sys.time()),
    paste("Marker producer job:", Sys.getenv("SLURM_JOB_ID", unset = "manual")),
    paste("Input (read-only):", input),
    paste("Subclustering producer:", subrun, "; pipeline parent:", parent_job),
    paste(
      "Subclustering record:",
      file.path(
        "results",
        parent_job,
        "clusters_plots",
        "subset",
        folder,
        subrun,
        ".INFO"
      )
    ),
    paste("Output:", directory),
    paste("Cluster column:", sel.clust),
    paste("Seed:", config$run$seed),
    paste(
      "FindAllMarkers: RNA data; test.use=",
      cfg$test_use,
      "; only.pos=",
      cfg$only_pos,
      "; min.pct=",
      cfg$min_pct,
      "; logfc.threshold=",
      cfg$logfc_threshold,
      sep = ""
    ),
    "min.diff.pct=-Inf; max.cells.per.ident=Inf; return.thresh=Inf.",
    paste(
      "Significance: adjusted p <",
      cfg$adjusted_p_cutoff,
      "; ranking: descending avg_log2FC, adjusted p, gene."
    ),
    "Reuse normalized/scaled RNA; join layers and translate existing feature IDs in memory only.",
    "Gene annotation: unflagged pig Ensembl symbols and unambiguous frozen MSigDB pig rescue; feature keys retained; no human renaming.",
    "No normalization, QC, reclustering or source RDS writes. Exploratory cluster markers, not condition DE.",
    "No exploratory dotplot or heatmap images exported; use the publication scripts explicitly.",
    paste(
      "Seurat:",
      packageVersion("Seurat"),
      "; presto:",
      packageVersion("presto")
    )
  )
  completed <- FALSE
  validation <- "Validation: pending"
  writeLines(info, info_path)
  on.exit(
    {
      info[info == "Status: running"] <- paste(
        "Status:",
        if (completed) "completed" else "failed/incomplete; see execution log"
      )
      writeLines(c(info, paste("Finished:", Sys.time()), validation), info_path)
    },
    add = TRUE
  )
  input_info <- file.info(input)[, c("size", "mtime")]

  # These annotated tables are reopened after writing.
  pipeline_csvs <- c(
    "markers_prefiltered.csv",
    "markers_significant.csv",
    "markers_top5.csv",
    "markers_top20.csv",
    "markers_top25.csv",
    "cluster_summary.csv",
    "sample_cluster_counts.csv"
  )

  # Retain normalized/scaled expression for discovery and heatmap plotting.
  message(
    "Loading subset marker input: ",
    input,
    " (resolution ",
    resolution,
    ")"
  )
  alldata <- readRDS(input)
  cells <- colnames(alldata)
  metadata <- alldata[[]]
  stopifnot(
    !anyDuplicated(cells),
    identical(rownames(metadata), cells),
    !anyNA(metadata$sample),
    !anyNA(metadata[[sel.clust]])
  )
  DefaultAssay(alldata) <- "RNA"
  alldata <- JoinLayers(
    alldata,
    assay = "RNA",
    layers = "^data($|\\.)",
    new = "data"
  )
  scale_layers <- Layers(alldata[["RNA"]], search = "^scale\\.data($|\\.)")
  if (length(scale_layers) > 1L) {
    alldata <- JoinLayers(
      alldata,
      assay = "RNA",
      layers = "^scale\\.data($|\\.)",
      new = "scale.data"
    )
  }
  message(
    "RNA layers kept: ",
    paste(
      Layers(alldata[["RNA"]], search = "^(data|scale\\.data)($|\\.)"),
      collapse = ", "
    )
  )

  # Translate IDs before discovery so tables and plots use the same names.
  renames_applied <- intersect(
    names(cfg$feature_renames),
    rownames(alldata[["RNA"]])
  )
  alldata <- rename_cc_features(alldata, cfg$feature_renames)
  message(
    "Feature renames applied: ",
    length(renames_applied),
    "/",
    length(cfg$feature_renames)
  )
  if (!sel.clust %in% colnames(alldata[[]])) {
    stop("Cluster column not found in object: ", sel.clust)
  }
  Idents(alldata) <- alldata[[sel.clust, drop = TRUE]]
  stopifnot(
    identical(colnames(alldata), cells),
    identical(alldata[[]], metadata)
  )
  info <- c(
    info,
    paste(
      "Actual:",
      ncol(alldata),
      "cells;",
      nlevels(Idents(alldata)),
      "clusters;",
      length(unique(alldata$sample)),
      "samples."
    ),
    paste(
      "Feature renames applied:",
      length(renames_applied),
      "/",
      length(cfg$feature_renames)
    ),
    capture.output(print(cfg$feature_renames[renames_applied]))
  )
  writeLines(info, info_path)

  # Discovery uses log-normalized expression, not regression residuals.
  set.seed(config$run$seed)
  markers <- FindAllMarkers(
    alldata,
    assay = "RNA",
    slot = "data",
    test.use = cfg$test_use,
    only.pos = cfg$only_pos,
    min.pct = cfg$min_pct,
    logfc.threshold = cfg$logfc_threshold,
    min.diff.pct = -Inf,
    max.cells.per.ident = Inf,
    return.thresh = Inf,
    verbose = FALSE
  )
  markers <- markers |>
    dplyr::mutate(
      cluster = factor(cluster, levels = levels(Idents(alldata))),
      pct_difference = pct.1 - pct.2
    ) |>
    dplyr::arrange(cluster, dplyr::desc(avg_log2FC), p_val_adj, gene)
  gene_mapping <- prepare_pig_gene_mapping(rownames(alldata[["RNA"]]))
  markers <- annotate_pig_genes(markers, gene_mapping)
  significant <- dplyr::filter(markers, p_val_adj < cfg$adjusted_p_cutoff)

  # Save annotation alongside the unchanged marker statistics.
  write.csv(
    markers,
    file.path(directory, "markers_prefiltered.csv"),
    row.names = FALSE
  )
  write.csv(
    significant,
    file.path(directory, "markers_significant.csv"),
    row.names = FALSE
  )
  saved_tables <- list(
    markers_prefiltered.csv = markers,
    markers_significant.csv = significant
  )
  for (n in c(5, 20, 25)) {
    top <- significant |>
      dplyr::group_by(cluster) |>
      dplyr::slice_head(n = n) |>
      dplyr::ungroup()
    write.csv(
      top,
      file.path(directory, paste0("markers_top", n, ".csv")),
      row.names = FALSE
    )
    saved_tables[[paste0("markers_top", n, ".csv")]] <- top
  }
  counts <- table(Idents(alldata))
  summary <- data.frame(
    cluster = names(counts),
    n_cells = as.integer(counts),
    n_prefiltered = as.integer(table(markers$cluster)),
    n_significant = as.integer(table(significant$cluster)),
    n_adjusted_p_zero = as.integer(table(markers$cluster[
      markers$p_val_adj == 0
    ]))
  )
  composition <- as.data.frame(
    table(sample = alldata$sample, cluster = alldata[[sel.clust, drop = TRUE]]),
    responseName = "n_cells"
  )
  write.csv(
    summary,
    file.path(directory, "cluster_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    composition,
    file.path(directory, "sample_cluster_counts.csv"),
    row.names = FALSE
  )
  saved_tables$cluster_summary.csv <- summary
  saved_tables$sample_cluster_counts.csv <- composition

  # Reopen the actual CSVs and compare every column with its source table.
  for (csv in pipeline_csvs) {
    expected <- as.data.frame(saved_tables[[csv]])
    # Unmapped-name columns may be entirely blank; retain them as character.
    reopened <- read.csv(
      file.path(directory, csv),
      stringsAsFactors = FALSE,
      check.names = FALSE,
      colClasses = ifelse(
        vapply(expected, is.character, logical(1)),
        "character",
        NA_character_
      )
    )
    for (column in names(expected)) {
      if (is.factor(expected[[column]]) || is.character(expected[[column]])) {
        expected[[column]] <- as.character(expected[[column]])
        reopened[[column]] <- as.character(reopened[[column]])
      }
    }
    stopifnot(
      identical(names(reopened), names(expected)),
      isTRUE(all.equal(
        as.list(reopened),
        as.list(expected),
        check.attributes = FALSE,
        tolerance = 1e-12
      ))
    )
  }
  stopifnot(
    sum(summary$n_cells) == ncol(alldata),
    sum(composition$n_cells) == ncol(alldata)
  )
  rm(saved_tables, reopened, expected)

  stopifnot(identical(file.info(input)[, c("size", "mtime")], input_info))
  validation <- c(
    "Validation: all seven annotated CSVs reopened and values matched; cluster/sample counts conserve cells; source RDS size/mtime unchanged. No plots generated.",
    paste(
      "Marker rows:",
      nrow(markers),
      "; significant rows:",
      nrow(significant)
    ),
    capture.output(print(counts))
  )
  completed <- TRUE
  message("Annotated subset marker tables completed: ", directory)
  invisible(directory)
}

# Sourcing defines the runner without executing marker discovery.
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) {
    stop(
      "Usage: Rscript R/clusters/subset_markers.R <data/job-XXX/subsets/<Folder>/<subrun>/subclustered.rds> <resolution>"
    )
  }
  run_subset_markers(args[[1L]], as.numeric(args[[2L]]))
}
