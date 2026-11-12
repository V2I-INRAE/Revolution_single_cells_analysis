# Paired pseudobulk DEG and GSEA for one subset resolution. Settings: R/DEG/config.R.
# Usage: Rscript R/DEG/main.R --input=<subclustered.rds> --subset=<label>
#   --cluster-col=harmony_snn_res.X [--clusters=all|1,2] [--output-root=<res0.X/DEG>]
#   [--gsea-exclusions=<approved CSV>]
# Every run replaces <output_root> and the DEG_<contrast> publication figures.
suppressPackageStartupMessages({
  library(Seurat)
  library(edgeR)
  library(ggplot2)
})
source("R/DEG/config.R")
source("R/utils/save_plots.R")
source("R/utils/pig_gene_mapping.R")
source("R/DEG/pseudobulk_de.R")
source("R/DEG/pseudobulk_gsea.R")
source("R/plots/publications/plot_DEG_for_pub.R")
source("R/plots/publications/publication_paths.R")

run_pseudobulk <- function(
  input,
  subset_label,
  cluster_col,
  clusters = "all",
  output_root = NULL,
  gsea_exclusions = list()
) {
  resolution <- sub("^harmony_snn_res[.]", "", cluster_col)
  if (is.null(output_root)) {
    output_root <- file.path(subset_resolution_dir(input, resolution), "DEG")
  }
  # output_root is deleted recursively below; never let it point at the resolution folder.
  if (basename(output_root) != "DEG") {
    stop("output_root must be a DEG folder: ", output_root)
  }
  publication_dir <- file.path(dirname(output_root), "publication_figures")
  contrasts <- pseudobulk_comparisons()$contrast

  object <- readRDS(input)
  metadata <- object[[]]
  observed <- as.character(sort(unique(as.integer(as.character(metadata[[
    cluster_col
  ]])))))
  clusters <- if (identical(clusters, "all")) {
    observed
  } else {
    as.character(clusters)
  }
  # A mistyped cluster ID would otherwise give a silently all-skipped cluster.
  if (!all(clusters %in% observed)) {
    stop(
      "Unknown clusters: ",
      paste(setdiff(clusters, observed), collapse = ", ")
    )
  }

  sample_panel <- unique(metadata[, c(
    "sample",
    "pig",
    "pressure",
    "time_point"
  )])
  sample_panel$sample <- as.character(sample_panel$sample)
  sample_panel$time_point_original <- as.character(sample_panel$time_point)
  # REVO31-N8 (T8H) is analysed with the T10H specimens.
  sample_panel$time_point <- sub(
    "^T8H$",
    "T10H",
    sample_panel$time_point_original
  )
  conditions <- c(
    T0H_None = "Control",
    T4H_positive = "P4",
    T4H_negative = "N4",
    T10H_positive = "P10",
    T10H_negative = "N10"
  )
  sample_panel$condition <- factor(
    unname(conditions[paste(
      sample_panel$time_point,
      sample_panel$pressure,
      sep = "_"
    )]),
    levels = unname(conditions)
  )
  # An unknown label would silently drop the specimen from every contrast.
  if (anyNA(sample_panel$condition)) {
    stop("Unexpected time/pressure combination in the sample panel.")
  }
  # The relabelling must not give one pig two specimens in the same condition.
  if (anyDuplicated(sample_panel[, c("pig", "condition")])) {
    stop("Two specimens map to one pig/condition.")
  }
  sample_panel$pig <- factor(sample_panel$pig)
  rownames(sample_panel) <- sample_panel$sample

  # Gene universe: features present in every original count layer.
  count_layers <- SeuratObject::Layers(object[["RNA"]], search = "^counts[.]")
  common_genes <- Reduce(
    intersect,
    lapply(count_layers, function(layer) {
      SeuratObject::Features(object[["RNA"]], layer = layer)
    })
  )
  gene_mapping <- prepare_pig_gene_mapping(
    rownames(object[["RNA"]]),
    reference_dir = deg_config$reference_dir
  )
  references <- lapply(
    setNames(nm = deg_config$collections),
    function(collection) {
      prepare_pseudobulk_gsea(common_genes, collection)
    }
  )

  # Idempotent: replace every output of a previous run.
  unlink(output_root, recursive = TRUE)
  unlink(file.path(
    publication_dir,
    as.vector(outer(
      paste0("DEG_", contrasts),
      c(".pdf", ".png", ".INFO"),
      paste0
    ))
  ))
  dir.create(output_root, recursive = TRUE)
  write.csv(
    sample_panel,
    file.path(output_root, "sample_panel.csv"),
    row.names = FALSE
  )
  write.csv(
    gene_mapping,
    file.path(output_root, "gene_display_mapping.csv"),
    row.names = FALSE
  )
  write.csv(
    references$hallmark$mapping,
    file.path(output_root, "candidate_gene_mapping.csv"),
    row.names = FALSE
  )

  summaries <- eligibility <- list()
  for (cluster in clusters) {
    directory <- file.path(
      output_root,
      sprintf("cluster_%02d", as.integer(cluster))
    )
    dir.create(directory)
    label <- paste(subset_label, "cluster", cluster)
    cat("\n", label, "started", format(Sys.time()), "\n")
    population <- prepare_pseudobulk(
      object,
      common_genes,
      cluster_col,
      cluster,
      sample_panel
    )
    write.csv(
      population$sample_audit,
      file.path(directory, "population_sample_audit.csv"),
      row.names = FALSE
    )
    rows <- statuses <- list()
    for (contrast in contrasts) {
      comparison_dir <- file.path(directory, contrast)
      dir.create(comparison_dir)
      pb <- prepare_paired_comparison(population, contrast)
      write.csv(
        pb$sample_audit,
        file.path(comparison_dir, "sample_audit.csv"),
        row.names = FALSE
      )
      write.csv(
        pb$sample_info,
        file.path(comparison_dir, "sample_metadata.csv"),
        row.names = FALSE
      )
      s <- data.frame(
        cluster,
        status = pb$status,
        reason = pb$reason,
        cells = sum(pb$sample_info$n_cells),
        samples = nrow(pb$sample_info),
        pigs = nlevels(pb$sample_info$pig),
        min_retained_cells = if (nrow(pb$sample_info)) {
          min(pb$sample_info$n_cells)
        } else {
          NA_integer_
        },
        below_min_cells_samples = sum(
          pb$sample_audit$in_comparison & !pb$sample_audit$passes_cell_filter
        ),
        unmatched_samples = sum(
          pb$sample_audit$exclusion_reason == "unmatched_pig"
        ),
        candidate_genes = length(common_genes),
        design_rank = pb$design_rank,
        residual_df = pb$residual_df,
        filter_min_group_size = pb$filter_min_group_size,
        pb$contrast_status,
        tested = NA_integer_,
        up = NA_integer_,
        down = NA_integer_,
        up_fdr_logfc = NA_integer_,
        down_fdr_logfc = NA_integer_,
        row.names = NULL
      )
      for (collection in deg_config$collections) {
        s[[paste0(collection, "_status")]] <- "not run: skipped contrast"
      }
      if (pb$status == "eligible") {
        pb <- fit_pseudobulk(pb, comparison_dir, gene_mapping)
        tab <- pb$results
        excluded <- unique(c(
          gsea_exclusions[[cluster]],
          gsea_exclusions[[paste(cluster, contrast, sep = "/")]]
        ))
        for (collection in deg_config$collections) {
          gsea <- run_pseudobulk_gsea(
            pb,
            references[[collection]],
            file.path(comparison_dir, "gsea", collection),
            label,
            excluded
          )
          s[[paste0(collection, "_status")]] <- gsea$status
          if (!is.null(gsea$summary)) {
            coverage <- gsea$summary
            coverage$contrast <- NULL
            names(coverage) <- paste0(collection, "_", names(coverage))
            s <- cbind(s, coverage)
          }
        }
        s$status <- "completed"
        s$contrast_status <- "tested"
        s$tested <- pb$summary$tested
        s$up <- pb$summary$up
        s$down <- pb$summary$down
        s$up_fdr_logfc <- sum(
          tab$FDR <= deg_config$deg_fdr & tab$logFC >= deg_config$deg_logfc
        )
        s$down_fdr_logfc <- sum(
          tab$FDR <= deg_config$deg_fdr & tab$logFC <= -deg_config$deg_logfc
        )
      }
      write.csv(
        pb$contrast_status,
        file.path(comparison_dir, "contrast_status.csv"),
        row.names = FALSE
      )
      rows[[contrast]] <- s
      statuses[[contrast]] <- pb$contrast_status
      cat(
        label,
        contrast,
        s$status,
        "; pairs",
        s$matched_pigs,
        "; genes",
        s$tested,
        "\n"
      )
    }
    cluster_summary <- dplyr::bind_rows(rows)
    write.csv(
      dplyr::bind_rows(statuses),
      file.path(directory, "contrast_status.csv"),
      row.names = FALSE
    )
    write.csv(
      cluster_summary,
      file.path(directory, "DEG_summary.csv"),
      row.names = FALSE
    )
    summaries[[cluster]] <- cluster_summary
    eligibility[[cluster]] <- data.frame(
      cluster,
      requested = nrow(cluster_summary),
      tested = sum(cluster_summary$contrast_status == "tested"),
      skipped = sum(cluster_summary$contrast_status == "skipped"),
      observed_cells = sum(population$sample_audit$n_cells),
      samples_passing_min_cells = sum(
        population$sample_audit$passes_cell_filter
      )
    )
  }
  summary <- dplyr::bind_rows(summaries)
  write.csv(
    summary,
    file.path(output_root, "cluster_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    dplyr::bind_rows(eligibility),
    file.path(output_root, "cluster_eligibility.csv"),
    row.names = FALSE
  )
  for (contrast in contrasts) {
    plot_DEG_for_pub(
      output_root,
      publication_dir,
      contrast,
      subset_label,
      resolution
    )
  }
  invisible(summary)
}

if (sys.nframe() == 0L) {
  usage <- paste(
    "Usage: Rscript R/DEG/main.R --input=<subclustered.rds> --subset=<label>",
    "--cluster-col=harmony_snn_res.X [--clusters=all|1,2] [--output-root=<res0.X/DEG>]",
    "[--gsea-exclusions=<approved CSV>]"
  )
  args <- commandArgs(trailingOnly = TRUE)
  keys <- sub("^--([^=]+)=.*$", "\\1", args)
  # A mistyped flag would otherwise be ignored, e.g. GSEA silently run without approved exclusions.
  allowed <- c(
    "input",
    "subset",
    "cluster-col",
    "clusters",
    "output-root",
    "gsea-exclusions"
  )
  if (
    !length(args) ||
      !all(grepl("^--[a-z-]+=.+$", args)) ||
      !all(keys %in% allowed)
  ) {
    stop(usage)
  }
  opt <- as.list(setNames(sub("^--[^=]+=", "", args), keys))
  exclusions <- list()
  if (!is.null(opt$`gsea-exclusions`)) {
    x <- read.csv(opt$`gsea-exclusions`, colClasses = "character")
    if (!"contrast" %in% names(x)) {
      x$contrast <- "*"
    }
    # The approved file lists every subset and resolution.
    if ("subset" %in% names(x)) {
      x <- x[
        x$subset == opt$subset &
          x$resolution == sub("^harmony_snn_res[.]", "", opt$`cluster-col`),
      ]
    }
    # A mistyped contrast would silently drop an approved exclusion.
    if (!all(x$contrast %in% c("*", pseudobulk_comparisons()$contrast))) {
      stop("Unknown contrast in GSEA exclusions.")
    }
    exclusions <- split(
      x$gene,
      ifelse(
        x$contrast == "*",
        x$cluster,
        paste(x$cluster, x$contrast, sep = "/")
      )
    )
  }
  clusters <- if (is.null(opt$clusters) || opt$clusters == "all") {
    "all"
  } else {
    strsplit(opt$clusters, ",")[[1]]
  }
  run_pseudobulk(
    opt$input,
    opt$subset,
    opt$`cluster-col`,
    clusters,
    opt$`output-root`,
    exclusions
  )
}
