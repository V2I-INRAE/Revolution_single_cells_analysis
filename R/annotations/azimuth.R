# Pan-human Azimuth cell-level annotation of one saved Seurat object, run locally
# (AzimuthAPI::ANNotate -> reticulate -> panhumanpy). Settings: R/annotations/config.R.
# Usage: Rscript R/annotations/azimuth.R --input=<object.rds> [--output-dir=<dir>]
# Default output: <input folder>/annotations/azimuth_ph_v1/. The input is read-only.
# Requires the uv environment (cd R/annotations && uv sync) and the orthology
# snapshot (Rscript R/annotations/build_orthology_reference.R).
# The model is human: pig features enter only through verified one-to-one orthologues;
# model genes without a mapped feature are zero-filled by panhumanpy.
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})
source("R/annotations/config.R")

# Pig feature -> human model gene. Every feature is kept with its status.
map_features <- function(features, orthology, panel) {
  ensembl <- grepl("^ENSSSCG[0-9]+$", features)
  # R/io.R prefixes pig mitochondrial genes with MT- on load; Ensembl names lack it.
  lookup <- sub("^MT-", "", features)
  named <- unique(orthology[
    nzchar(orthology$pig_gene_name),
    c("pig_gene_id", "pig_gene_name")
  ])
  genes_per_name <- table(named$pig_gene_name)
  pig_id <- ifelse(ensembl, features, NA_character_)
  unique_name <- !ensembl &
    lookup %in% names(genes_per_name)[genes_per_name == 1L]
  pig_id[unique_name] <- named$pig_gene_id[match(
    lookup[unique_name],
    named$pig_gene_name
  )]

  hits <- orthology[nzchar(orthology$human_gene_id), ]
  hits_per_gene <- table(hits$pig_gene_id)
  row <- match(pig_id, hits$pig_gene_id)
  human <- hits$human_gene_name[row]
  status <- rep("used", length(features))
  status[
    !ensembl & !lookup %in% names(genes_per_name)
  ] <- "pig_symbol_not_in_ensembl"
  status[
    !ensembl & lookup %in% names(genes_per_name)[genes_per_name > 1L]
  ] <- "pig_symbol_multiple_genes"
  status[
    ensembl & !features %in% orthology$pig_gene_id
  ] <- "pig_id_not_in_ensembl"
  pending <- status == "used"
  status[pending & is.na(row)] <- "no_human_orthologue"
  pending <- status == "used"
  status[
    pending & as.integer(hits_per_gene[pig_id]) > 1L
  ] <- "multiple_human_orthologues"
  pending <- status == "used"
  status[
    pending & hits$orthology_type[row] != "ortholog_one2one"
  ] <- "not_one2one"
  pending <- status == "used"
  status[pending & !ensembl & human != features] <- "orthologue_symbol_differs"
  pending <- status == "used"
  status[pending & !human %in% panel] <- "not_in_model_panel"
  pending <- status == "used"
  claimed <- human[pending]
  status[
    pending & human %in% claimed[duplicated(claimed)]
  ] <- "human_gene_collision"
  human[status != "used"] <- NA_character_
  data.frame(
    feature = features,
    pig_gene_id = pig_id,
    human_gene_name = human,
    status = status,
    stringsAsFactors = FALSE
  )
}

# Normalized expression of the mapped features from every RNA data layer, all cells.
# Features absent from a sample layer are zero, as in Seurat's JoinLayers.
mapped_expression <- function(object, mapping) {
  used <- mapping[mapping$status == "used", ]
  layers <- Layers(object[["RNA"]], search = "^data($|[.])")
  parts <- lapply(layers, function(layer) {
    m <- LayerData(object, assay = "RNA", layer = layer)
    # Guards against counts, scaled data or a second log transform reaching the model.
    check <- colSums(expm1(m[, seq_len(min(50L, ncol(m))), drop = FALSE]))
    if (any(abs(check - 1e4) > 1)) {
      stop(layer, " is not LogNormalize (scale factor 10000) expression")
    }
    missing <- setdiff(used$feature, rownames(m))
    if (length(missing)) {
      m <- rbind(
        m,
        Matrix(
          0,
          length(missing),
          ncol(m),
          sparse = TRUE,
          dimnames = list(missing, colnames(m))
        )
      )
    }
    list(data = m[used$feature, , drop = FALSE], absent = length(missing))
  })
  absent <- setNames(vapply(parts, `[[`, integer(1), "absent"), layers)
  x <- do.call(cbind, lapply(parts, `[[`, "data"))
  if (anyDuplicated(colnames(x)) || !setequal(colnames(x), colnames(object))) {
    stop("Data layers do not cover every cell exactly once")
  }
  x <- x[, colnames(object)]
  rownames(x) <- used$human_gene_name
  list(data = x, absent = absent)
}

run_azimuth <- function(input, output_dir = NULL) {
  cfg <- annotation_config$azimuth
  if (is.null(output_dir)) {
    output_dir <- file.path(dirname(input), "annotations", cfg$output_folder)
  }
  if (dir.exists(output_dir)) {
    stop("Output exists, remove it first: ", output_dir)
  }
  input_info <- file.info(input)[, c("size", "mtime")]

  reticulate::use_python(cfg$python, required = TRUE)
  panhumanpy <- reticulate::import("panhumanpy")
  panel <- readLines(file.path(
    dirname(panhumanpy$`__file__`),
    "_tools",
    cfg$model_version,
    "inference_feature_panel",
    "inference_feature_panel.txt"
  ))
  panel <- panel[nzchar(panel)]
  orthology <- read.delim(
    annotation_config$orthology$file,
    colClasses = "character",
    quote = "",
    na.strings = character()
  )

  message("Loading (read-only): ", input)
  object <- readRDS(input)
  mapping <- map_features(rownames(object[["RNA"]]), orthology, panel)
  expression <- mapped_expression(object, mapping)
  n_used <- sum(mapping$status == "used")
  message(
    "Model genes with a mapped feature: ",
    n_used,
    "/",
    length(panel),
    sprintf(" (%.1f%%)", 100 * n_used / length(panel))
  )
  print(table(mapping$status))
  message(
    "Mapped features absent per sample layer (zero-filled): ",
    paste(range(expression$absent), collapse = "-")
  )

  # A data-only query has no counts layer for nCount/nFeature; Seurat warns, as in
  # AzimuthAPI::CloudAzimuth, which builds its query the same way.
  query <- suppressWarnings(CreateSeuratObject(CreateAssay5Object(
    data = expression$data
  )))
  warn <- getOption("warn") # ANNotate sets options(warn = -1) globally
  started <- Sys.time()
  # Expression was verified above, so panhumanpy's integer-count heuristic is bypassed;
  # process_obj = FALSE keeps every per-cell label (no rare-label merging, no Idents).
  annotated <- AzimuthAPI::ANNotate(
    query,
    assay = "RNA",
    model_version = cfg$model_version,
    normalization_override = TRUE,
    extract_embeddings = FALSE,
    umap_embeddings = FALSE,
    process_obj = FALSE
  )
  options(warn = warn)
  message("Inference time: ", format(round(Sys.time() - started, 1)))

  returned <- setdiff(colnames(annotated[[]]), colnames(query[[]]))
  renamed <- c(
    full_hierarchical_labels = "hierarchy",
    final_level_labels = "label",
    final_level_confidence = "confidence",
    full_consistent_hierarchy = "consistent",
    azimuth_broad = "broad",
    azimuth_medium = "medium",
    azimuth_fine = "fine"
  )
  if (!all(names(renamed) %in% returned)) {
    stop(
      "Missing Azimuth outputs: ",
      paste(setdiff(names(renamed), returned), collapse = ", ")
    )
  }
  predictions <- annotated[[]][colnames(object), returned, drop = FALSE]
  names(predictions) <- paste0(
    cfg$column_prefix,
    ifelse(returned %in% names(renamed), renamed[returned], returned)
  )
  predictions <- data.frame(
    cell = colnames(object),
    predictions,
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  label <- paste0(cfg$column_prefix, "label")
  confidence <- predictions[[paste0(cfg$column_prefix, "confidence")]]
  if (
    anyNA(predictions[[label]]) ||
      anyNA(confidence) ||
      any(confidence < 0 | confidence > 1)
  ) {
    stop("Missing labels or confidence outside [0, 1]")
  }

  dir.create(output_dir, recursive = TRUE)
  saveRDS(predictions, file.path(output_dir, "predictions.rds"))
  write.csv(
    predictions,
    file.path(output_dir, "predictions.csv"),
    row.names = FALSE
  )
  write.csv(
    mapping,
    file.path(output_dir, "gene_mapping.csv"),
    row.names = FALSE
  )

  # Reopen the saved outputs; predictions must align with the input's cells.
  reopened <- readRDS(file.path(output_dir, "predictions.rds"))
  csv <- read.csv(file.path(output_dir, "predictions.csv"), check.names = FALSE)
  stopifnot(
    identical(reopened, predictions),
    identical(reopened$cell, colnames(object)),
    nrow(csv) == ncol(object),
    identical(csv$cell, colnames(object)),
    identical(file.info(input)[, c("size", "mtime")], input_info)
  )
  message(
    "Saved ",
    nrow(reopened),
    " cells to ",
    output_dir,
    "; input unchanged."
  )
  # "False" in broad/medium/fine marks an inconsistent hierarchy (consistent == FALSE).
  print(sort(
    table(reopened[[paste0(cfg$column_prefix, "medium")]]),
    decreasing = TRUE
  ))
  invisible(predictions)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (
    !length(args) ||
      length(args) > 2L ||
      any(!grepl("^--(input|output-dir)=.+$", args)) ||
      anyDuplicated(sub("=.*$", "", args)) ||
      !any(startsWith(args, "--input="))
  ) {
    stop(
      "Usage: Rscript R/annotations/azimuth.R --input=<object.rds> [--output-dir=<dir>]"
    )
  }
  value <- function(name) {
    flag <- grep(paste0("^--", name, "="), args, value = TRUE)
    if (length(flag)) sub("^[^=]+=", "", flag) else NULL
  }
  run_azimuth(value("input"), value("output-dir"))
}
