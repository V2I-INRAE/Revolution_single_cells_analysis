# Raw-count reading, run-directory layout and .INFO run-record lifecycle.
# The pipeline owns exactly one data/job-<id>/.INFO (authoritative) plus a
# mirrored copy in results/job-<id>/.INFO; both are rewritten from a single
# in-memory record so they never disagree.

suppressPackageStartupMessages({
  library(Seurat)
})

# --- Raw input paths ---------------------------------------------------------

qc_input_path <- function(sample_id, input_source) {
  folder <- file.path("data", "raw_data", input_source, sample_id)
  switch(input_source,
    rhapsody = file.path(folder, paste0(sample_id, "_RSEC_MolsPerCell_MEX.zip")),
    soupx = file.path(folder, "job-44409807", "corrected_counts.rds"),
    cellbender = file.path(folder, paste0(sample_id, "_cellbender_FPR_0.01_filtered.h5")),
    stop("Unknown QC input source: ", input_source)
  )
}

input_source_notes <- function(input_source) {
  switch(input_source,
    soupx = "Ambient correction producer: job-44409807; corrected counts reused, not rerun.",
    cellbender = "FPR 0.01 filtered counts; source runs and unresolved warnings: data/raw_data/cellbender/.INFO.",
    NULL
  )
}

read_qc_counts <- function(sample_id, input_source) {
  path <- qc_input_path(sample_id, input_source)
  if (!file.exists(path)) stop("QC input not found: ", path)
  switch(input_source,
    rhapsody = {
      mex_dir <- tempfile(paste0(sample_id, "_filtered_MEX_"))
      on.exit(unlink(mex_dir, recursive = TRUE))
      utils::unzip(path, exdir = mex_dir)
      Seurat::Read10X(data.dir = mex_dir)
    },
    soupx = readRDS(path),
    cellbender = {
      # Read only expression counts: Read10X_h5 also tries to read latent groups.
      h5 <- hdf5r::H5File$new(path, mode = "r")
      on.exit(h5$close_all())
      matrix <- h5[["matrix"]]
      Matrix::sparseMatrix(
        i = matrix[["indices"]][] + 1L,
        p = matrix[["indptr"]][],
        x = as.numeric(matrix[["data"]][]),
        dims = matrix[["shape"]][],
        dimnames = list(
          make.unique(matrix[["features/name"]][]),
          matrix[["barcodes"]][]
        )
      )
    }
  )
}

# --- Seurat construction -----------------------------------------------------

build_seurat_obj_from_counts <- function(counts, sample_id, qc_config) {
  # sample id reaches the object three ways: the barcodes (unique cell
  # names across samples, required for the final merge), orig.ident
  # (project) and the explicit `sample` metadata column set below
  colnames(counts) <- paste(sample_id, colnames(counts), sep = "_")

  # --- Pig mitochondrial gene names do not have the MT-prefix
  mito_genes <- qc_config$metrics$mitochondrial_genes
  idx <- match(mito_genes, rownames(counts))
  stopifnot(!anyNA(idx))
  rownames(counts)[idx] <- paste0("MT-", mito_genes)
  stopifnot(anyDuplicated(rownames(counts)) == 0)

  seurat_obj <- Seurat::CreateSeuratObject(counts = counts, project = sample_id)

  # --- Sample metadata parsed from the sample id (REVOxx-CODE)
  seurat_obj$sample <- sample_id
  sample_parts <- strsplit(sample_id, "-", fixed = TRUE)[[1]]
  condition <- sample_parts[[2]]
  seurat_obj$pig <- sample_parts[[1]]
  seurat_obj$pressure <- switch(substr(condition, 1, 1),
    C = "None",
    P = "positive",
    N = "negative",
    stop("Unknown pressure code in sample ID: ", sample_id)
  )
  seurat_obj$time_point <- if (condition == "C") {
    "T0H"
  } else {
    paste0("T", substring(condition, 2), "H")
  }

  return(seurat_obj)
}

concatenate_qc_objects <- function(seurat_objs, qc_config) {
  concatenated <- merge(seurat_objs[[1]], y = seurat_objs[-1], merge.data = FALSE)
  concatenated@misc$qc <- list(
    sample_order = names(seurat_objs),
    params = qc_config,
    samples = lapply(seurat_objs, function(x) x@misc$qc),
    session = sessionInfo(),
    created = Sys.time()
  )
  return(concatenated)
}

# --- Run layout ---------------------------------------------------------------

plot_subdirs <- c("qc_plots", "norm_plots", "integration_plots", "clusters_plots")

# Creates data/job-<id>, results/job-<id>/{plot subdirs}; refuses existing
# run directories so a requeue can never overwrite a previous run.
create_run_layout <- function(run_id) {
  data_dir <- file.path("data", run_id)
  results_dir <- file.path("results", run_id)
  for (directory in c(data_dir, results_dir)) {
    if (dir.exists(directory)) {
      stop("Run directory already exists: ", directory)
    }
  }
  plot_dirs <- setNames(
    file.path(results_dir, plot_subdirs),
    plot_subdirs
  )
  for (directory in c(data_dir, unname(plot_dirs))) {
    if (!dir.create(directory, recursive = TRUE, showWarnings = FALSE)) {
      stop("Cannot create run directory: ", directory)
    }
  }
  list(
    run_id = run_id,
    data_dir = data_dir,
    results_dir = results_dir,
    plot_dirs = plot_dirs,
    info_paths = c(
      data = file.path(data_dir, ".INFO"),
      results = file.path(results_dir, ".INFO")
    )
  )
}

# --- Run record (.INFO) -------------------------------------------------------

new_run_record <- function(config, input_paths) {
  git_commit <- tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[[1]],
    error = function(e) NA_character_
  )
  list(
    run_id = NULL, # set by the caller once the layout exists
    status = "running",
    failed_stage = NULL,
    error = NULL,
    started = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    finished = NULL,
    code = paste(
      if (is.na(git_commit)) "git revision unavailable" else paste0("git ", git_commit),
      paste0("R ", getRversion()),
      paste0("Seurat ", utils::packageVersion("Seurat")),
      paste0("harmony ", utils::packageVersion("harmony")),
      sep = "; "
    ),
    config = config,
    input_paths = input_paths,
    input_notes = input_source_notes(config$run$input_source),
    stages = list(),
    final_object = NULL
  )
}

record_stage <- function(record, stage, lines) {
  record$stages[[stage]] <- c(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    lines
  )
  record
}

finish_run_record <- function(record) {
  record$status <- "completed"
  record$finished <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  record
}

fail_run_record <- function(record, stage, error) {
  record$status <- paste0("failed at stage: ", stage)
  record$failed_stage <- stage
  record$finished <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  record$error <- conditionMessage(error)
  record
}

format_value <- function(value) {
  if (is.data.frame(value)) {
    paste(apply(value, 1, paste, collapse = ":"), collapse = ", ")
  } else if (is.list(value)) {
    paste(names(value), vapply(value, format_value, character(1)),
      sep = "=", collapse = "; "
    )
  } else if (!is.null(names(value)) && any(nzchar(names(value)))) {
    paste(names(value), vapply(value, format_value, character(1)),
      sep = "=", collapse = ", "
    )
  } else if (length(value) == 1L) {
    as.character(value)
  } else {
    paste(value, collapse = ", ")
  }
}

format_config <- function(config) {
  sections <- c(
    # dims is a bare vector, unlike the other (named-list) sections.
    list(shared = list(dims = config$dims)),
    config[c("run", "qc", "norm", "cell_cycle", "integration", "clustering")]
  )
  lines <- unlist(lapply(names(sections), function(name) {
    c(paste0("### ", name),
      paste0("- ", names(sections[[name]]), ": ",
        vapply(sections[[name]], format_value, character(1))))
  }))
  # One blank line between sections for scannability.
  lines <- unlist(lapply(seq_along(lines), function(i) {
    if (i > 1L && startsWith(lines[[i]], "###")) c("", lines[[i]]) else lines[[i]]
  }))
  lines
}

render_info <- function(record) {
  with(record, {
    c(
      paste0("# Pipeline run ", run_id),
      "",
      "## Run",
      paste0("- Status: ", status),
      paste0("- Started: ", started),
      if (!is.null(finished)) paste0("- Finished: ", finished) else NULL,
      paste0("- Input source: ", config$run$input_source),
      paste0("- Code: ", code),
      if (!is.null(error)) paste0("- Error: ", error) else NULL,
      "",
      "## Inputs",
      paste0("- ", names(input_paths), ": ", input_paths),
      if (!is.null(input_notes)) paste0("- ", input_notes) else NULL,
      "",
      "## Settings",
      format_config(config),
      "",
      "## Stages",
      if (length(stages) == 0L) {
        "- (none completed yet)"
      } else {
        unlist(lapply(names(stages), function(name) {
          c(paste0("### ", name), paste0("- ", stages[[name]]))
        }))
      },
      "",
      "## Outputs",
      paste0("- Final object: ",
        if (is.null(final_object)) "pending" else final_object),
      paste0("- Figures: results/", run_id, "/{",
        paste(plot_subdirs, collapse = ","), "}")
    )
  })
}

# Rewrite both .INFO copies from the single in-memory record. The data/ copy
# is authoritative; a crash between the two writes can at worst leave the
# results copy one update behind.
write_info <- function(record, info_paths) {
  lines <- render_info(record)
  for (path in info_paths) writeLines(lines, path)
  invisible(lines)
}

# --- Final object --------------------------------------------------------------

# Atomic save: a crash during saveRDS leaves no half-written final object.
save_final_object <- function(seurat_obj, path) {
  tmp <- paste0(path, ".tmp")
  saveRDS(seurat_obj, tmp)
  if (!file.rename(tmp, path)) {
    stop("Could not finalize final object save: ", path)
  }
  invisible(path)
}
