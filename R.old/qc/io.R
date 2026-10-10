# io operations

suppressPackageStartupMessages({
  library(Seurat)
})

qc_input_path <- function(sample_id, input_source) {
  folder <- file.path("data", "raw_data", input_source, sample_id)
  switch(input_source,
    rhapsody = file.path(folder, paste0(sample_id, "_RSEC_MolsPerCell_MEX.zip")),
    soupx = file.path(folder, "job-44409807", "corrected_counts.rds"),
    cellbender = file.path(folder, paste0(sample_id, "_cellbender_FPR_0.01_filtered.h5")),
    stop("Unknown QC input source: ", input_source)
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

build_seurat_obj_from_counts <- function(counts, sample_id) {
  # sample id reaches the object three ways: the barcodes (unique cell
  # names across samples, required for the final merge), orig.ident
  # (project) and the explicit `sample` metadata column set below
  colnames(counts) <- paste(sample_id, colnames(counts), sep = "_")

  # --- Pig mitochondrial gene names do not have the MT-prefix
  mito_genes <- qc_params$metrics$mitochondrial_genes
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

concatenate_qc_objects <- function(seurat_objs) {
  concatenated <- merge(seurat_objs[[1]], y = seurat_objs[-1], merge.data = FALSE)
  concatenated@misc$qc <- list(
    sample_order = names(seurat_objs),
    params = qc_params,
    samples = lapply(seurat_objs, function(x) x@misc$qc),
    session = sessionInfo(),
    created = Sys.time()
  )
  return(concatenated)
}

write_qc_object <- function(seurat_obj, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(seurat_obj, path)
}
