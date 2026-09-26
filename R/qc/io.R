# io operations

suppressPackageStartupMessages({
  library(Seurat)
})

read_filtered_matrix <- function(folder_raw, sample_id) {
  # --- get path to the zipped filtered matrix
  mex_zip <- file.path(
    folder_raw,
    sample_id,
    paste0(sample_id, "_RSEC_MolsPerCell_MEX.zip")
  )
  stopifnot("MEX zip not found" = file.exists(mex_zip))

  # --- unzip the matrix to barcodes, features and matrix count
  mex_dir <- file.path(tempdir(), paste0(sample_id, "_filtered_MEX"))
  utils::unzip(mex_zip, exdir = mex_dir)
  stopifnot(all(file.exists(file.path(
    mex_dir,
    c("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz")
  ))))

  return(mex_dir)
}

build_seurat_obj <- function(unzipped_mex_dir, sample_id) {
  counts <- Seurat::Read10X(data.dir = unzipped_mex_dir)

  # sample id reaches the object three ways: the barcodes (unique cell
  # names across samples, required for the final merge), orig.ident
  # (project) and the explicit `sample` metadata column set below
  colnames(counts) <- paste(sample_id, colnames(counts), sep = "_")

  # --- Pig mitochondrial gene names do not have the MT-prefix
  mito_genes <- c(
    "ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB",
    "ND1", "ND2", "ND3", "ND4", "ND5", "ND6"
  )
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

write_concatenated_obj <- function(seurat_objs) {
  concat_dir <- file.path("data", "clean_concatenated_data")
  dir.create(concat_dir, showWarnings = FALSE, recursive = TRUE)
  concatenated <- merge(seurat_objs[[1]], y = seurat_objs[-1])
  saveRDS(concatenated, file.path(concat_dir, "clean_concatenated.rds"))
}
