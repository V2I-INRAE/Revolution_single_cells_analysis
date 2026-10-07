read_soupx_mex <- function(archive, directory) {
  utils::unzip(archive, exdir = directory)
  features <- read.delim(gzfile(file.path(directory, "features.tsv.gz")),
    header = FALSE, colClasses = "character", quote = "", check.names = FALSE)
  barcodes <- readLines(gzfile(file.path(directory, "barcodes.tsv.gz")))
  connection <- gzfile(file.path(directory, "matrix.mtx.gz"), "rt")
  on.exit(close(connection))
  counts <- as(Matrix::readMM(connection), "CsparseMatrix")
  stopifnot(
    nrow(counts) == nrow(features), ncol(counts) == length(barcodes),
    ncol(features) == 3L, all(features[[3]] == "Gene Expression"),
    !anyDuplicated(features[[1]]), !anyDuplicated(barcodes),
    all(nzchar(features[[1]])), all(nzchar(features[[2]])),
    all(nzchar(barcodes)), all(is.finite(counts@x)),
    all(counts@x >= 0), all(counts@x == floor(counts@x))
  )
  dimnames(counts) <- list(features[[1]], barcodes)
  list(counts = counts, features = features)
}

align_soupx_counts <- function(raw, filtered) {
  rows <- match(rownames(filtered$counts), rownames(raw$counts))
  columns <- match(colnames(filtered$counts), colnames(raw$counts))
  stopifnot(!anyNA(rows), !anyNA(columns),
    identical(filtered$features[[2]], raw$features[[2]][rows]))
  # The raw called-cell columns are the zero-padded filtered matrix only
  # after both shared counts and every omitted row have been verified.
  toc <- raw$counts[, columns, drop = FALSE]
  difference <- Matrix::drop0(toc[rows, , drop = FALSE] - filtered$counts)
  omitted <- setdiff(seq_len(nrow(toc)), rows)
  stopifnot(length(difference@x) == 0L,
    sum(toc[omitted, , drop = FALSE]) == 0)
  # Apply symbol repair once, on the full feature universe, after ID matching.
  feature_names <- make.unique(raw$features[[2]])
  rownames(raw$counts) <- rownames(toc) <- feature_names
  list(tod = raw$counts, toc = toc, feature_map = data.frame(
    feature_id = raw$features[[1]], symbol = raw$features[[2]],
    soupx_name = feature_names, present_in_filtered = seq_len(nrow(toc)) %in% rows),
    zero_padded_features = length(omitted))
}

read_putative_cell_count <- function(path) {
  lines <- readLines(path)
  start <- match("#Cells#", lines)
  stopifnot(!is.na(start))
  section <- lines[(start + 1L):length(lines)]
  end <- which(!nzchar(section) | startsWith(section, "#"))[1]
  if (!is.na(end)) section <- section[seq_len(end - 1L)]
  cells <- read.csv(text = paste(section, collapse = "\n"))
  count <- cells$Putative_Cell_Count[cells$Bioproduct_Type == "mRNA"]
  stopifnot(length(count) == 1L, count > 0)
  as.integer(count)
}

soupx_background_check <- function(tod, toc, expected_cells) {
  totals <- Matrix::colSums(tod)
  in_range <- totals > 0 & totals < 100
  knee <- sort(totals, decreasing = TRUE)[expected_cells]
  n_background <- sum(in_range & !names(totals) %in% colnames(toc))
  list(soup_range = c(0, 100), expected_cells = expected_cells,
    molecules_at_called_cell_rank = unname(knee),
    knee_margin = unname(knee / 100),
    barcodes_in_range = sum(in_range), noncalled_barcodes_in_range = n_background,
    called_barcodes_in_range = sum(in_range & names(totals) %in% colnames(toc)),
    passed = isTRUE(knee > 1000 && n_background >= 2000 &&
      !any(in_range & names(totals) %in% colnames(toc))))
}

validate_soupx_counts <- function(corrected, original) {
  stopifnot(inherits(corrected, "sparseMatrix"),
    identical(dimnames(corrected), dimnames(original)),
    all(is.finite(corrected@x)), all(corrected@x >= 0),
    all(corrected@x == floor(corrected@x)),
    all((corrected - original)@x <= 0))
}
