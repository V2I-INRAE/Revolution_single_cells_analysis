# Small artificial software fixtures only; these are not biological analysis data.
source("scripts/soupx/io.R")
suppressPackageStartupMessages(library(Matrix))
expect_error <- function(expression) {
  stopifnot(inherits(tryCatch(force(expression), error = identity), "error"))
}

# Exercise actual ZIP/gzip/MatrixMarket IO, including numeric-looking barcodes.
local({
  directory <- tempfile("soupx-io-test-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  contents <- list(
    "matrix.mtx.gz" = c("%%MatrixMarket matrix coordinate integer general", "%",
      "2 2 3", "1 1 7", "2 1 2", "2 2 11"),
    "features.tsv.gz" = c("id1\tsymbol1\tGene Expression",
      "id2\tsymbol2\tGene Expression"),
    "barcodes.tsv.gz" = c("01", "02"))
  for (name in names(contents)) {
    connection <- gzfile(file.path(directory, name), "wt")
    writeLines(contents[[name]], connection)
    close(connection)
  }
  archive <- file.path(directory, "fixture.zip")
  utils::zip(archive, file.path(directory, names(contents)), flags = "-j -q")
  imported <- read_soupx_mex(archive, file.path(directory, "extracted"))
  stopifnot(identical(colnames(imported$counts), c("01", "02")),
    identical(rownames(imported$counts), c("id1", "id2")),
    identical(as.numeric(imported$counts), c(7, 2, 0, 11)))
})

# Asymmetric counts and reversed filtered orders expose positional alignment bugs.
counts <- sparseMatrix(i = c(1, 2, 1, 2, 3), j = c(1, 1, 2, 2, 3),
  x = c(7, 2, 3, 11, 1), dims = c(3, 3),
  dimnames = list(c("id1", "id2", "id3"), c("01", "02", "03")))
features <- data.frame(id = c("id1", "id2", "id3"),
  symbol = c("same", "same", "ambient"), type = "Gene Expression")
raw <- list(counts = counts, features = features)
filtered <- list(counts = counts[c(2, 1), c(2, 1)], features = features[c(2, 1), ])
aligned <- align_soupx_counts(raw, filtered)
stopifnot(identical(colnames(aligned$toc), c("02", "01")),
  identical(rownames(aligned$toc), c("same", "same.1", "ambient")),
  identical(as.numeric(aligned$toc[, 1]), c(3, 11, 0)),
  aligned$zero_padded_features == 1L)
bad <- raw
bad$counts[3, 1] <- 1
expect_error(align_soupx_counts(bad, filtered))
bad <- filtered
bad$counts[1, 1] <- 10
expect_error(align_soupx_counts(raw, bad))

corrected <- aligned$toc
corrected[1, 1] <- 1
validate_soupx_counts(corrected, aligned$toc)
bad <- corrected
bad[2, 1] <- 12
expect_error(validate_soupx_counts(bad, aligned$toc))
bad <- corrected
bad[1, 1] <- 0.5
expect_error(validate_soupx_counts(bad, aligned$toc))
expect_error(validate_soupx_counts(corrected[, 2:1], aligned$toc))

# Exact boundaries: 0 and 100 excluded, 99 included, >1,000 knee and >=2,000 wells.
totals <- c(1001, rep(99, 2000), 100, 0)
tod <- sparseMatrix(i = rep(1L, length(totals)), j = seq_along(totals),
  x = totals, dims = c(1, length(totals)),
  dimnames = list("gene", paste0("barcode", seq_along(totals))))
toc <- tod[, 1, drop = FALSE]
check <- soupx_background_check(tod, toc, 1L)
stopifnot(check$passed, check$noncalled_barcodes_in_range == 2000,
  check$molecules_at_called_cell_rank == 1001)
bad <- tod
bad[1, 1] <- 1000
stopifnot(!soupx_background_check(bad, toc, 1L)$passed)
bad <- tod
bad[1, 2] <- 100
stopifnot(!soupx_background_check(bad, toc, 1L)$passed)
stopifnot(!soupx_background_check(tod, tod[, 2, drop = FALSE], 1L)$passed)
stopifnot(read_putative_cell_count(
  "data/raw_data/REVO30-P4/REVO30-P4_Metrics_Summary.csv") == 13422L)
cat("SoupX MEX reader, alignment, count-integrity, background-boundary and metrics checks passed.\n")
