# Run or source from the project root. The source RDS is never overwritten.
# Prepares the working object and counts cells; no gene filtering or DE.
suppressPackageStartupMessages(library(SeuratObject))

input <- "data/job-44542861/subsets/Endothelial/job-44831031/subclustered.rds"
cluster_col <- "harmony_snn_res.0.3"
endothelial <- readRDS(input)
metadata <- endothelial[[]]
stopifnot(
  identical(rownames(metadata), colnames(endothelial)),
  !anyDuplicated(rownames(metadata)),
  !anyNA(metadata[, c("sample", "pig", "time_point", cluster_col)])
)

# Preserve sample IDs, including REVO31-N8; only the analysis time changes.
endothelial$time_point[endothelial$time_point == "T8H"] <- "T10H"
cell_counts <- table(
  sample = endothelial$sample,
  cluster = endothelial[[cluster_col, drop = TRUE]]
)
stopifnot(sum(cell_counts) == ncol(endothelial))

cat("Cells relabeled from T8H to T10H:", sum(metadata$time_point == "T8H"), "\n")
cat("Cell counts at", cluster_col, "(including zero combinations):\n")
print(addmargins(cell_counts))
