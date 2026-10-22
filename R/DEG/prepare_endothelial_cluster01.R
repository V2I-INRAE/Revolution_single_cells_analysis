# Shared count preparation for the cluster 1 diagnostics and paired DE pilot.
# Expects input and cluster_col; creates counts, sample_info, keep and TMM-normalized y.
# No model fitting, DE tests, plots or file writes.
endothelial <- readRDS(input)
endothelial$time_point[endothelial$time_point == "T8H"] <- "T10H"

# Define the approved gene universe before removing the late samples.
count_layers <- Layers(endothelial[["RNA"]], search = "^counts\\.")
common_genes <- Reduce(intersect, lapply(count_layers, function(layer) {
  Features(endothelial[["RNA"]], layer = layer)
}))
cells <- colnames(endothelial)[
  endothelial[[cluster_col, drop = TRUE]] == "1" &
    endothelial$time_point %in% c("T0H", "T4H")
]
pilot <- subset(endothelial, cells = cells, features = common_genes)
counts <- AggregateExpression(pilot, assays = "RNA", group.by = "sample",
  return.seurat = FALSE, verbose = FALSE)$RNA

sample_info <- unique(pilot[[]][, c("sample", "pig", "pressure", "time_point")])
sample_info <- sample_info[match(colnames(counts), sample_info$sample), ]
rownames(sample_info) <- sample_info$sample
sample_info$pig <- factor(sample_info$pig)
sample_info$condition <- factor(
  ifelse(sample_info$time_point == "T0H", "Control",
    ifelse(sample_info$pressure == "positive", "P4", "N4")),
  levels = c("Control", "P4", "N4")
)
sample_info$n_cells <- as.integer(table(pilot$sample)[sample_info$sample])
design <- model.matrix(~ pig + condition, data = sample_info)

y <- DGEList(counts = counts, samples = sample_info)
keep <- filterByExpr(y, group = sample_info$condition)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- normLibSizes(y, method = "TMM")
