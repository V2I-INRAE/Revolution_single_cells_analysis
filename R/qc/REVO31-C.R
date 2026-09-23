# Import and clean filtered BD Rhapsody MEX for REVO31-C into a Seurat object

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

sample_id <- "REVO31-C"
mex_zip <- file.path(
  "raw_data",
  sample_id,
  paste0(sample_id, "_RSEC_MolsPerCell_MEX.zip")
)

stopifnot("MEX zip not found" = file.exists(mex_zip))

# --- Extract MEX files (Read10X reads a directory, not a zip)

mex_dir <- file.path(tempdir(), paste0(sample_id, "_filtered_MEX"))
utils::unzip(mex_zip, exdir = mex_dir)
stopifnot(all(file.exists(file.path(
  mex_dir,
  c("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz")
))))

# --- Read the filtered matrix and build the Seurat object
counts <- Seurat::Read10X(data.dir = mex_dir) # genes x cells, RSEC counts
colnames(counts) <- paste(sample_id, colnames(counts), sep = "_")

# The BD pig reference names mitochondrial genes without the MT- prefix
# (ATP6, ATP8, COX1-3, CYTB, ND1-6); prepend MT- so "^MT-" regex works.
mito_genes <- c(
  "ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB",
  "ND1", "ND2", "ND3", "ND4", "ND5", "ND6"
)
idx <- match(mito_genes, rownames(counts))
stopifnot(!anyNA(idx))
rownames(counts)[idx] <- paste0("MT-", mito_genes)
stopifnot(anyDuplicated(rownames(counts)) == 0)

obj <- Seurat::CreateSeuratObject(counts = counts, project = sample_id)

# --- QC metrics (inspection only, no filtering)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[LS]") # descriptive only, not used for filtering
obj$log10GenesPerUMI <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)

qs <- function(x) round(quantile(x, c(0, 0.25, 0.5, 0.75, 1)), 2)
summary_tbl <- rbind(
  nCount_RNA = qs(obj$nCount_RNA),
  nFeature_RNA = qs(obj$nFeature_RNA),
  percent.mt = qs(obj$percent.mt),
  percent.ribo = qs(obj$percent.ribo),
  log10GenesPerUMI = qs(obj$log10GenesPerUMI)
)
colnames(summary_tbl) <- c("min", "q25", "median", "q75", "max")
cat("\nQC metric summary:\n")
print(summary_tbl)

# --- QC plots before filtering
violin_plots <- VlnPlot(
  obj,
  features = c(
    "nFeature_RNA", "nCount_RNA", "percent.mt",
    "percent.ribo", "log10GenesPerUMI"
  ),
  ncol = 3, pt.size = 0
) + plot_layout(ncol = 3)

scatter_count_feature <- FeatureScatter(
  obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"
)
scatter_feature_mt <- FeatureScatter(
  obj, feature1 = "nFeature_RNA", feature2 = "percent.mt"
)

metrics <- c(
  "nFeature_RNA", "nCount_RNA", "percent.mt",
  "percent.ribo", "log10GenesPerUMI"
)
md_long <- do.call(rbind, lapply(metrics, function(m) {
  data.frame(metric = m, values = obj[[m]][, 1])
}))
qc_histograms <- ggplot(md_long, aes(values)) +
  geom_histogram(
    bins = 60, fill = "grey25", colour = "white",
    linewidth = 0.1
  ) +
  facet_wrap(~metric, scales = "free") +
  theme_classic(base_size = 12) +
  labs(
    x = NULL, y = "Cells",
    title = paste(sample_id, "- QC metric distributions")
  )

# --- Save plots -> results/qc/ (file names prefixed with the sample name)
plot_dir <- file.path("results", "qc", "before_filtering", sample_id)

save_plot <- function(p, name, width, height) {
  png(file.path(plot_dir, paste0(sample_id, "_", name, ".png")),
    width = width, height = height, res = 300
  )
  print(p)
  dev.off()
}

save_plot(violin_plots, "qc_violin", 3600, 2400)
save_plot(scatter_count_feature, "scatter_count_vs_feature", 2400, 2000)
save_plot(scatter_feature_mt, "scatter_feature_vs_mt", 2400, 2000)
save_plot(qc_histograms, "qc_histograms", 3600, 2000)


# --- Persist for downstream QC
# out_rds <- file.path("output", "import", paste0(sample_id, ".rds"))
# dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
# saveRDS(obj, out_rds)
# cat("Saved:", out_rds, "\n")
