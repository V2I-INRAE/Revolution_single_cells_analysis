# Coarse-subset UMAP: all saved subsets on one square figure, PNG and PDF.
# Usage: Rscript R/plots/plots_subsets.R <annotated_clustered.rds>
#
# Read-only: the object is loaded for plotting only; nothing is written back.
# Identity handling: Cluster_Highlight_Plot matches `cluster_name` against the
# active identity, so identities are set in memory to the saved subset column.
# The file on disk is never modified.
#
# PDF size note: large cell layers as vector circles can overwhelm vector
# editors, so point layers are rasterized at 300 dpi via ggrastr
# (text, axes and legends stay vector). PNG is rendered by save_figure.
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(scCustomize)
})
source("R/utils/save_plots.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript R/plots/plots_subsets.R <data/job-XXX/annotated_clustered.rds>")
}
input <- args[[1L]]
reduction <- "umap_harmony"
seed <- 1234L # production run seed (results/job-44493014/.INFO)

out_dir <- file.path("results", basename(dirname(input)),
  "clusters_plots", "subset")

if (!file.exists(input)) stop("Input not found: ", input)
existing <- list.files(out_dir, pattern = "^coarse_subsets_all\\.(png|pdf)$")
if (length(existing) > 0L) {
  stop("Combined subset plots already exist in ", out_dir, ": ",
    paste(existing, collapse = ", "),
    "; preserve them before running again.")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
in_mtime <- file.mtime(input)
in_size <- file.size(input)

message("Loading (read-only): ", input)
obj <- readRDS(input)

# Reuse approved annotations; do not redefine subsets from cluster numbers.
SeuratObject::Idents(obj) <- "subset"
subsets <- levels(SeuratObject::Idents(obj))
subset_colors <- stats::setNames(viz_palette[seq_along(subsets)], subsets)
message("Plotting all subsets: ", paste(subsets, collapse = ", "))

set.seed(seed)
p_all <- scCustomize::Cluster_Highlight_Plot(
  seurat_object = obj,
  cluster_name = subsets,
  highlight_color = unname(subset_colors[subsets]),
  background_color = "lightgray",
  reduction = reduction,
  pt.size = 0.1,
  stroke.size = 0,
  raster = FALSE,
  label = FALSE
) +
  ggplot2::labs(
    title = "Coarse subsets, all on one UMAP",
    caption = "All cells on umap_harmony, coloured by saved subset. No re-clustering."
  ) +
  ggplot2::theme(legend.position = "right")

p_all <- ggrastr::rasterise(p_all, layers = "Point", dpi = 300)

save_figure(p_all, out_dir, "coarse_subsets_all", width = 12, height = 12)

out_pdf_all <- file.path(out_dir, "coarse_subsets_all.pdf")

ggplot2::ggsave(out_pdf_all, plot = p_all, width = 12, height = 12,
  units = "in", device = grDevices::pdf, compress = TRUE,
  useDingbats = TRUE, bg = "white")

message("  Saved: ", out_pdf_all, " (", file.size(out_pdf_all), " bytes)")

stopifnot(
  "Input file changed during run" =
    file.exists(input) && file.mtime(input) == in_mtime &&
      file.size(input) == in_size
)
message("Combined subset PNG/PDF completed: ", out_dir)
