# Coarse-subset UMAPs: one panel per subset, wrapped into one PNG and PDF.
# Usage: Rscript R/plots/plots_subsets.R <annotated_clustered.rds>
#
# Read-only: the object is loaded for plotting only; nothing is written back.
# Identity handling: Cluster_Highlight_Plot matches `cluster_name` against the
# active identity, so identities are set in memory to the saved subset column.
# The file on disk is never modified.
#
# PDF size note: 8 panels x 256k cells as vector circles (~14 MB) crashes
# vector editors, so point layers are rasterized at 300 dpi via ggrastr
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
existing <- list.files(out_dir, pattern = "^coarse_subsets\\.(png|pdf)$")
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
message("Plotting subset panels: ", paste(subsets, collapse = ", "))

set.seed(seed)
panels <- lapply(seq_along(subsets), function(i) {
  panel <- scCustomize::Cluster_Highlight_Plot(seurat_object = obj,
    cluster_name = subsets[[i]], highlight_color = viz_palette[[i]],
    background_color = "lightgray", reduction = reduction,
    pt.size = 0.1, stroke.size = 0, raster = FALSE, label = FALSE) +
    ggplot2::labs(title = subsets[[i]]) +
    ggplot2::theme(legend.position = "none")
  # Rasterize point layers only (300 dpi); labels, axes, legends stay
  # vector so the PDF stays editable and light.
  ggrastr::rasterise(panel, layers = "Point", dpi = 300)
})
p <- patchwork::wrap_plots(panels, ncol = 3) +
  patchwork::plot_annotation(title = "Coarse subsets",
    caption = "All cells on umap_harmony; each panel highlights one saved subset, with all other cells grey. No re-clustering.")

height <- 4 * ceiling(length(subsets) / 3)
save_figure(p, out_dir, "coarse_subsets", width = 12, height = height)
out_pdf <- file.path(out_dir, "coarse_subsets.pdf")
# Compressed vector circles use the standard Dingbats font to reduce file size.
ggplot2::ggsave(out_pdf, plot = p, width = 12, height = height,
  units = "in", device = grDevices::pdf, compress = TRUE,
  useDingbats = TRUE, bg = "white")
message("  Saved: ", out_pdf, " (", file.size(out_pdf), " bytes)")

stopifnot(
  "Input file changed during run" =
    file.exists(input) && file.mtime(input) == in_mtime &&
      file.size(input) == in_size
)
message("Combined subset PNG/PDF completed: ", out_dir)
