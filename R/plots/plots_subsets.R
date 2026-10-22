# Faceted UMAPs from a saved annotated object: 12 subset panels and 24
# cluster panels (res 0.7). PNG and PDF.
# Usage: Rscript R/plots/plots_subsets.R <data/job-XXX/annotated_clustered.rds>
#
# Read-only: the object is loaded for plotting only; nothing is written back.
# Identity handling: Cluster_Highlight_Plot matches `cluster_name` against the
# active identity, so identities are set in memory to the saved columns.
# The file on disk is never modified.
#
# PDF size note: large cell layers as vector circles can overwhelm vector
# editors, so point layers are rasterized at 300 dpi via ggrastr
# (text, axes and legends stay vector). PNGs are rendered by save_figure;
# the cluster PDF uses a multipage pdf() device loop.
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
cluster_col <- "harmony_snn_res.0.7"
seed <- 1234L # production run seed (results/job-44493014/.INFO)

out_dir <- file.path("results", basename(dirname(input)),
  "clusters_plots", "subset")

if (!file.exists(input)) stop("Input not found: ", input)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}
planned <- c("coarse_subsets_faceted.png", "coarse_subsets_faceted.pdf",
  "clusters_res0.7_faceted_p1.png", "clusters_res0.7_faceted_p2.png",
  "clusters_res0.7_faceted.pdf")
existing <- planned[file.exists(file.path(out_dir, planned))]
if (length(existing) > 0L) {
  stop("Plots already exist in ", out_dir, ": ",
    paste(existing, collapse = ", "),
    "; preserve them before running again.")
}
in_mtime <- file.mtime(input)
in_size <- file.size(input)

message("Loading (read-only): ", input)
obj <- readRDS(input)
stopifnot(
  "Saved UMAP missing" = reduction %in% SeuratObject::Reductions(obj),
  "subset column missing" = "subset" %in% names(obj[[]]),
  "cluster column missing" = cluster_col %in% names(obj[[]])
)

# One shared panel: highlight a single group, everything else grey.
highlight_panel <- function(group, color) {
  panel <- scCustomize::Cluster_Highlight_Plot(
    seurat_object = obj,
    cluster_name = group,
    highlight_color = color,
    background_color = "lightgray",
    reduction = reduction,
    pt.size = 0.1,
    stroke.size = 0,
    raster = FALSE,
    label = FALSE,
    seed = seed
  ) +
    ggplot2::labs(title = group) +
    ggplot2::theme(legend.position = "none")
  ggrastr::rasterise(panel, layers = "Point", dpi = 300)
}

# Subset facets (12 panels, one page).
SeuratObject::Idents(obj) <- "subset"
subsets <- levels(SeuratObject::Idents(obj))
subset_colors <- stats::setNames(viz_palette[seq_along(subsets)], subsets)
message("Plotting subset panels: ", paste(subsets, collapse = ", "))
subset_panels <- lapply(subsets, function(s) highlight_panel(s, subset_colors[[s]]))
p_subsets <- patchwork::wrap_plots(subset_panels, ncol = 3) +
  patchwork::plot_annotation(title = "Coarse subsets",
    caption = "All cells on umap_harmony; each panel highlights one saved subset, with all other cells grey. No re-clustering.")
save_figure(p_subsets, out_dir, "coarse_subsets_faceted", width = 12, height = 16)
out_pdf_sub <- file.path(out_dir, "coarse_subsets_faceted.pdf")
ggplot2::ggsave(out_pdf_sub, plot = p_subsets, width = 12, height = 16,
  units = "in", device = grDevices::pdf, compress = TRUE,
  useDingbats = TRUE, bg = "white")
message("  Saved: ", out_pdf_sub, " (", file.size(out_pdf_sub), " bytes)")

# Cluster facets (24 panels, two pages of 12; colours assigned once).
SeuratObject::Idents(obj) <- cluster_col
clusters <- levels(SeuratObject::Idents(obj))
cluster_colors <- stats::setNames(viz_palette[seq_along(clusters)], clusters)
message("Plotting cluster panels: ", paste(clusters, collapse = ", "))
cluster_panels <- lapply(clusters, function(cl) highlight_panel(cl, cluster_colors[[cl]]))
pages <- split(seq_along(cluster_panels), ceiling(seq_along(cluster_panels) / 12))
page_plots <- lapply(seq_along(pages), function(i) {
  patchwork::wrap_plots(cluster_panels[pages[[i]]], ncol = 3) +
    patchwork::plot_annotation(
      title = paste0("Clusters res 0.7 (page ", i, " of ", length(pages), ")"),
      caption = "All cells on umap_harmony; each panel highlights one saved cluster, with all other cells grey. No re-clustering.")
})
for (i in seq_along(page_plots)) {
  save_figure(page_plots[[i]], out_dir,
    paste0("clusters_res0.7_faceted_p", i), width = 12, height = 16)
}
out_pdf_cl <- file.path(out_dir, "clusters_res0.7_faceted.pdf")
grDevices::pdf(out_pdf_cl, width = 12, height = 16,
  onefile = TRUE, compress = TRUE, useDingbats = TRUE, bg = "white")
for (page in page_plots) print(page)
grDevices::dev.off()
message("  Saved: ", out_pdf_cl, " (", file.size(out_pdf_cl), " bytes)")

stopifnot(
  "Input file changed during run" =
    file.exists(input) && file.mtime(input) == in_mtime &&
      file.size(input) == in_size
)
message("Faceted subset/cluster PNG/PDF completed: ", out_dir)
