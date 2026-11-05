# Replot a saved marker panel locally without rerunning marker discovery.
# Usage: Rscript R/clusters/replot_marker_dotplot.R <object.rds> <markers_top5.csv> <cluster-column> <parent|subset> <output-dir>
suppressPackageStartupMessages(library(Seurat))
source("R/config.R")
source("R/cell_cycle.R")
source("R/plots/dotplot_colour.R")
run_replot <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 5L) stop("Expected object, marker CSV, cluster column, parent|subset, output directory")
  input <- args[[1L]]; panel <- args[[2L]]; cluster <- args[[3L]]; kind <- args[[4L]]; out <- args[[5L]]
  if (!kind %in% c("parent", "subset")) stop("Invalid plot kind")
  if (!file.exists(input) || !file.exists(panel)) stop("Input object or marker table missing")
  if (dir.exists(out)) stop("Output directory exists; refusing to overwrite: ", out)
  dir.create(out, recursive = TRUE)
  info_path <- file.path(out, ".INFO")
  info <- c("# Marker dotplot replot", "Status: running", paste("Started:", Sys.time()),
    paste("Input object (read-only):", input), paste("Saved marker panel (read-only):", panel),
    paste("Cluster column:", cluster), paste("Plot kind:", kind),
    paste("RNA LogNormalize; Seurat DotPlot; gene-wise scaling:", kind == "parent"),
    "YlOrRd; per-plot 75th percentile cap; dots at or above the cap share the maximum colour.",
    "See execution log for distribution and cap. No marker discovery or object writes.")
  writeLines(info, info_path)
  completed <- FALSE
  on.exit({
    info[info == "Status: running"] <- if (completed) "Status: completed" else "Status: failed/incomplete; see log"
    writeLines(c(info, paste("Finished:", Sys.time())), info_path)
  }, add = TRUE)
  before <- file.info(c(input, panel))[, c("size", "mtime")]
  top5 <- read.csv(panel, stringsAsFactors = FALSE)
  obj <- readRDS(input)
  metadata <- obj[[]]
  if (!cluster %in% names(metadata) || anyNA(metadata[[cluster]]) ||
      !identical(rownames(metadata), colnames(obj)) || anyDuplicated(colnames(obj))) stop("Invalid cluster metadata")
  DefaultAssay(obj) <- "RNA"
  obj <- JoinLayers(obj, assay = "RNA", layers = "^data($|[.])", new = "data")
  obj <- rename_cc_features(obj, config$markers$feature_renames)
  missing <- setdiff(unique(top5$gene), rownames(obj[["RNA"]]))
  if (length(missing)) stop("Saved marker features missing: ", paste(missing, collapse = ", "))
  Idents(obj) <- obj[[cluster, drop = TRUE]]
  if (!setequal(as.character(unique(top5$cluster)), levels(Idents(obj)))) stop("Marker and object clusters differ")
  genes <- rev(as.character(unique(top5$gene)))
  dots <- Seurat::DotPlot(obj, features = genes, group.by = cluster, assay = "RNA", scale = kind == "parent")
  dat <- dots$data
  idx <- match(paste(top5$cluster, top5$gene), paste(dat$id, dat$features.plot))
  if (anyNA(idx) || any(abs(dat$pct.exp[idx] / 100 - top5$pct.1) > 0.0011)) stop("Detection fractions disagree with saved markers")
  dots <- dots + dotplot_colour_scale(dots) + ggplot2::coord_flip() + ggplot2::labs(
    title = paste(cluster, "| top 5 markers per cluster"),
    subtitle = if (kind == "parent") "Dot size: detection; colour: gene-wise scaled average expression" else "Dot size: detection; colour: normalized average expression",
    caption = "Colour capped at plot 75th percentile. Exploratory cluster markers, not condition DE.", x = NULL, y = "Cluster") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9), axis.text.x = ggplot2::element_text(size = 9))
  side <- max(12, nlevels(Idents(obj)) * 0.35, length(genes) * 0.22 + 2)
  pdf <- file.path(out, "dotplot_top5.pdf"); png <- file.path(out, "dotplot_top5.png")
  ggplot2::ggsave(pdf, dots, width = side, height = side, units = "in", bg = "white")
  ggplot2::ggsave(png, dots, width = side, height = side, units = "in", dpi = 300, bg = "white")
  stopifnot(identical(rawToChar(readBin(pdf, "raw", 5L)), "%PDF-"),
    all(dim(png::readPNG(png, native = TRUE)) > 0),
    identical(before, file.info(c(input, panel))[, c("size", "mtime")]))
  info <- c(info, paste("Validation: PDF header and PNG reopened; inputs unchanged; detection agrees with marker table;", ncol(obj), "cells;", length(levels(Idents(obj))), "clusters;", length(genes), "marker genes."))
  completed <- TRUE
}
run_replot()
