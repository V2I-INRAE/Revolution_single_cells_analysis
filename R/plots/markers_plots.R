suppressPackageStartupMessages(library(ggplot2))

save_marker_plot <- function(plot, stem, directory, width, height, extension = "png") {
  path <- file.path(directory, paste0(stem, ".", extension))
  ggsave(path, plot, width = width, height = height, dpi = 300, bg = "white")
}

plot_markers_dotplot <- function(alldata, top5, sel.clust, directory) {
  dots <- DotPlot(alldata,
    features = rev(as.character(unique(top5$gene))),
    group.by = sel.clust, assay = "RNA", cols = c("#FFF7EC", "#B2182B")
  ) + coord_flip()
  dots <- dots + labs(
    title = paste(sel.clust, "| top 5 markers per cluster"),
    subtitle = "Dot size: detection; colour: scaled average expression",
    caption = "All cells; RNA LogNormalize expression. Exploratory cluster markers, not condition DE.",
    x = NULL, y = "Cluster"
  ) + theme(
    axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 9)
  )
  save_marker_plot(dots, "dotplot_top5", directory,
    width = max(12, nlevels(Idents(alldata)) * 0.35),
    height = max(6, length(unique(top5$gene)) * 0.22 + 2)
  )
}

plot_marker_features <- function(
  alldata, features, sel.clust, reduction, directory, label
) {

  present <- features[features %in% rownames(alldata[["RNA"]])]
  missing <- setdiff(features, present)
  if (length(missing) > 0L) {
    warning(
      label, ": not in object, skipped: ", paste(missing, collapse = ", ")
    )
  }
  if (length(present) == 0L) {
    return(invisible(NULL))
  }

  directory <- file.path(directory, label)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  caption <- paste(strwrap(paste(
    "All cells; RNA LogNormalize data on", reduction,
    "(no z-scores); viridis, low purple to high yellow; zeros in grey.",
    "Unavailable sample features also plot as zero after layer joining.",
    "Configured panels, not cluster annotations."
  ), width = 65), collapse = "\n")

  # RNA data is already log-normalized. Do not log it again or use scale.data.
  for (gene in present) {
    plot <- scCustomize::FeaturePlot_scCustom(
      seurat_object = alldata,
      features = gene,
      reduction = reduction,
      layer = "data",
      colors_use = viridisLite::viridis(256),
      pt.size = 0.1,
      stroke.size = 0,
      raster = FALSE,
      order = TRUE,
      combine = FALSE
    )[[1L]] + labs(
      subtitle = paste(strwrap(paste(label, "|", sel.clust), width = 45), collapse = "\n"),
      caption = caption
    ) + theme(
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(size = 7, hjust = 0)
    )
    for (extension in c("png", "svg")) {
      save_marker_plot(plot, gene, directory, width = 5.2, height = 5.8,
        extension = extension)
    }
  }
}
