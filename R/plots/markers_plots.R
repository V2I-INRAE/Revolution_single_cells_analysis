suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# Paginated per-panel PDF layout: square 12 x 12 in pages,
# 2 columns x 2 rows = 4 gene tiles per page. Titles stay selectable.
panel_pdf_ncol <- 2L
panel_pdf_nrow <- 2L
panel_pdf_genes_per_page <- panel_pdf_ncol * panel_pdf_nrow
# Raster resolution of the embedded cell layer. Frame elements (titles,
# axes, legend text) stay vector; only the point layer is rasterized via
# ggrastr::rasterise() below. Retain 600 dpi for print and Illustrator.
panel_pdf_raster_dpi <- 600
# Retain the existing format-specific point sizes.
marker_png_pt.size <- 0.1
panel_pdf_pt.size <- 0.04

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
  side <- max(12, nlevels(Idents(alldata)) * 0.35,
    length(unique(top5$gene)) * 0.22 + 2)
  for (extension in c("png", "pdf")) {
    save_marker_plot(dots, "dotplot_top5", directory,
      width = side, height = side, extension = extension)
  }
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
    "(no z-scores); viridis, zero in low purple to high yellow.",
    "Unavailable sample features also plot as zero after layer joining.",
    "Configured panels, not cluster annotations."
  ), width = 65), collapse = "\n")

  # RNA data is already log-normalized. Do not log it again or use scale.data.
  #
  # Per-gene square PNGs are written one file per gene. The panel PDF is
  # one multipage square file: 4 gene tiles per page (2 x 2 grid).
  # Cell layers are rasterized at panel_pdf_raster_dpi (600) so the PDF
  # stays small and Illustrator-friendly; titles, axes and legend
  # text stay vector. The page header/footer carry the panel subtitle and
  # methods caption once per page instead of repeating them on every tile.
  page_header <- paste(label, "|", sel.clust)
  page_footer <- paste(
    "All cells; RNA LogNormalize data on", reduction,
    "(no z-scores); viridis, zero in low purple to high yellow.",
    "Unavailable sample features also plot as zero after layer joining.",
    "Configured panels, not cluster annotations."
  )
  # Print one full or partial page to the panel PDF device. Drawing goes
  # through patchworkGrob + grid.draw (not print.patchwork) so rendering
  # errors propagate instead of being swallowed by patchwork 1.3.2's
  # print method, which catches grid.draw() errors without rethrowing.
  print_panel_page <- function(plots) {
    page <- patchwork::wrap_plots(plots,
      ncol = panel_pdf_ncol, nrow = panel_pdf_nrow) +
      patchwork::plot_annotation(
        title = page_header,
        caption = paste(strwrap(page_footer, width = 110), collapse = "\n"),
        theme = theme(
          plot.title = element_text(size = 12, face = "bold"),
          plot.caption = element_text(size = 7, hjust = 0))
      )
    grid::grid.newpage()
    grid::grid.draw(patchwork::patchworkGrob(page))
  }

  pdf_path <- file.path(directory, paste0(label, ".pdf"))

  grDevices::pdf(pdf_path, width = 12, height = 12, bg = "white")

  pdf_device <- grDevices::dev.cur()

  on.exit(grDevices::dev.off(pdf_device), add = TRUE)

  page_plots <- list()

  for (gene in present) {
    plot <- scCustomize::FeaturePlot_scCustom(
      seurat_object = alldata,
      features = gene,
      reduction = reduction,
      layer = "data",
      colors_use = viridisLite::viridis(256),
      na_cutoff = 0,
      pt.size = marker_png_pt.size,
      stroke.size = 0,
      raster = FALSE,
      order = TRUE,
      combine = FALSE
    )[[1L]]
    # Keep the panel|cluster subtitle and methods caption on each PNG.
    save_marker_plot(plot + labs(
      subtitle = paste(
        strwrap(paste(label, "|", sel.clust), width = 45),
        collapse = "\n"),
      caption = caption
    ) + theme(
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(size = 7, hjust = 0)
    ), gene, directory, width = 5.2, height = 5.2, extension = "png")
    # Separate PDF tile copy: gene title only (page header/footer carry the
    # panel context once per page). display="rectangles" renders the
    # colourbar as vector coloured rectangles with vector text, so the whole
    # legend can be removed or restyled in Illustrator.
    pdf_plot <- plot + labs(
      title = gene,
      subtitle = NULL,
      caption = NULL
    ) + theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5)
    ) + guides(colour = guide_colourbar(display = "rectangles"))
    # Apply the existing PDF-only point size; the PNG above is unaffected.
    # The base plot uses a constant size, hence aes_params.
    pdf_point_scale <- panel_pdf_pt.size / marker_png_pt.size

    for (i in seq_along(pdf_plot$layers)) {
      if (inherits(pdf_plot$layers[[i]]$geom, "GeomPoint") &&
        !is.null(pdf_plot$layers[[i]]$aes_params$size)) {
        pdf_plot$layers[[i]]$aes_params$size <-
          pdf_plot$layers[[i]]$aes_params$size * pdf_point_scale
      }
    }

    # Rasterize only the point layer; frame elements stay vector (Seurat
    # 5.5.1's own raster=TRUE path ignores the dpi setting, so rasterising
    # the vector plot explicitly preserves geometry at panel_pdf_raster_dpi).
    page_plots[[gene]] <- ggrastr::rasterise(pdf_plot,
      layers = "Point", dpi = panel_pdf_raster_dpi)

    if (length(page_plots) == panel_pdf_genes_per_page) {
      print_panel_page(page_plots)
      page_plots <- list()
    }
  }

  if (length(page_plots) > 0L) {
    print_panel_page(page_plots)
  }
}
