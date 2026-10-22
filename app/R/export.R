# Match the hybrid PDFs in R/plots/markers_plots.R: raster cells, vector frame.
save_explorer_pdf <- function(path, data, foreground, view, gene = NULL,
  colors = character(0)) {
  expression <- isTRUE(foreground$marker$showscale)
  # One point layer preserves the black-under-color drawing order and becomes
  # one raster image, rather than hundreds of thousands of editable circles.
  points <- data.frame(
    umap_1 = c(data$cells$umap_1, foreground$x),
    umap_2 = c(data$cells$umap_2, foreground$y),
    color = c(rep(if (expression) NA_real_ else "#000000", nrow(data$cells)),
      foreground$marker$color))
  plot <- ggplot2::ggplot(points,
    ggplot2::aes(x = umap_1, y = umap_2, colour = color)) +
    ggplot2::geom_point(size = view$size * 25.4 / 96, stroke = 0) +
    ggplot2::coord_fixed(xlim = unlist(view$ranges$umap_1),
      ylim = unlist(view$ranges$umap_2), expand = FALSE) +
    ggplot2::labs(x = "UMAP 1", y = "UMAP 2",
      title = if (is.null(gene)) "UMAP" else gene$label,
      subtitle = data$label) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank())
  if (expression) {
    plot <- plot + ggplot2::scale_colour_gradientn(
      colours = vapply(foreground$marker$colorscale, `[[`, character(1), 2L),
      values = vapply(foreground$marker$colorscale, `[[`, numeric(1), 1L),
      limits = c(foreground$marker$cmin, foreground$marker$cmax),
      na.value = "#000000", name = "RNA\nLogNormalize",
      guide = ggplot2::guide_colourbar(display = "rectangles"))
  } else {
    plot <- plot + ggplot2::scale_colour_identity(name = "Cluster",
      breaks = unname(colors), labels = names(colors),
      guide = if (length(colors)) "legend" else "none")
  }
  ggplot2::ggsave(path,
    plot = ggrastr::rasterise(plot, layers = "Point", dpi = 600),
    device = grDevices::pdf, width = view$width / 96,
    height = view$height / 96, units = "in", bg = "white",
    useDingbats = FALSE)
}
