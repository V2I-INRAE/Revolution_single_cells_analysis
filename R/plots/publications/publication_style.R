# DEG-only final-size styling. Dimensions are physical millimetres; place at
# 100% in Affinity. Do not enlarge the canvas independently of the aesthetics.
deg_theme <- function() {
  ggplot2::theme_classic(base_size = 10, base_family = "sans") +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "black"),
      axis.text = ggplot2::element_text(size = 9, colour = "black"),
      axis.title = ggplot2::element_text(size = 10),
      axis.line = ggplot2::element_line(linewidth = 0.25),
      axis.ticks = ggplot2::element_line(linewidth = 0.25),
      axis.ticks.length = grid::unit(1.5, "mm"),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.title = ggplot2::element_text(size = 11, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9),
      plot.caption = ggplot2::element_text(size = 9, hjust = 0),
      plot.tag = ggplot2::element_text(size = 14, face = "bold"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 9),
      legend.text = ggplot2::element_text(size = 9),
      legend.key.size = grid::unit(3.5, "mm"),
      legend.spacing.x = grid::unit(1, "mm"),
      plot.margin = ggplot2::margin(3, 3, 3, 3, unit = "mm"))
}

deg_wrap <- function(x, width = 95) {
  paste(vapply(strsplit(x, "\n", fixed = TRUE)[[1]], function(line)
    paste(strwrap(line, width = width), collapse = "\n"), character(1)), collapse = "\n")
}

save_deg_figure <- function(figure, directory, stem, width = 170, height) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  # Cairo embeds fonts while keeping text as text (not outlines).
  ggplot2::ggsave(file.path(directory, paste0(stem, ".pdf")), figure,
    width = width, height = height, units = "mm", device = grDevices::cairo_pdf,
    family = "sans", bg = "white")
  save_figure(figure, directory, stem, width = width / 25.4, height = height / 25.4)
}
