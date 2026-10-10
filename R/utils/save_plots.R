# Shared plotting utilities: one save function for every figure, and the
# house colour palette used by sample/cluster visualizations.

# The visualization palette: select from it per plot (first N for samples,
# as many as needed for clusters). Never reorder: colours are part of the
# visual record of past runs.
viz_palette <- c(
  "#3A5BA0", "#F5A623", "#1ABC9C", "#90141A", "#EBA5AB", "#64B024",
  "#2E5111", "#8EC9EB", "#5C6C6B", "#D35400", "#8C6D31", "#9B59B6",
  "#E46571", "#C5A9D8", "#5A8F5A", "#8B4513", "#4E9A9A", "#6B8E23",
  "#D4753C", "#97767A", "#2E86C1", "#203161", "#BEACAE", "#C44E52",
  "#7B5EA7", "#E8A838", "#46878F", "#B07AA1", "#D98880", "#86714D",
  "#6E2F84", "#7BC657", "#708090", "#B8860B", "#40569A", "#931419",
  "#C04736", "#F09F2E", "#D83746", "#DA5E28", "#C24935", "#E2DBA4",
  "#C0D0CB", "#BDD2CB", "#DDDDDB"
)

# Save one figure as PNG: 300 dpi, white background, inches. Accepts ggplot
# objects and gtables (pheatmap output).
save_figure <- function(figure, directory, stem, width, height, dpi = 300) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(directory, paste0(stem, ".png"))
  ggplot2::ggsave(
    path, plot = figure,
    width = width, height = height, units = "in", dpi = dpi, bg = "white"
  )
  message("  Saved: ", path)
  invisible(path)
}
