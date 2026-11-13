# Apply a per-plot percentile cap (default 75th) to the values already mapped to colour.
# This affects only the display; the DotPlot data and detection sizes are unchanged.
dotplot_colour_scale <- function(plot, cap = 0.75) {
  values <- plot$data$avg.exp.scaled
  if (is.null(values) || any(!is.finite(values)) || !length(values)) {
    stop("DotPlot colour values are missing or non-finite")
  }
  upper <- unname(stats::quantile(values, cap))
  lower <- min(0, values)
  if (upper <= lower) {
    stop("DotPlot colour cap does not exceed its lower colour limit")
  }
  message(sprintf(
    "DotPlot colour distribution: min %.4f; p%g %.4f; max %.4f; capped %d/%d",
    min(values),
    100 * cap,
    upper,
    max(values),
    sum(values >= upper),
    length(values)
  ))
  ggplot2::scale_colour_gradientn(
    colours = RColorBrewer::brewer.pal(9L, "YlOrRd"),
    limits = c(lower, upper),
    oob = scales::squish
  )
}
