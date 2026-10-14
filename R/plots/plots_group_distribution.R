# Group composition figures: one 100% stacked-bar panel per grouping, showing
# the experimental-group mix within each cluster. Bar heights are
# row-normalized composition, not absolute counts: they describe sampling
# across groups, not differential abundance, and no group is an annotation.
#
# Depends on R/utils/save_plots.R (save_figure(), viz_palette), sourced by
# the caller as in every other R/plots/*.R file.

suppressPackageStartupMessages({
  library(ggplot2)
})

# Canonical display levels. REVO31-N8 (T8H) is kept as its own condition and
# time level so no sample is silently dropped.
condition_from_code <- c(
  C = "Control (C)",
  P4 = "Positive pressure (P4)",
  N4 = "Negative pressure (N4)",
  N8 = "Negative pressure (N8)",
  P10 = "Positive pressure (P10)",
  N10 = "Negative pressure (N10)"
)
condition_order <- unname(condition_from_code[c("C", "P4", "N4", "N8", "P10", "N10")])

pressure_from_value <- c(
  None = "Control (C)",
  positive = "Positive pressure",
  negative = "Negative pressure"
)
pressure_order <- unname(pressure_from_value)

time_from_value <- c(
  T0H = "Control (C or T0H)",
  T4H = "T4H",
  T8H = "T8H",
  T10H = "T10H"
)
time_order <- unname(time_from_value)

# Add the three grouping columns from the standard pipeline metadata
# (sample, pressure, time_point as set in R/io.R). Unknown codes or values
# stop loudly: a new sample suffix must get an explicit level, never a
# silent recode or drop.
derive_group_columns <- function(meta) {
  missing <- setdiff(c("sample", "pressure", "time_point"), names(meta))
  if (length(missing) > 0L) {
    stop("Metadata lacks columns: ", paste(missing, collapse = ", "))
  }
  codes <- sub("^.*-", "", as.character(meta$sample))
  unknown_codes <- setdiff(unique(codes), names(condition_from_code))
  if (length(unknown_codes) > 0L) {
    stop("Unknown condition codes in sample IDs: ",
      paste(unknown_codes, collapse = ", "))
  }
  unknown_pressure <- setdiff(
    unique(as.character(meta$pressure)), names(pressure_from_value))
  if (length(unknown_pressure) > 0L) {
    stop("Unknown pressure values: ", paste(unknown_pressure, collapse = ", "))
  }
  unknown_time <- setdiff(
    unique(as.character(meta$time_point)), names(time_from_value))
  if (length(unknown_time) > 0L) {
    stop("Unknown time_point values: ", paste(unknown_time, collapse = ", "))
  }
  meta$condition_group <- factor(unname(condition_from_code[codes]),
    levels = condition_order)
  meta$pressure_group <- factor(
    unname(pressure_from_value[as.character(meta$pressure)]),
    levels = pressure_order)
  meta$time_group <- factor(
    unname(time_from_value[as.character(meta$time_point)]),
    levels = time_order)
  meta
}

# One 100% stacked-bar plot: clusters on x, group mix on y. Full canonical
# levels are kept (drop = FALSE) so colours stay comparable across runs.
# The caption is wrapped to fit the figure width (long single-line captions
# would run past the plot margin and clip).
plot_group_stacked <- function(meta, cluster_col, group_col,
                               group_title, legend_title, caption) {
  caption <- paste(strwrap(caption, width = 120), collapse = "\n")
  stopifnot(
    "Cluster column missing from metadata" = cluster_col %in% names(meta),
    "Group column missing from metadata; run derive_group_columns() first" =
      group_col %in% names(meta)
  )
  groups <- levels(meta[[group_col]])
  colors <- stats::setNames(viz_palette[seq_along(groups)], groups)
  plot_data <- data.frame(
    cluster = factor(meta[[cluster_col]]),
    group = factor(meta[[group_col]], levels = groups)
  )
  ggplot(plot_data, aes(x = cluster, fill = group)) +
    geom_bar(position = "fill") +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    scale_fill_manual(values = colors, drop = FALSE) +
    labs(title = group_title, x = "Cluster", y = "Within-cluster proportion",
      fill = legend_title, caption = caption) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
      plot.caption = element_text(hjust = 0))
}

plot_condition_distribution <- function(meta, cluster_col) {
  plot_group_stacked(meta, cluster_col, "condition_group",
    group_title = "Cluster composition by condition",
    legend_title = "Condition",
    caption = paste(
      "100% stacked bars: condition mix within each cluster (row-normalized).",
      "REVO31-N8 kept as its own level.",
      "Composition describes sampling, not differential abundance."))
}

plot_pressure_distribution <- function(meta, cluster_col) {
  plot_group_stacked(meta, cluster_col, "pressure_group",
    group_title = "Cluster composition by pressure",
    legend_title = "Pressure",
    caption = paste(
      "100% stacked bars: pressure mix within each cluster (row-normalized).",
      "Control is unventilated (C); positive/negative pool all time points.",
      "Composition describes sampling, not differential abundance."))
}

plot_time_distribution <- function(meta, cluster_col) {
  plot_group_stacked(meta, cluster_col, "time_group",
    group_title = "Cluster composition by time point",
    legend_title = "Time point",
    caption = paste(
      "100% stacked bars: time-point mix within each cluster (row-normalized).",
      "T8H comes from REVO31-N8 only and is kept as its own level.",
      "Composition describes sampling, not differential abundance."))
}

# Full set: derive groups from the object metadata, render the three panels,
# save them as PNGs. Returns the saved paths invisibly.
plot_group_distributions <- function(obj, cluster_col = "harmony_snn_res.0.7",
                                     plot_dir) {
  meta <- derive_group_columns(obj[[]])
  plots <- list(
    condition = plot_condition_distribution(meta, cluster_col),
    pressure = plot_pressure_distribution(meta, cluster_col),
    time = plot_time_distribution(meta, cluster_col)
  )
  stems <- c(condition = "group_condition_per_cluster",
    pressure = "group_pressure_per_cluster",
    time = "group_time_per_cluster")
  width <- max(10, length(unique(meta[[cluster_col]])) * 0.5)
  paths <- mapply(function(panel, stem) {
    save_figure(panel, plot_dir, stem, width = width, height = 6)
  }, plots, stems, USE.NAMES = TRUE)
  invisible(paths)
}
