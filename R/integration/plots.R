suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggprism)
  library(rlang)
})

plot_grouping <- function(
  metadata,
  group_by = c("sample", "pressure", "time", "pressure_time")
) {
  group_by <- match.arg(group_by)
  required <- switch(
    group_by,
    sample = "sample",
    pressure = "pressure",
    time = "time_point",
    pressure_time = c("pressure", "time_point")
  )
  missing_columns <- setdiff(required, colnames(metadata))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing plotting metadata: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  pressure <- as.character(metadata$pressure)
  time <- as.character(metadata$time_point)
  if (group_by != "sample") {
    unknown_pressure <- setdiff(
      unique(pressure),
      c("None", "C", "positive", "negative")
    )
    if (length(unknown_pressure) > 0L || anyNA(pressure)) {
      stop("Unexpected or missing pressure values in PCA metadata.")
    }
    pressure[pressure == "None"] <- "C"
  }
  if (group_by %in% c("time", "pressure_time")) {
    unknown_time <- setdiff(unique(time), c("T0H", "T4H", "T8H", "T10H"))
    if (length(unknown_time) > 0L || anyNA(time)) {
      stop("Unexpected or missing time values in PCA metadata.")
    }
    time[time == "T8H"] <- "T10H"
  }

  if (group_by == "sample") {
    values <- as.character(metadata$sample)
    if (anyNA(values) || any(values == "")) {
      stop("Missing sample values in PCA metadata.")
    }
    values <- factor(values, levels = unique(values))
    legend_title <- "Sample"
  } else if (group_by == "pressure") {
    values <- factor(pressure, levels = c("C", "positive", "negative"))
    legend_title <- "Pressure"
  } else if (group_by == "time") {
    values <- factor(time, levels = c("T0H", "T4H", "T10H"))
    legend_title <- "Time"
  } else {
    if (any(time != "T0H" & pressure == "C")) {
      stop("Control pressure is only expected at T0H.")
    }
    values <- ifelse(time == "T0H", "C", paste(pressure, time, sep = "_"))
    values <- factor(
      values,
      levels = c(
        "C", "positive_T4H", "negative_T4H",
        "positive_T10H", "negative_T10H"
      )
    )
    legend_title <- "Pressure and time"
  }
  if (anyNA(values)) {
    stop("PCA grouping produced missing labels.")
  }

  list(values = values, title = legend_title)
}

plot_integration_umaps <- function(
  seurat_obj,
  group_by = c("sample", "pressure", "time", "pressure_time"),
  output_dir = NULL,
  seed = 1234,
  filename = NULL,
  facet_samples = FALSE,
  method_reduction = NULL,
  reductions = NULL
) {
  group_by <- match.arg(group_by)
  stopifnot("Sample faceting requires group_by = 'sample'" = !facet_samples || group_by == "sample")
  grouping <- plot_grouping(seurat_obj[[]], group_by)
  plot_obj <- seurat_obj
  plot_obj$integration_group <- grouping$values
  # By default, plot the checkpoint's own method beside its unintegrated baseline.
  method <- seurat_obj@misc$integration$method
  stopifnot(
    "Expected an integration checkpoint with a declared method" =
      method %in% c("scvi", "harmony", "cca")
  )
  method_label <- switch(method, scvi = "scVI", harmony = "Harmony", cca = "CCA")
  if (is.null(method_reduction)) method_reduction <- paste0("umap_", method)
  if (is.null(reductions)) {
    reductions <- setNames(
      c("umap", method_reduction),
      c("Unintegrated", method_label)
    )
  }
  if (!all(reductions %in% Reductions(seurat_obj))) {
    stop("Checkpoint is missing its saved UMAP reductions: ",
      paste(setdiff(reductions, Reductions(seurat_obj)), collapse = ", "))
  }
  for (reduction in reductions) {
    coordinates <- Embeddings(seurat_obj, reduction)
    stopifnot(setequal(rownames(coordinates), colnames(seurat_obj)),
      ncol(coordinates) == 2L, all(is.finite(coordinates)))
  }
  if (group_by == "sample") {
    palette <- c(
      "#3A5BA0", "#F5A623", "#1ABC9C", "#90141A", "#EBA5AB", "#64B024",
      "#2E5111", "#8EC9EB", "#5C6C6B", "#D35400", "#8C6D31", "#9B59B6",
      "#E46571", "#C5A9D8", "#5A8F5A", "#8B4513", "#4E9A9A", "#6B8E23",
      "#D4753C", "#97767A", "#2E86C1", "#203161", "#BEACAE"
    )
    stopifnot(nlevels(grouping$values) <= length(palette))
  } else {
    levels(plot_obj$integration_group) <- switch(group_by,
      pressure = c("Control", "Positive", "Negative"),
      time = c("0 h", "4 h", "10 h"),
      pressure_time = c("Control", "Positive 4 h", "Negative 4 h", "Positive 10 h", "Negative 10 h")
    )
    palette <- switch(group_by,
      pressure = c("#5A8F5A", "#D4753C", "#3A5BA0"),
      time = c("#5A8F5A", "#D4753C", "#7B5EA7"),
      pressure_time = c("#5A8F5A", "#D4753C", "#3A5BA0", "#C44E52", "#7B5EA7")
    )
  }
  categories <- levels(plot_obj$integration_group)
  stopifnot("Every requested category must have cells" = all(table(plot_obj$integration_group) > 0))
  colours <- setNames(palette[seq_along(categories)], categories)
  panels <- lapply(names(reductions), function(method) {
    args <- list(object = plot_obj, reduction = reductions[[method]],
      group_by = "integration_group", palcolor = colours,
      theme = "theme_blank", theme_args = list(base_size = 13),
      title = method, seed = seed, show_stat = FALSE, label = FALSE,
      xlab = "UMAP 1", ylab = "UMAP 2",
      raster = TRUE, raster_dpi = c(1200, 1200), pt_size = 2, pt_alpha = 0.5,
      bg_color = "#707070", order = "random", legend.position = "none")
    if (group_by == "sample" && !facet_samples) {
      args$legend.position <- "bottom"
      return(do.call(scplotter::CellDimPlot, args) + labs(colour = "Sample") +
        guides(colour = guide_legend(ncol = 6, override.aes = list(size = 2.5, alpha = 1))))
    }
    args$raster_dpi <- c(1500, 1500)
    args$pt_size <- 3
    args$pt_alpha <- 1
    args$highlight <- TRUE
    args$highlight_color <- "black"
    args$highlight_size <- 3
    args$highlight_stroke <- 1
    args$highlight_alpha <- 1
    if (group_by == "pressure") {
      args$pt_alpha <- 0.5
      pressure_panels <- lapply(c("Positive", "Negative"), function(category) {
        args$highlight <- sprintf('integration_group == "%s"', category)
        args$palcolor <- setNames(rep("#707070", length(categories)), categories)
        args$palcolor[category] <- colours[category]
        args$title <- paste(method, "\u2014", category)
        do.call(scplotter::CellDimPlot, args)
      })
      return(patchwork::wrap_plots(pressure_panels, ncol = 2))
    }
    args$facet_by <- "integration_group"
    args$facet_ncol <- if (facet_samples) 5 else length(categories)
    p <- do.call(scplotter::CellDimPlot, args)
    if (facet_samples) p <- p + theme(plot.margin = margin(12, 18, 24, 18))
    p
  })
  if (facet_samples) {
    names(panels) <- names(reductions)
    if (!is.null(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      if (is.null(filename)) filename <- "umap_sample_facets"
      for (method in names(panels)) {
        stem <- paste0(tools::file_path_sans_ext(filename), "_", tolower(method), "_lognorm")
        ggsave(file.path(output_dir, paste0(stem, ".png")), panels[[method]],
          width = 25, height = 5 * ceiling(length(categories) / 5), dpi = 300, bg = "white")
      }
    }
    return(invisible(panels))
  }
  p <- patchwork::wrap_plots(panels,
    ncol = if (group_by == "sample") length(reductions) else 1) +
    patchwork::plot_layout(guides = "collect")
  p <- p & theme(legend.position = if (group_by == "sample") "bottom" else "none")
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (is.null(filename)) filename <- paste0("umap_", group_by, "_lognorm")
    width <- switch(group_by,
      sample = 7 * length(reductions), pressure = 10, time = 15, pressure_time = 25
    )
    height <- if (group_by == "sample") 9 else 3.75 * length(reductions)
    ggsave(file.path(output_dir, paste0(tools::file_path_sans_ext(filename), ".png")), p,
      width = width, height = height, dpi = 300, bg = "white")
  }
  invisible(p)
}

plot_integration_figures <- function(seurat_obj, output_dir, seed = 1234,
  method_reduction = NULL) {
  for (group_by in c("sample", "pressure", "time", "pressure_time")) {
    plot_integration_umaps(seurat_obj, group_by = group_by,
      output_dir = output_dir, seed = seed, method_reduction = method_reduction)
  }
  plot_integration_umaps(seurat_obj, group_by = "sample", facet_samples = TRUE,
    output_dir = output_dir, seed = seed, method_reduction = method_reduction)

  method <- seurat_obj@misc$integration$method
  expected <- c(
    paste0("umap_", c("sample", "pressure", "time", "pressure_time"), "_lognorm.png"),
    paste0("umap_sample_facets_", c("unintegrated", method), "_lognorm.png")
  )
  figures <- list.files(output_dir, pattern = "[.]png$")
  stopifnot(setequal(figures, expected))
  invisible(figures)
}

plot_integration_metrics <- function(
  diagnostics, output_dir = NULL, filename = "integration_metrics.png",
  width = 12, height = 5
) {
  scores <- diagnostics$scores
  scores$method <- factor(scores$reduction,
    levels = c("pca", "harmony", "integrated_cca", "integrated_scvi"),
    labels = c("Unintegrated", "Harmony", "CCA", "scVI")
  )
  ilisi <- ggplot(
    scores, aes(x = .data$method, y = .data$ilisi, fill = .data$method)
  ) +
    geom_violin(show.legend = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, show.legend = FALSE) +
    theme_prism() +
    labs(
      x = NULL, y = "Raw sample iLISI",
      title = "Local sample mixing"
    )
  silhouette <- ggplot(
    scores, aes(
      x = .data$method, y = .data$batch_silhouette,
      fill = .data$method
    )
  ) +
    geom_violin(show.legend = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, show.legend = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_prism() +
    labs(
      x = NULL, y = "Raw sample silhouette width",
      title = "Sample separation (mixing near zero)"
    )
  p <- patchwork::wrap_plots(ilisi, silhouette, ncol = 2) +
    patchwork::plot_annotation(
      caption = paste(
        "Shared sampled cells; sample mixing is not proof",
        "of biological preservation."
      )
    )
  if (!is.null(output_dir)) {
    ggsave(file.path(output_dir, filename), p,
      width = width, height = height, dpi = 300
    )
  }
  invisible(p)
}
