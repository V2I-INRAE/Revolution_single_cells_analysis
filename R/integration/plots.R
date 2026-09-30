suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggprism)
  library(rlang)
})

# plot_grouping() is shared with the existing PCA plots in norm_feat/plots.R.
plot_integration_umaps <- function(
  seurat_obj,
  group_by = c("sample", "pressure", "time", "pressure_time"),
  output_dir = NULL,
  seed = 1234,
  filename = NULL,
  facet_samples = FALSE
) {
  group_by <- match.arg(group_by)
  stopifnot("Sample faceting requires group_by = 'sample'" = !facet_samples || group_by == "sample")
  grouping <- plot_grouping(seurat_obj[[]], group_by)
  plot_obj <- seurat_obj
  plot_obj$integration_group <- grouping$values
  reductions <- c(Unintegrated = "umap", scVI = "umap_scvi", Harmony = "umap_harmony")
  stopifnot(all(reductions %in% Reductions(seurat_obj)))
  for (reduction in reductions) {
    coordinates <- Embeddings(seurat_obj, reduction)
    stopifnot(setequal(rownames(coordinates), colnames(seurat_obj)),
      ncol(coordinates) == 2L, all(is.finite(coordinates)))
  }
  if (group_by == "sample") {
    # Maximin CIELAB selection from the scientific-visualization palette,
    # excluding very pale colors (L* >= 80). Fixed mapping across all methods.
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
    # Raster sizes are rounded by plotthis: 3/1500 is 20% larger than 2/1200.
    args$raster_dpi <- c(1500, 1500)
    args$pt_size <- 3
    args$pt_alpha <- 1
    args$highlight <- TRUE
    args$highlight_color <- "black"
    args$highlight_size <- 3
    args$highlight_stroke <- 1
    args$highlight_alpha <- 1
    if (group_by == "pressure") {
      # Native faceting cannot hide Control's facet while retaining its cells.
      # Two native highlight calls keep every cell as context in both columns.
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
        for (extension in c("png", "svg")) {
          ggsave(file.path(output_dir, paste0(stem, ".", extension)), panels[[method]],
            width = 25, height = 5 * ceiling(length(categories) / 5), dpi = 300, bg = "white")
        }
      }
    }
    return(invisible(panels))
  }
  p <- patchwork::wrap_plots(panels, ncol = if (group_by == "sample") 3 else 1) +
    patchwork::plot_layout(guides = "collect")
  p <- p & theme(legend.position = if (group_by == "sample") "bottom" else "none")
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (is.null(filename)) filename <- paste0("umap_", group_by, "_lognorm")
    width <- switch(group_by, sample = 21, pressure = 10, time = 15, pressure_time = 25)
    height <- if (group_by == "sample") 9 else 15
    for (extension in c("png", "svg")) {
      ggsave(file.path(output_dir, paste0(tools::file_path_sans_ext(filename), ".", extension)), p,
        width = width, height = height, dpi = 300, bg = "white")
    }
  }
  invisible(p)
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
