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
  width = 15,
  height = 6,
  filename = NULL
) {
  group_by <- match.arg(group_by)
  grouping <- plot_grouping(seurat_obj[[]], group_by)
  plot_obj <- seurat_obj
  plot_obj$integration_group <- grouping$values
  reductions <- c(
    Unintegrated = "umap", Harmony = "umap_harmony", CCA = "umap_cca"
  )
  reductions <- reductions[reductions %in% Reductions(seurat_obj)]
  if (length(reductions) == 0L) stop("No comparison UMAPs found in the object.")
  panels <- lapply(names(reductions), function(method) {
    DimPlot(
      plot_obj,
      reduction = reductions[[method]],
      group.by = "integration_group", shuffle = TRUE, seed = seed
    ) + NoAxes() + labs(title = method, colour = grouping$title)
  })
  p <- patchwork::wrap_plots(panels, ncol = length(reductions)) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  if (!is.null(output_dir)) {
    if (is.null(filename)) filename <- paste0("umap_", group_by, ".png")
    ggsave(file.path(output_dir, filename), p,
      width = width, height = height, dpi = 300
    )
  }
  invisible(p)
}

plot_integration_metrics <- function(
  diagnostics, output_dir = NULL, filename = "integration_metrics.png",
  width = 12, height = 5
) {
  scores <- diagnostics$scores
  scores$method <- factor(scores$reduction,
    levels = c("pca", "harmony", "integrated_cca"),
    labels = c("Unintegrated", "Harmony", "CCA")
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
