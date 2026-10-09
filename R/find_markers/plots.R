suppressPackageStartupMessages(library(ggplot2))

marker_colours <- c(
  "#3A5BA0", "#F5A623", "#1ABC9C", "#90141A", "#EBA5AB", "#64B024",
  "#2E5111", "#8EC9EB", "#5C6C6B", "#D35400", "#8C6D31", "#9B59B6",
  "#E46571", "#C5A9D8", "#5A8F5A", "#8B4513", "#4E9A9A", "#6B8E23",
  "#D4753C", "#97767A", "#2E86C1", "#203161", "#BEACAE", "#C44E52",
  "#7B5EA7", "#E8A838", "#46878F", "#B07AA1", "#D98880", "#86714D",
  "#6E2F84", "#7BC657", "#708090", "#B8860B", "#40569A", "#931419"
)

save_marker_plot <- function(plot, stem, directory, width, height) {
  ggsave(file.path(directory, paste0(stem, ".png")), plot,
    width = width, height = height, dpi = 300, bg = "white"
  )
}

plot_marker_bars <- function(markers, clusters, resolution, directory) {
  pages <- split(clusters, ceiling(seq_along(clusters) / 6))
  xmax <- max(c(0.25, markers$avg_log2FC)) * 1.05
  for (page in seq_along(pages)) {
    panel_clusters <- pages[[page]]
    rows <- ceiling(length(panel_clusters) / 3)
    height <- rows * 5.5
    path <- file.path(directory, sprintf("barplot_top25_%02d.png", page))
    png(path, width = 21, height = height, units = "in", res = 300, type = "cairo")
    tryCatch(
      {
        # Wide left margin keeps long gene labels clear of the neighbouring panel;
        # a roomy top margin keeps each panel title clear of the row above.
        par(
          mfrow = c(rows, 3), cex = 1, mar = c(2.2, 11, 4.2, 1.5),
          oma = c(4.4, 0, 3.6, 0), mgp = c(2.2, 0.6, 0)
        )
        for (cluster in panel_clusters) {
          selected <- markers[as.character(markers$cluster) == cluster, ]
          if (nrow(selected) == 0L) {
            plot.new()
            title(main = paste(cluster, "vs. rest"))
            text(0.5, 0.5, "No significant markers")
          } else {
            barplot(sort(setNames(selected$avg_log2FC, selected$gene)),
              horiz = TRUE, las = 1, main = paste(cluster, "vs. rest"),
              col = "grey65", border = "white", yaxs = "i", cex.names = 0.7,
              xlim = c(0, xmax)
            )
            abline(v = c(0, 0.25), lty = c(1, 2))
          }
        }
        mtext(paste("CCA resolution", resolution, "| all cells | top positive markers"),
          outer = TRUE, side = 3, line = 1
        )
        # One outer axis label for the whole page instead of one per panel.
        mtext("Average log2 fold change", outer = TRUE, side = 1, line = 1.2, cex = 0.95)
        mtext("Ranked by log2FC. Dashed line: log2FC = 0.25.",
          outer = TRUE, side = 1, line = 2.9, cex = 0.8
        )
      },
      finally = dev.off()
    )
  }
}

plot_markers_dotplot <- function(obj, markers, resolution, directory) {
  features <- unique(markers$gene)
  if (length(features) == 0L) {
    return(invisible(NULL))
  }
  Idents(obj) <- obj[[paste0("cca_snn_res.", resolution), drop = TRUE]]
  pages <- split(features, ceiling(seq_along(features) / 50))
  for (page in seq_along(pages)) {
    genes <- pages[[page]]
    title <- paste("CCA resolution", resolution, "| top 5 markers per cluster")
    caption <- "All cells; RNA LogNormalize expression. Exploratory cluster markers, not condition DE."
    dots <- DotPlot(
      obj,
      features = rev(genes),
      assay = "RNA",
      col.min = -2.5,
      col.max = 2.5
    ) +
      scale_colour_gradient2(
        low = "#2166AC",
        mid = "#F7F7F7",
        high = "#B2182B",
        midpoint = 0,
        limits = c(-2.5, 2.5)
      ) +
      guides(colour = guide_colourbar(title = "Scaled average\nexpression")) +
      coord_flip() +
      labs(
        title = title,
        subtitle = "Dot size: detection; colour: scaled average expression",
        caption = caption,
        x = NULL,
        y = "Cluster"
      ) +
      theme(
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 9)
      )
    save_marker_plot(dots, sprintf("dotplot_top5_%02d", page), directory,
      width = max(12, nlevels(Idents(obj)) * 0.35), height = max(6, length(genes) * 0.22 + 2)
    )
  }
  invisible(NULL)
}

plot_top_markers_heatmap <- function(obj, markers, resolution, directory) {
  features <- unique(markers$gene)
  if (length(features) == 0L) {
    return(invisible(NULL))
  }
  Idents(obj) <- obj[[paste0("cca_snn_res.", resolution), drop = TRUE]]
  # Fresh scaling over all cells, without the PCA regression covariates.
  obj <- ScaleData(obj,
    assay = "RNA", features = features,
    vars.to.regress = NULL, verbose = FALSE
  )
  colours <- setNames(marker_colours[seq_along(levels(Idents(obj)))], levels(Idents(obj)))
  pages <- split(features, ceiling(seq_along(features) / 50))
  for (page in seq_along(pages)) {
    genes <- pages[[page]]
    heatmap <- DoHeatmap(obj,
      features = genes, assay = "RNA",
      group.colors = colours, label = FALSE, raster = TRUE, lines.width = 1,
      disp.min = -2.5, disp.max = 2.5
    ) +
      scale_fill_gradient2(
        low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
        midpoint = 0, limits = c(-2.5, 2.5), name = "Scaled\nexpression"
      ) +
      labs(
        title = paste("CCA resolution", resolution, "| top 5 markers per cluster"),
        subtitle = "All cells, grouped by cluster; unregressed gene scaling",
        caption = "All cells; RNA LogNormalize expression. Exploratory cluster markers, not condition DE."
      ) +
      theme(
        axis.text.y = element_text(size = 9),
        legend.position = "right",
        legend.box = "vertical",
        legend.justification = "center"
      )
    save_marker_plot(heatmap, sprintf("heatmap_top5_%02d", page), directory,
      width = 20, height = max(7, length(genes) * 0.22 + 2)
    )
  }
  invisible(NULL)
}

plot_marker_features <- function(obj, features, resolution, directory, label) {
  present <- features[features %in% rownames(obj[["RNA"]])]
  missing <- setdiff(features, present)
  if (length(missing) > 0L) {
    warning(label, ": not in object, skipped: ", paste(missing, collapse = ", "))
  }
  if (length(present) == 0L) {
    return(invisible(NULL))
  }

  cluster_panel <- scplotter::CellDimPlot(
    obj,
    group_by = paste0("cca_snn_res.", resolution),
    reduction = "umap_cca",
    label = TRUE,
    # Without this, UMAP labels are sequential 1..n instead of cluster IDs 0..24.
    label_insitu = TRUE
  )
  panels <- c(
    list(cluster_panel),
    lapply(present, function(gene) {
      scplotter::FeatureStatPlot(
        obj,
        features = gene,
        plot_type = "dim",
        reduction = "umap_cca",
        highlight = TRUE,
        theme = "theme_blank"
      ) +
        # theme_blank strips panel titles; restore the gene name per panel.
        ggplot2::labs(title = gene)
    })
  )
  combined <- patchwork::wrap_plots(panels, ncol = 3) +
    patchwork::plot_annotation(
      title = paste(label, "lineage markers | CCA resolution", resolution),
      caption = "All cells; RNA LogNormalize expression on the CCA UMAP. Exploratory cluster markers, not condition DE."
    )
  cols <- 3
  rows <- ceiling(length(panels) / cols)
  save_marker_plot(combined, paste0("features_", label), directory,
    width = cols * 5.2, height = rows * 4.8
  )
  invisible(NULL)
}

write_marker_detection <- function(obj, features, resolution, directory, label) {
  present <- features[features %in% rownames(obj[["RNA"]])]
  missing <- setdiff(features, present)
  if (length(missing) > 0L) {
    warning(label, ": not in object, skipped: ", paste(missing, collapse = ", "))
  }
  if (length(present) == 0L) return(invisible(NULL))

  expression <- FetchData(
    obj[["RNA"]],
    vars = present,
    cells = rownames(obj@meta.data),
    layer = "data",
    clean = FALSE
  )
  if (!identical(rownames(expression), rownames(obj@meta.data)) ||
      !all(present %in% colnames(expression))) {
    stop(label, ": FetchData did not return the requested cells and genes")
  }
  cells <- obj@meta.data[[paste0("cca_snn_res.", resolution)]]
  detected <- as.matrix(expression) > 0
  if (anyNA(detected) || anyNA(cells)) {
    stop(label, ": missing values in expression or cluster identities")
  }
  pct <- sapply(present, function(gene) tapply(detected[, gene], cells, mean))
  pct <- pct[order(as.numeric(rownames(pct))), , drop = FALSE]
  if (anyNA(pct)) {
    stop(label, ": detection undefined for some cluster-gene combination")
  }
  write.csv(round(pct * 100, 2), file.path(directory, paste0("detection_", label, ".csv")))
  invisible(pct)
}
