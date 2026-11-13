# Publication marker dot plot. Uses the existing subset plot's RNA / unscaled
# average / red palette convention, with Seurat's DotPlot (scCustomize is not
# required). One vector PDF; dimensions and text are sized for 5x reduction.
# Called by R/plots/publications/main.R.
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})
source('R/clusters/markers.R') # feature names used by the marker producer
source('R/cell_cycle.R') # collision-checked rename_cc_features()
source('R/plots/dotplot_colour.R')
source('R/utils/pig_gene_mapping.R')

plot_top_pub <- function(obj, panel_path, cluster_col, output_dir, cap = 0.9) {
  panel <- read.csv(panel_path, colClasses = 'character', check.names = FALSE)
  genes <- unique(panel$gene)
  metadata <- obj[[]]
  cluster_ids <- sort(
    unique(as.character(metadata[[cluster_col]])),
    method = 'radix'
  )
  # A panel from another resolution would silently plot the wrong clusters.
  if (!setequal(cluster_ids, panel$cluster)) {
    stop('Panel and object clusters differ')
  }
  DefaultAssay(obj) <- 'RNA'
  obj <- JoinLayers(obj, assay = 'RNA', layers = '^data($|[.])', new = 'data')
  # Apply only recorded marker-name translations required by this panel, in memory.
  missing <- setdiff(genes, rownames(obj[['RNA']]))
  rename_map <- marker_feature_renames[marker_feature_renames %in% missing]
  if (length(rename_map)) {
    obj <- rename_cc_features(obj, rename_map)
  }
  missing <- setdiff(genes, rownames(obj[['RNA']]))
  if (length(missing)) {
    stop('Panel features absent from RNA: ', paste(missing, collapse = ', '))
  }
  obj[[cluster_col]] <- factor(
    as.character(metadata[[cluster_col]]),
    levels = cluster_ids
  )
  plot <- Seurat::DotPlot(
    obj,
    features = genes,
    group.by = cluster_col,
    assay = 'RNA',
    scale = FALSE,
    dot.scale = 12
  )
  # Labels are separate from the expression feature keys, including collisions.
  gene_mapping <- prepare_pig_gene_mapping(rownames(obj[['RNA']]))
  axis_labels <- setNames(gene_mapping$display_gene, gene_mapping$gene)
  dat <- plot$data
  # Match producer's displayed pct.1 (rounded to 3 decimals) on every panel row.
  # This detects incorrect feature translations, cluster columns and input objects.
  index <- match(
    paste(panel$cluster, panel$gene),
    paste(dat$id, dat$features.plot)
  )
  diff <- abs(dat$pct.exp[index] / 100 - as.numeric(panel[['pct.1']]))
  if (anyNA(diff) || any(diff > 0.0011)) {
    print(panel[is.na(diff) | diff > 0.0011, c('cluster', 'gene', 'pct.1')])
    stop('Detection fractions disagree with saved marker table')
  }

  # Seurat uses log1p(mean(expm1(log-normalized expression))) when scale=FALSE.
  plot <- plot +
    dotplot_colour_scale(plot, cap = cap) +
    scale_x_discrete(labels = axis_labels) +
    labs(
      colour = 'Average expression\n(log1p mean)',
      caption = sprintf(
        'YlOrRd colour capped at the %gth percentile of plotted values.',
        100 * cap
      )
    ) +
    scale_size_area(
      max_size = 12,
      limits = c(0, 100),
      breaks = c(25, 50, 75, 100),
      name = 'Cells expressing (%)'
    ) +
    labs(x = NULL, y = 'Cluster') +
    guides(
      colour = guide_colourbar(
        title = 'Average expression\n(log1p mean)',
        display = 'rectangles'
      ),
      size = guide_legend(title = 'Cells expressing (%)')
    ) +
    theme_classic(base_size = 45) +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 40,
        colour = 'black'
      ),
      axis.text.y = element_text(size = 45, colour = 'black'),
      axis.title.y = element_text(size = 45),
      axis.line = element_line(linewidth = 1.5),
      axis.ticks = element_line(linewidth = 1.2),
      legend.title = element_text(size = 40),
      legend.text = element_text(size = 36),
      legend.key.height = grid::unit(0.6, 'in'),
      plot.margin = margin(35, 45, 35, 35)
    )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(output_dir, 'dotplot.pdf')
  ggsave(
    out,
    plot,
    width = max(33.46, length(genes) * 0.58 + 7),
    height = 18,
    units = 'in',
    limitsize = FALSE,
    bg = 'white'
  )
  invisible(out)
}
