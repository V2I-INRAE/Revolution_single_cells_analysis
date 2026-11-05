# Publication marker dot plot. Uses the existing subset plot's RNA / unscaled
# average / red palette convention, with Seurat's DotPlot (scCustomize is not
# required). One vector PDF; dimensions and text are sized for 5x reduction.
# Usage: Rscript R/plots/publications/plot_top_pub.R <subclustered.rds> <panel.csv> <resolution> [output_dir]
suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
source('R/clusters/markers.R') # feature names used by the marker producer
source('R/cell_cycle.R')       # collision-checked rename_cc_features()
source('R/plots/dotplot_colour.R')

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L || length(args) > 4L) stop('Usage: plot_top_pub.R <subclustered.rds> <panel.csv> <resolution> [output_dir]')
input <- args[[1L]]
panel_path <- args[[2L]]
resolution <- args[[3L]]
cluster_col <- paste0('harmony_snn_res.', resolution)
panel <- read.csv(panel_path, colClasses = 'character', check.names = FALSE)
if (anyNA(panel$gene) || anyDuplicated(panel[c('cluster', 'gene')])) stop('Invalid marker panel')
genes <- unique(panel$gene)
input_info <- file.info(c(input, panel_path))[, c('size', 'mtime')]
obj <- readRDS(input)
metadata <- obj[[]]
if (!cluster_col %in% colnames(metadata) || anyNA(metadata[[cluster_col]]) ||
    !identical(rownames(metadata), colnames(obj))) stop('Cluster metadata or cell alignment invalid')
cluster_ids <- sort(unique(as.character(metadata[[cluster_col]])), method = 'radix')
if (!setequal(cluster_ids, panel$cluster)) stop('Panel and object clusters differ')
DefaultAssay(obj) <- 'RNA'
obj <- JoinLayers(obj, assay = 'RNA', layers = '^data($|[.])', new = 'data')
# Apply only recorded marker-name translations required by this panel, in memory.
missing <- setdiff(genes, rownames(obj[['RNA']]))
rename_map <- marker_feature_renames[marker_feature_renames %in% missing]
if (length(rename_map)) obj <- rename_cc_features(obj, rename_map)
missing <- setdiff(genes, rownames(obj[['RNA']]))
if (length(missing)) stop('Panel features absent from RNA: ', paste(missing, collapse = ', '))
cat('Feature translations:', paste(names(rename_map), rename_map, sep = '=', collapse = ', '), '\n')
obj[[cluster_col]] <- factor(as.character(metadata[[cluster_col]]), levels = cluster_ids)
plot <- Seurat::DotPlot(obj, features = genes, group.by = cluster_col,
  assay = 'RNA', scale = FALSE, dot.scale = 12)
dat <- plot$data
if (nrow(dat) != length(genes) * length(cluster_ids) || anyNA(dat$pct.exp) ||
    anyNA(dat$avg.exp.scaled) || !setequal(as.character(dat$features.plot), genes) ||
    !setequal(as.character(dat$id), cluster_ids)) stop('DotPlot omitted features or clusters')
# Match producer's displayed pct.1 (rounded to 3 decimals) on every panel row.
# This detects incorrect feature translations, cluster columns and input objects.
index <- match(paste(panel$cluster, panel$gene), paste(dat$id, dat$features.plot))
if (anyNA(index)) stop('DotPlot lacks a panel gene/cluster combination')
diff <- abs(dat$pct.exp[index] / 100 - as.numeric(panel[['pct.1']]))
if (any(diff > 0.0011)) {
  bad <- panel[diff > 0.0011, c('cluster','gene','pct.1'), drop = FALSE]
  print(bad)
  stop('Detection fractions disagree with saved marker table; max difference ', max(diff))
}
cat('Cells:', ncol(obj), '; clusters:', length(cluster_ids), '; genes:', length(genes),
    '; max detection difference:', max(diff), '\n')

# Seurat uses log1p(mean(expm1(log-normalized expression))) when scale=FALSE.
plot <- plot + dotplot_colour_scale(plot) +
  labs(colour = 'Average expression\n(log1p mean)',
    caption = 'YlOrRd colour capped at the 75th percentile of plotted values.') +
  scale_size_area(max_size = 12, limits = c(0, 100), breaks = c(25, 50, 75, 100),
    name = 'Cells expressing (%)') +
  labs(x = NULL, y = 'Cluster') +
  guides(colour = guide_colourbar(title = 'Average expression\n(log1p mean)', display = 'rectangles'),
    size = guide_legend(title = 'Cells expressing (%)')) +
  theme_classic(base_size = 45) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 40, colour = 'black'),
    axis.text.y = element_text(size = 45, colour = 'black'),
    axis.title.y = element_text(size = 45),
    axis.line = element_line(linewidth = 1.5),
    axis.ticks = element_line(linewidth = 1.2),
    legend.title = element_text(size = 40), legend.text = element_text(size = 36),
    legend.key.height = grid::unit(0.6, 'in'),
    plot.margin = margin(35, 45, 35, 35))
out_dir <- if (length(args) == 4L) args[[4L]] else file.path(dirname(dirname(panel_path)), 'publications_figures')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stem <- paste0('dotplot_top5_plus_canonical_res', resolution, '_pub_ylorrd_p75')
out <- file.path(out_dir, paste0(stem, '.pdf'))
info_path <- file.path(out_dir, paste0(stem, '.INFO'))
if (file.exists(out) || file.exists(info_path)) stop('Output already exists; refusing to overwrite: ', out)
writeLines(c('# Publication marker dot plot', 'Status: running', paste('Started:', Sys.time()),
  paste('Input object:', input), paste('Marker panel:', panel_path), paste('Cluster column:', cluster_col),
  'RNA unscaled average; YlOrRd, per-plot 75th percentile colour cap; detection fraction controls dot size.',
  'No marker discovery, normalization, or input writes.'), info_path)
width <- max(33.46, length(genes) * 0.58 + 7)
ggsave(out, plot, width = width, height = 18, units = 'in', limitsize = FALSE, bg = 'white')
cat('Page inches:', width, 'x 18; after 5x reduction:', width / 5, 'x 3.6\n')
if (!identical(rawToChar(readBin(out, 'raw', 5L)), '%PDF-')) stop('Output is not a PDF')
if (!identical(input_info, file.info(c(input, panel_path))[, c('size', 'mtime')])) stop('Input changed')
cat('Output:', out, '\n')
info <- readLines(info_path)
info[info == 'Status: running'] <- 'Status: completed'
writeLines(c(info, paste('Finished:', Sys.time()),
  'Validation: PDF header checked; inputs unchanged; detection fractions matched.'), info_path)
