suppressPackageStartupMessages(library(Seurat))

prepare_marker_umaps <- function(obj, panels, resolution) {
  config <- obj@misc$clustering
  column <- config$partitions$column[config$partitions$resolution == resolution]
  stopifnot("Requested clustering resolution/UMAP must exist" =
    length(column) == 1L && column %in% colnames(obj[[]]) &&
      config$umap %in% Reductions(obj))
  genes <- unique(unlist(panels, use.names = FALSE))
  available <- genes %in% rownames(obj[["RNA"]])
  report <- data.frame(gene = genes, status = ifelse(available, "available", "missing"),
    mean_lognorm = NA_real_, max_lognorm = NA_real_)
  if (any(!available)) warning("Skipping missing markers: ", paste(genes[!available], collapse = ", "))
  stopifnot("No requested markers found in RNA" = any(available))

  cells <- colnames(obj)
  genes <- genes[available]
  obj <- subset(obj, features = genes)
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]], layers = "data", new = "data")
  expression <- as.matrix(LayerData(obj, assay = "RNA", layer = "data")[genes, cells, drop = FALSE])
  stopifnot("Joined normalized RNA must cover every cell with finite nonnegative values" =
    identical(colnames(obj), cells) && all(is.finite(expression)) && all(expression >= 0))
  report$mean_lognorm[available] <- rowMeans(expression)
  report$max_lognorm[available] <- apply(expression, 1, max)
  list(object = obj, report = report, column = column, resolution = resolution)
}

plot_marker_umaps <- function(prepared, panels, output_dir) {
  obj <- prepared$object
  config <- obj@misc$clustering
  report <- prepared$report
  write.csv(report, file.path(output_dir, "marker_availability.csv"), row.names = FALSE)
  title <- paste(config$method, "— resolution", prepared$resolution)
  labels <- levels(obj[[]][, prepared$column])
  stopifnot("More clusters than reference palette colours" = length(labels) <= length(clustering_colours))
  reference <- scplotter::CellDimPlot(obj, reduction = config$umap,
    group_by = prepared$column, label = TRUE, label_repel = TRUE,
    palcolor = setNames(clustering_colours[seq_along(labels)], labels),
    pt_size = 3, pt_alpha = 1, raster = TRUE, seed = config$settings$seed,
    xlab = "UMAP 1", ylab = "UMAP 2",
    theme = "theme_blank", title = paste(title, "— cluster reference"))
  save_clustering_plot(reference, "cluster_reference", output_dir, 9, 6)

  for (compartment in names(panels)) {
    for (cell_type in names(panels[[compartment]])) {
      requested <- panels[[compartment]][[cell_type]]
      genes <- requested[requested %in% report$gene[report$status == "available"]]
      if (!length(genes)) next
      message("Marker UMAPs: ", config$method, " / ", cell_type, " / ", paste(genes, collapse = ", "))
      plot <- scplotter::FeatureStatPlot(obj, features = genes,
        reduction = config$umap, plot_type = "dim", assay = "RNA", layer = "data",
        pos_only = "no", bg_cutoff = -Inf, lower_cutoff = 0,
        upper_cutoff = max(report$max_lognorm[match(genes, report$gene)]),
        lower_quantile = 0, upper_quantile = 1, palette = "magma",
        color_name = "Log-normalized RNA", pt_size = 3, pt_alpha = 1, order = "random",
        seed = config$settings$seed, raster = TRUE, raster_dpi = c(600, 600),
        xlab = "UMAP 1", ylab = "UMAP 2", ncol = length(genes), theme = "theme_blank",
        title = if (length(genes) == 1L) genes else NULL)
      skipped <- setdiff(requested, genes)
      caption <- "Colour: per-cell log-normalized RNA, from zero to this figure's maximum. Fixed-size points; no rescaling or averaging."
      if (length(skipped)) caption <- paste(caption, "Skipped:", paste(skipped, collapse = ", "))
      width <- max(9, 6 * length(genes))
      caption <- paste(strwrap(caption, width = floor(width * 10)), collapse = "\n")
      plot <- patchwork::wrap_plots(plot) +
        patchwork::plot_annotation(title = paste(title, "—", cell_type), caption = caption,
          theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0)))
      stem <- tolower(gsub("[^[:alnum:]]+", "_", paste(compartment, cell_type, sep = "_")))
      save_clustering_plot(plot, stem, output_dir, width, 6)
    }
  }
  invisible(report)
}
