plot_soupx_markers <- function(matrices, clusters, sample_id, output,
  genes = c("CD3D", "CD3E", "GZMA", "PRF1", "MS4A1", "CD79B", "JCHAIN", "MZB1",
    "CD68", "C1QA", "CSF3R", "TREM1", "PECAM1", "CLDN5", "VWF", "COL1A1",
    "DCN", "LUM", "FOXJ1", "RSPH1", "SCGB1A1", "SCGB3A2", "SFTPC", "SFTPB",
    "AGER", "POU2F3"), labels = names(matrices)) {
  stopifnot(identical(names(clusters), colnames(matrices$Raw)))
  missing <- setdiff(genes, rownames(matrices$Raw))
  genes <- intersect(genes, rownames(matrices$Raw))
  writeLines(missing, file.path(output, "missing_panel_genes.txt"))
  if (!length(genes)) stop("No candidate marker genes available for plotting")

  # Only plotting copies are normalized; saved corrected counts remain integers.
  obj <- Seurat::CreateSeuratObject(matrices$Raw, min.cells = 0, min.features = 0)
  obj$preliminary_cluster <- factor(clusters,
    levels = as.character(sort(unique(as.integer(clusters)))))
  obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize",
    scale.factor = 10000)
  for (method in setdiff(names(matrices), "Raw")) {
    obj[[method]] <- SeuratObject::CreateAssay5Object(counts = matrices[[method]])
    obj <- Seurat::NormalizeData(obj, assay = method,
      normalization.method = "LogNormalize", scale.factor = 10000)
  }

  summaries <- list()
  for (method in names(matrices)) {
    assay <- if (method == "Raw") "RNA" else method
    values <- SeuratObject::GetAssayData(obj, assay = assay, layer = "data")[genes, , drop = FALSE]
    for (cluster in levels(obj$preliminary_cluster)) {
      cells <- obj$preliminary_cluster == cluster
      summaries[[length(summaries) + 1L]] <- data.frame(method = method,
        cluster = cluster, gene = genes, cells = sum(cells),
        molecules = Matrix::rowSums(matrices[[method]][genes, cells, drop = FALSE]),
        mean_log1p_normalized_expression = Matrix::rowMeans(values[, cells, drop = FALSE]),
        expressing_fraction = Matrix::rowMeans(values[, cells, drop = FALSE] > 0))
    }
  }
  summaries <- do.call(rbind, summaries)
  write.csv(summaries, file.path(output, "marker_expression_by_cluster.csv"), row.names = FALSE)
  upper <- max(summaries$mean_log1p_normalized_expression)
  for (index in seq_along(matrices)) {
    method <- names(matrices)[index]
    assay <- if (method == "Raw") "RNA" else method
    plot <- scplotter::FeatureStatPlot(obj, features = genes, ident = "preliminary_cluster",
      assay = assay, layer = "data", plot_type = "dot", center_zero = FALSE,
      cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = TRUE,
      show_column_names = TRUE, palette = "magma", lower_cutoff = 0,
      upper_cutoff = upper, name = "Mean log1p normalized expression",
      # Explicit default avoids partial matching dot_size to dot_size_name upstream.
      dot_size = function(x) sum(x > 0, na.rm = TRUE) / length(x),
      dot_size_name = "Fraction expressing", rows_name = "",
      title = paste(sample_id, labels[index], "— fixed pre-correction clusters"),
      row_names_side = "left")
    ggplot2::ggsave(file.path(output, paste0("markers_", method, ".png")),
      plot, width = 13, height = 11, dpi = 300, bg = "white")
  }
  list(genes = genes, missing_genes = missing, cluster_levels = levels(obj$preliminary_cluster),
    paths = file.path(output, paste0("markers_", names(matrices), ".png")),
    layer = "data", normalization = "LogNormalize, scale factor 10000; plotting copies only",
    color_limits = c(0, upper), dot_size = "fraction expressing (0–1)",
    note = "Candidate marker panel only; no cell-type assignments or automatic biological acceptance.")
}
