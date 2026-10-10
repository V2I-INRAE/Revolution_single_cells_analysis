top_cluster_markers <- function(markers, n) {
  markers |>
    dplyr::group_by(cluster) |>
    dplyr::slice_head(n = n) |>
    dplyr::ungroup()
}

find_all_cluster_markers <- function(
  obj,
  resolution,
  only_pos = TRUE,
  min_pct = 0.25,
  logfc_threshold = 0.25,
  test_use = "wilcox"
) {
  Idents(obj) <- obj[[paste0("cca_snn_res.", resolution), drop = TRUE]]
  markers <- FindAllMarkers(
    obj,
    assay = "RNA",
    slot = "data",
    test.use = test_use,
    only.pos = only_pos,
    min.pct = min_pct,
    logfc.threshold = logfc_threshold,
    min.diff.pct = -Inf,
    max.cells.per.ident = Inf,
    return.thresh = Inf,
    verbose = FALSE
  )
  if (nrow(markers) == 0L) {
    stop("No markers returned; check the Seurat warnings in the log.")
  }

  markers |>
    dplyr::mutate(
      cluster = factor(cluster, levels = levels(Idents(obj))),
      pct_difference = pct.1 - pct.2
    ) |>
    dplyr::arrange(
      cluster,
      dplyr::desc(avg_log2FC),
      p_val_adj, gene
    )
}
