# Publication heatmap of the top-25 cluster markers, all clusters of one subset.
# Same figure as the marker pipeline (plot_cluster_heatmap), with display names.
# Called by R/plots/publications/main.R.
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})
source("R/clusters/markers.R") # feature names used by the marker producer
source("R/cell_cycle.R") # collision-checked rename_cc_features()
source("R/utils/save_plots.R")
source("R/plots/marker_cluster_plots.R")

plot_heatmap_for_pub <- function(obj, top25_path, cluster_col, output_dir) {
  top25 <- read.csv(top25_path, colClasses = "character", check.names = FALSE)
  # Cluster order, then marker rank; a gene shared by clusters appears once (first).
  top25 <- top25[order(as.integer(top25$cluster)), ]
  genes <- unique(top25$gene)
  gene_labels <- setNames(
    top25$display_gene[match(genes, top25$gene)],
    genes
  )
  DefaultAssay(obj) <- "RNA"
  obj <- JoinLayers(obj, assay = "RNA", layers = "^data($|[.])", new = "data")
  # Apply only recorded marker-name translations required by these genes, in memory.
  missing <- setdiff(genes, rownames(obj[["RNA"]]))
  rename_map <- marker_feature_renames[marker_feature_renames %in% missing]
  if (length(rename_map)) {
    obj <- rename_cc_features(obj, rename_map)
  }
  missing <- setdiff(genes, rownames(obj[["RNA"]]))
  if (length(missing)) {
    stop("Top-25 genes absent from RNA: ", paste(missing, collapse = ", "))
  }
  # Stored scale.data covers HVGs only; rescale exactly the plotted genes with the
  # settings recorded in the object (gene-wise, so HVG values are unchanged).
  recorded <- Seurat::Command(obj, "ScaleData.RNA")@params
  obj <- Seurat::ScaleData(
    obj,
    assay = "RNA",
    features = genes,
    vars.to.regress = recorded$vars.to.regress,
    model.use = recorded$model.use,
    use.umi = recorded$use.umi,
    do.scale = recorded$do.scale,
    do.center = recorded$do.center,
    scale.max = recorded$scale.max,
    verbose = FALSE
  )
  ids <- sort(unique(as.integer(as.character(obj[[cluster_col, drop = TRUE]]))))
  obj[[cluster_col]] <- factor(
    as.character(obj[[cluster_col, drop = TRUE]]),
    levels = as.character(ids)
  )
  Idents(obj) <- cluster_col
  plot_cluster_heatmap(
    obj,
    genes,
    cluster_col,
    output_dir,
    seed = deg_config$seed,
    level_label = "Cluster",
    gene_labels = gene_labels
  )
  invisible(file.path(output_dir, "heatmap_top25.pdf"))
}
