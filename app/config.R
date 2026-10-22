# Paths are resolved from the project root, not from the caller's home folder.
# Only these views appear in the app or are prepared by default, in this order.
explorer_active_views <- c("parent", "dc", "endothelial", "epithelial",
  "granulocytes", "macrophages", "nk", "stromal")

explorer_config <- function(project_dir) {
  views <- list(
    parent = list(label = "Cleaned parent (resolution 0.7)",
      input = "data/job-44542861/clean_annotated_object.rds",
      embedding_input = "data/job-44542861/annotated_clustered.rds",
      cache = "clean-44542861-original-umap-res0.7", cluster_column = "harmony_snn_res.0.7",
      producer = "44827588 (cleanup), 44611569 (annotation), 44542861 (pipeline)",
      parent = TRUE),
    dc = list(label = "DC (resolution 0.2)",
      input = "data/job-44542861/subsets/DC/job-44831038/subclustered.rds",
      cache = "dc-44831038-res0.2", cluster_column = "harmony_snn_res.0.2",
      producer = "44831038 (subclustering), 44827588 (cleanup), 44542861 (pipeline)",
      parent = FALSE),
    endothelial = list(label = "Endothelial (resolution 0.3)",
      input = "data/job-44542861/subsets/Endothelial/job-44831031/subclustered.rds",
      cache = "endothelial-44831031-res0.3", cluster_column = "harmony_snn_res.0.3",
      producer = "44831031 (subclustering), 44827588 (cleanup), 44542861 (pipeline)",
      parent = FALSE),
    epithelial = list(label = "Epithelial (resolution 0.3)",
      input = "data/job-44542861/subsets/Epithelial/job-44831032/subclustered.rds",
      cache = "epithelial-44831032-res0.3", cluster_column = "harmony_snn_res.0.3",
      producer = "44831032 (subclustering), 44827588 (cleanup), 44542861 (pipeline)",
      parent = FALSE),
    granulocytes = list(label = "Granulocytes (resolution 0.2)",
      input = "data/job-44542861/subsets/Granulocytes/job-44831033/subclustered.rds",
      cache = "granulocytes-44831033-res0.2", cluster_column = "harmony_snn_res.0.2",
      producer = "44831033 (subclustering), 44827588 (cleanup), 44542861 (pipeline)",
      parent = FALSE),
    macrophages = list(label = "Macrophages (resolution 0.4)",
      input = "data/job-44542861/subsets/Macrophage_Monocytes/job-44831034/subclustered.rds",
      cache = "macrophages-44831034-res0.4", cluster_column = "harmony_snn_res.0.4",
      producer = "44831034 (subclustering), 44827588 (cleanup), 44542861 (pipeline)",
      parent = FALSE),
    nk = list(label = "NK cells (resolution 0.3)",
      input = "data/job-44542861/subsets/NK_cells/job-44831035/subclustered.rds",
      cache = "nk-44831035-res0.3", cluster_column = "harmony_snn_res.0.3",
      producer = "44831035 (subclustering), 44827588 (cleanup), 44542861 (pipeline)",
      parent = FALSE),
    stromal = list(label = "Stromal (resolution 0.3)",
      input = "data/job-44542861/subsets/Stromal/job-44831036/subclustered.rds",
      cache = "stromal-44831036-res0.3", cluster_column = "harmony_snn_res.0.3",
      producer = "44831036 (subclustering), 44827588 (cleanup), 44542861 (pipeline)",
      parent = FALSE)
  )
  lapply(views, function(view) {
    view$input <- file.path(project_dir, view$input)
    if (!is.null(view$embedding_input)) {
      view$embedding_input <- file.path(project_dir, view$embedding_input)
    }
    view$cache <- file.path(project_dir, "app/cache", view$cache)
    view$reduction <- "umap_harmony"
    view
  })
}
