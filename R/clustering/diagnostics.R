sample_clustering_cells <- function(metadata, n_cells, seed) {
  stopifnot("Diagnostic subset cannot exceed the input cell count" = n_cells <= nrow(metadata))
  metadata <- metadata[order(rownames(metadata), method = "radix"), , drop = FALSE]
  groups <- split(rownames(metadata), as.character(metadata$sample))
  expected <- n_cells * lengths(groups) / nrow(metadata)
  allocation <- floor(expected)
  remaining <- n_cells - sum(allocation)
  if (remaining > 0L) {
    extra <- order(expected - allocation, decreasing = TRUE)[seq_len(remaining)]
    allocation[extra] <- allocation[extra] + 1L
  }
  stopifnot("Diagnostic subset must represent every sample" = all(allocation > 0))
  set.seed(seed)
  cells <- unlist(lapply(seq_along(groups), function(i) {
    sample(groups[[i]], allocation[i])
  }), use.names = FALSE)
  data.frame(cell = cells, sample = as.character(metadata[cells, "sample"]))
}

evaluate_clustering <- function(obj, diagnostic_cells) {
  config <- obj@misc$clustering
  metadata <- obj[[]]
  cells <- diagnostic_cells$cell
  embedding <- Seurat::Embeddings(obj, config$reduction)[
    cells, clustering_dims(config$settings, config$method), drop = FALSE]
  message("Silhouettes: ", config$reduction, " on ", length(cells), " cells")
  # One packed distance vector per method, reused across all resolutions.
  distances <- stats::dist(embedding, method = config$settings$silhouette_distance)
  scores <- sizes <- composition <- summaries <- list()
  for (i in seq_len(nrow(config$partitions))) {
    resolution <- config$partitions$resolution[i]
    column <- config$partitions$column[i]
    labels <- metadata[, column]
    subset_labels <- as.character(metadata[cells, column])
    groups <- factor(subset_labels)
    valid <- nlevels(groups) > 1L && nlevels(groups) < length(cells)
    widths <- rep(NA_real_, length(cells))
    if (valid) {
      widths <- cluster::silhouette(as.integer(groups), distances)[, "sil_width"]
      stopifnot("Invalid silhouette values" =
        all(is.finite(widths) & widths >= -1 & widths <= 1))
    }
    scores[[i]] <- data.frame(cell = cells, resolution = resolution,
      cluster = subset_labels, silhouette = widths)
    counts <- table(labels)
    sampled <- as.integer(table(factor(subset_labels, levels = names(counts))))
    cluster_means <- tapply(widths, factor(subset_labels, levels = names(counts)), mean)
    sizes[[i]] <- data.frame(resolution = resolution, cluster = names(counts),
      n_cells = as.integer(counts), n_diagnostic = sampled,
      mean_silhouette = as.numeric(cluster_means),
      diagnostic_status = ifelse(sampled == 0, "absent",
        ifelse(sampled == 1, "singleton", "represented")))
    composition[[i]] <- as.data.frame(table(cluster = labels,
      sample = as.character(metadata$sample)), stringsAsFactors = FALSE)
    names(composition[[i]])[3] <- "n_cells"
    composition[[i]]$resolution <- resolution
    summaries[[i]] <- data.frame(resolution = resolution, n_clusters = length(counts),
      min_cluster_size = min(counts), n_small_clusters = sum(counts < 10),
      absent_clusters = sum(sampled == 0), singleton_clusters = sum(sampled == 1),
      mean_silhouette = if (valid) mean(widths) else NA_real_,
      mean_cluster_silhouette = if (valid) mean(cluster_means, na.rm = TRUE) else NA_real_,
      silhouette_status = if (valid) "available" else "requires 2 to n-1 sampled clusters")
    gc()
  }
  # Singleton widths retain the standard zero convention; no cells are excluded.
  list(diagnostic_cells = diagnostic_cells, silhouettes = do.call(rbind, scores),
    cluster_sizes = do.call(rbind, sizes), sample_composition = do.call(rbind, composition),
    summary = do.call(rbind, summaries))
}

compare_clustering <- function(runs) {
  # Use every cell for agreement; cluster numbers themselves need not match.
  reference_cells <- sort(rownames(runs[[1]]$assignments), method = "radix")
  partitions <- do.call(rbind, lapply(runs, function(run) {
    data.frame(method = run$method, run$partitions)
  }))
  agreements <- overlaps <- list()
  pairs <- utils::combn(seq_len(nrow(partitions)), 2, simplify = FALSE)
  for (pair in pairs) {
    a <- partitions[pair[1], ]
    b <- partitions[pair[2], ]
    same_method <- a$method == b$method
    adjacent <- same_method && pair[2] == pair[1] + 1L
    if (!adjacent && (same_method || a$resolution != b$resolution)) next
    x <- runs[[a$method]]$assignments[reference_cells, a$column]
    y <- runs[[b$method]]$assignments[reference_cells, b$column]
    key <- data.frame(method_a = a$method, resolution_a = a$resolution,
      method_b = b$method, resolution_b = b$resolution)
    agreements[[length(agreements) + 1L]] <- data.frame(key,
      n_clusters_a = length(unique(x)), n_clusters_b = length(unique(y)),
      adjusted_rand = mclust::adjustedRandIndex(x, y))
    tab <- as.data.frame(table(cluster_a = x, cluster_b = y))
    overlaps[[length(overlaps) + 1L]] <- cbind(key, tab[tab$Freq > 0, ])
  }
  list(agreement = do.call(rbind, agreements), overlap = do.call(rbind, overlaps))
}
