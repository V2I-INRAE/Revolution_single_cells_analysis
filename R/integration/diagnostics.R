suppressPackageStartupMessages({
  library(Seurat)
})

sample_diagnostic_cells <- function(metadata, n_cells = 20000, seed = 1234) {
  stopifnot(!anyNA(metadata$sample), anyDuplicated(rownames(metadata)) == 0L)
  groups <- split(rownames(metadata), as.character(metadata$sample))
  n_cells <- min(n_cells, nrow(metadata))
  # Largest-remainder allocation preserves proportions and the exact total.
  expected <- n_cells * lengths(groups) / nrow(metadata)
  allocation <- floor(expected)
  remaining <- n_cells - sum(allocation)
  if (remaining > 0L) {
    extra <- order(expected - allocation, decreasing = TRUE)[seq_len(remaining)]
    allocation[extra] <- allocation[extra] + 1L
  }
  stopifnot(
    "Diagnostic subset must represent every sample" = all(allocation > 0)
  )
  set.seed(seed)
  cells <- unlist(lapply(seq_along(groups), function(i) {
    sample(groups[[i]], size = allocation[i])
  }), use.names = FALSE)
  metadata[cells, "sample", drop = FALSE]
}

compute_integration_metrics <- function(
  seurat_obj, cells, dims = 1:30, perplexity = 30,
  reductions = c("pca", "harmony", "integrated_cca")
) {
  stopifnot(!anyDuplicated(cells), all(cells %in% colnames(seurat_obj)))
  metadata <- seurat_obj[[]][cells, "sample", drop = FALSE]
  labels <- as.integer(factor(metadata$sample))
  stopifnot(
    !anyNA(labels),
    "LISI needs more cells than three times perplexity" =
      length(cells) > 3 * perplexity,
    "Silhouettes require between two and n-1 groups" =
      length(unique(labels)) > 1L && length(unique(labels)) < length(cells)
  )
  scores <- lapply(reductions, function(reduction) {
    message("Diagnostics: ", reduction, " on ", length(cells), " cells")
    embedding <- Embeddings(seurat_obj, reduction)[cells, dims, drop = FALSE]
    stopifnot(
      identical(rownames(embedding), rownames(metadata)),
      all(is.finite(embedding))
    )
    ilisi <- lisi::compute_lisi(
      embedding, metadata, "sample",
      perplexity = perplexity
    )$sample
    # Never construct distances for the full production object.
    distances <- stats::dist(embedding)
    silhouette <- cluster::silhouette(labels, distances)[, "sil_width"]
    stopifnot(all(is.finite(ilisi)), all(is.finite(silhouette)))
    data.frame(
      cell = cells, sample = metadata$sample, reduction = reduction,
      ilisi = ilisi, batch_silhouette = silhouette
    )
  })
  summary <- do.call(rbind, lapply(scores, function(x) {
    data.frame(
      reduction = x$reduction[1], n_cells = nrow(x),
      mean_ilisi = mean(x$ilisi), median_ilisi = median(x$ilisi),
      ilisi_q25 = unname(quantile(x$ilisi, 0.25)),
      ilisi_q75 = unname(quantile(x$ilisi, 0.75)),
      batch_asw = mean(x$batch_silhouette),
      mean_absolute_batch_silhouette = mean(abs(x$batch_silhouette))
    )
  }))
  # These are raw, global sample-label scores, not cell-type conservation or
  # scIB scores. Greater mixing alone does not demonstrate better integration.
  list(
    scores = do.call(rbind, scores), summary = summary,
    settings = list(
      dims = dims, perplexity = perplexity, reductions = reductions,
      batch_col = "sample", distance = "euclidean",
      population = "shared proportionally stratified subset"
    )
  )
}
