# Cell alignment and summaries for doublet monitoring.

collect_doublet_cells <- function(
  seurat_obj,
  sample_id,
  coordinates,
  doubletfinder,
  scrublet
) {
  expected_cells <- colnames(seurat_obj)
  stopifnot(
    anyDuplicated(expected_cells) == 0,
    anyDuplicated(coordinates$cell) == 0,
    anyDuplicated(doubletfinder$cells$cell) == 0,
    anyDuplicated(scrublet$cells$cell) == 0,
    setequal(expected_cells, coordinates$cell),
    setequal(expected_cells, doubletfinder$cells$cell),
    setequal(expected_cells, scrublet$cells$cell)
  )

  coordinate_idx <- match(expected_cells, coordinates$cell)
  df_idx <- match(expected_cells, doubletfinder$cells$cell)
  scrublet_idx <- match(expected_cells, scrublet$cells$cell)
  stopifnot(!anyNA(c(coordinate_idx, df_idx, scrublet_idx)))

  data.frame(
    sample_id = sample_id,
    cell = expected_cells,
    UMAP1 = coordinates$UMAP1[coordinate_idx],
    UMAP2 = coordinates$UMAP2[coordinate_idx],
    nFeature_RNA = seurat_obj$nFeature_RNA,
    nCount_RNA = seurat_obj$nCount_RNA,
    doubletfinder_score = doubletfinder$cells$score[df_idx],
    doubletfinder_class = doubletfinder$cells$class[df_idx],
    scrublet_score = scrublet$cells$score[scrublet_idx],
    scrublet_current_class = scrublet$cells$current_class[scrublet_idx],
    scrublet_automatic_class = scrublet$cells$automatic_class[scrublet_idx],
    row.names = NULL
  )
}

collect_nfeature_summary <- function(cell_data) {
  classification <- rbind(
    data.frame(
      sample_id = cell_data$sample_id,
      method = "DoubletFinder",
      rule = "Selected pK",
      class = cell_data$doubletfinder_class,
      nFeature_RNA = cell_data$nFeature_RNA
    ),
    data.frame(
      sample_id = cell_data$sample_id,
      method = "Scrublet",
      rule = "Current score > 0.15",
      class = cell_data$scrublet_current_class,
      nFeature_RNA = cell_data$nFeature_RNA
    ),
    data.frame(
      sample_id = cell_data$sample_id,
      method = "Scrublet",
      rule = "Automatic threshold",
      class = cell_data$scrublet_automatic_class,
      nFeature_RNA = cell_data$nFeature_RNA
    )
  )
  classification <- classification[!is.na(classification$class), ]
  groups <- split(
    classification,
    interaction(
      classification$sample_id,
      classification$method,
      classification$rule,
      classification$class,
      drop = TRUE
    )
  )
  summaries <- lapply(groups, function(group) {
    quartiles <- quantile(group$nFeature_RNA, c(0.25, 0.5, 0.75))
    data.frame(
      sample_id = group$sample_id[[1]],
      method = group$method[[1]],
      rule = group$rule[[1]],
      class = group$class[[1]],
      n = nrow(group),
      q1 = unname(quartiles[[1]]),
      median = unname(quartiles[[2]]),
      q3 = unname(quartiles[[3]]),
      row.names = NULL
    )
  })
  do.call(rbind, summaries)
}

collect_call_summary <- function(doubletfinder_settings, scrublet_settings) {
  df <- doubletfinder_settings
  scrub <- scrublet_settings
  stopifnot(
    identical(df$sample_id, scrub$sample_id),
    df$n_cells == scrub$n_cells
  )

  data.frame(
    sample_id = df$sample_id,
    n_cells = df$n_cells,
    bd_expected_rate = df$expected_rate,
    doubletfinder_n_exp = df$n_exp,
    doubletfinder_homotypic = df$homotypic,
    doubletfinder_n_exp_adjusted = df$n_exp_adjusted,
    doubletfinder_pK = df$selected_pK,
    doubletfinder_pK_at_boundary = df$pK_at_boundary,
    doubletfinder_calls = df$called_doublets,
    doubletfinder_call_rate = df$called_doublets / df$n_cells,
    scrublet_automatic_threshold = scrub$automatic_threshold,
    scrublet_automatic_available = scrub$automatic_available,
    scrublet_automatic_calls = scrub$automatic_calls,
    scrublet_automatic_call_rate = if (scrub$automatic_available) {
      scrub$automatic_calls / scrub$n_cells
    } else {
      NA_real_
    },
    scrublet_current_threshold = scrub$current_threshold,
    scrublet_current_calls = scrub$current_calls,
    scrublet_current_call_rate = scrub$current_calls / scrub$n_cells,
    simulated_fraction_above_automatic = (
      scrub$simulated_fraction_above_automatic
    ),
    simulated_fraction_above_current = scrub$simulated_fraction_above_current,
    row.names = NULL
  )
}
