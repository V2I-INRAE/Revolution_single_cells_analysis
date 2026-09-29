# Doublet-method runners used by the standalone monitoring workflow.

run_doubletfinder_monitor <- function(
  seurat_obj,
  sample_id,
  params = qc_params$doublets
) {
  n_cells <- ncol(seurat_obj)
  rate <- if (is.null(params$expected_rate)) {
    bd_multiplet_rate(n_cells, params$bd_multiplet_table)
  } else {
    params$expected_rate
  }
  n_exp <- round(rate * n_cells)
  homotypic <- estimate_homotypic(seurat_obj, params)
  n_exp_adj <- round(n_exp * (1 - homotypic))
  n_exp_adj <- max(n_exp_adj, params$min_expected_doublets)

  bcmvn <- find_optimal_pk(seurat_obj, params)
  selected_idx <- which.max(bcmvn$BCmetric)
  selected_pk <- as.numeric(as.character(bcmvn$pK[selected_idx]))
  stopifnot(!is.na(selected_pk), selected_pk > 0, selected_pk <= 1)

  run_at_pk <- function(pk) {
    set.seed(params$seed)
    result <- DoubletFinder::doubletFinder(
      seurat_obj,
      PCs = params$pcs,
      pN = params$pN,
      pK = pk,
      nExp = n_exp_adj,
      sct = params$sct
    )
    pann_col <- grep("^pANN_", colnames(result[[]]), value = TRUE)
    class_col <- grep(
      "^DF.classifications_",
      colnames(result[[]]),
      value = TRUE
    )
    stopifnot(length(pann_col) == 1, length(class_col) == 1)

    data.frame(
      cell = colnames(result),
      score = result[[pann_col]][[1]],
      class = as.character(result[[class_col]][[1]]),
      row.names = NULL
    )
  }

  cells <- run_at_pk(selected_pk)
  at_boundary <- selected_idx %in% c(1L, nrow(bcmvn))
  sensitivity <- NULL
  if (at_boundary && nrow(bcmvn) > 1L) {
    adjacent_idx <- if (selected_idx == 1L) 2L else nrow(bcmvn) - 1L
    adjacent_pk <- as.numeric(as.character(bcmvn$pK[adjacent_idx]))
    adjacent_cells <- run_at_pk(adjacent_pk)
    adjacent_idx_by_cell <- match(cells$cell, adjacent_cells$cell)
    stopifnot(!anyNA(adjacent_idx_by_cell))
    adjacent_class <- adjacent_cells$class[adjacent_idx_by_cell]
    selected_doublets <- cells$cell[cells$class == "Doublet"]
    adjacent_doublets <- cells$cell[adjacent_class == "Doublet"]
    union_calls <- union(selected_doublets, adjacent_doublets)
    overlap <- if (length(union_calls) == 0L) {
      NA_real_
    } else {
      length(intersect(selected_doublets, adjacent_doublets)) /
        length(union_calls)
    }
    sensitivity <- data.frame(
      sample_id = sample_id,
      selected_pK = selected_pk,
      adjacent_pK = adjacent_pk,
      selected_calls = length(selected_doublets),
      adjacent_calls = length(adjacent_doublets),
      changed_calls = sum(cells$class != adjacent_class),
      jaccard = overlap,
      row.names = NULL
    )
  }

  bcmvn$sample_id <- sample_id
  bcmvn$pK_numeric <- as.numeric(as.character(bcmvn$pK))
  bcmvn$selected <- seq_len(nrow(bcmvn)) == selected_idx
  bcmvn$selected_at_boundary <- bcmvn$selected & at_boundary

  list(
    cells = cells,
    bcmvn = bcmvn,
    sensitivity = sensitivity,
    settings = data.frame(
      sample_id = sample_id,
      n_cells = n_cells,
      expected_rate = rate,
      n_exp = n_exp,
      homotypic = homotypic,
      n_exp_adjusted = n_exp_adj,
      selected_pK = selected_pk,
      pK_at_boundary = at_boundary,
      called_doublets = sum(cells$class == "Doublet"),
      row.names = NULL
    )
  )
}

run_scrublet_monitor <- function(
  seurat_obj,
  sample_id,
  current_threshold = 0.15
) {
  python_home <- Sys.getenv(
    "RETICULATE_PYTHON",
    unset = file.path(getwd(), ".venv-scrublet", "bin", "python")
  )
  scrublet_obj <- scrubletR::get_init_scrublet(
    seurat_obj = seurat_obj,
    python_home = python_home,
    expected_doublet_rate = bd_multiplet_rate(ncol(seurat_obj)),
    min_counts = 3L,
    n_prin_comps = 30L
  )

  observed_scores <- as.numeric(scrublet_obj$doublet_scores_obs_)
  simulated_scores <- as.numeric(scrublet_obj$doublet_scores_sim_)
  automatic_calls <- scrublet_obj$predicted_doublets_
  if (is.null(automatic_calls)) {
    automatic_calls <- rep(NA, length(observed_scores))
    automatic_threshold <- NA_real_
    automatic_available <- FALSE
    warning(
      "Scrublet automatic threshold unavailable for ", sample_id,
      "; retaining scores and current score > 0.15 calls"
    )
  } else {
    automatic_calls <- as.logical(automatic_calls)
    automatic_threshold <- as.numeric(scrublet_obj$threshold_)
    if (length(automatic_threshold) != 1L || !is.finite(automatic_threshold)) {
      stop("Scrublet returned invalid automatic threshold for ", sample_id)
    }
    automatic_available <- TRUE
  }
  stopifnot(length(observed_scores) == ncol(seurat_obj))

  current_calls <- observed_scores > current_threshold
  list(
    cells = data.frame(
      cell = colnames(seurat_obj),
      score = observed_scores,
      automatic_class = ifelse(
        is.na(automatic_calls),
        NA_character_,
        ifelse(automatic_calls, "Doublet", "Singlet")
      ),
      current_class = ifelse(current_calls, "Doublet", "Singlet"),
      row.names = NULL
    ),
    scores = data.frame(
      sample_id = sample_id,
      population = c(
        rep("Observed cells", length(observed_scores)),
        rep("Simulated doublets", length(simulated_scores))
      ),
      score = c(observed_scores, simulated_scores),
      row.names = NULL
    ),
    settings = data.frame(
      sample_id = sample_id,
      n_cells = ncol(seurat_obj),
      expected_rate = bd_multiplet_rate(ncol(seurat_obj)),
      automatic_threshold = automatic_threshold,
      current_threshold = current_threshold,
      automatic_calls = if (automatic_available) {
        sum(automatic_calls)
      } else {
        NA_integer_
      },
      automatic_available = automatic_available,
      current_calls = sum(current_calls),
      simulated_fraction_above_automatic = if (is.na(automatic_threshold)) {
        NA_real_
      } else {
        mean(simulated_scores > automatic_threshold)
      },
      simulated_fraction_above_current = mean(
        simulated_scores > current_threshold
      ),
      row.names = NULL
    )
  )
}
