# Cell-cycle scoring, tail labeling and tail removal (normalization stage).

suppressPackageStartupMessages({
  library(Seurat)
})

# --- Score S/G2M sur les data nomralisees. Metadonnees. On utilise une copie
# de l'object avec les layers jointes pour faire le calcul puis on applique
# les metadata / labels a l'obj d'origine.
score_cell_cycle <- function(seurat_obj, cc_config) {
  obj <- rename_cc_features(seurat_obj, cc_config$rename_map)

  message("Scoring cell cycle on a temporary pooled (JoinLayers) copy")

  copy <- JoinLayers(obj)
  s_features <- intersect(Seurat::cc.genes$s.genes, rownames(copy))
  g2m_features <- intersect(Seurat::cc.genes$g2m.genes, rownames(copy))

  message(
    "  Gene sets resolved: S ",
    length(s_features),
    "/43, G2M ",
    length(g2m_features),
    "/54"
  )

  copy <- CellCycleScoring(
    copy,
    s.features = s_features,
    g2m.features = g2m_features,
    ctrl = cc_config$ctrl,
    seed = cc_config$seed,
    verbose = FALSE
  )

  stopifnot(identical(colnames(obj), colnames(copy)))

  obj$S.Score <- copy$S.Score
  obj$G2M.Score <- copy$G2M.Score
  obj$Phase <- copy$Phase

  stopifnot(!anyNA(obj$S.Score), !anyNA(obj$G2M.Score))

  obj <- label_cc_tails(obj, cc_config$tail_cutoffs)

  per_sample <- data.frame(
    sample = names(table(obj$sample)),
    cells = as.integer(table(obj$sample)),
    labeled = as.integer(tapply(obj$cc_tail, obj$sample, sum))
  )
  per_sample$pct <- round(100 * per_sample$labeled / per_sample$cells, 2)
  obj@misc$cell_cycle <- list(
    method = paste(
      "common-reference scoring on a temporary JoinLayers copy;",
      "tail labels applied per sample"
    ),
    ctrl = cc_config$ctrl,
    seed = cc_config$seed,
    n_features = c(S = length(s_features), G2M = length(g2m_features)),
    s_features = s_features,
    g2m_features = g2m_features,
    cutoffs = cc_config$tail_cutoffs,
    per_sample = per_sample
  )
  message(
    "  Cell-cycle labels: ",
    sum(obj$cc_tail),
    "/",
    ncol(obj),
    " cells (",
    sprintf("%.2f%%", 100 * mean(obj$cc_tail)),
    ")"
  )
  obj
}

# Sur 97 gene du cycle cellulaire 90 etaient aussi dans la liste
# cochon. 7 maquants. En recupérant les ENSMBLID j'en ai trouvé 4
# de plus. Mais ils ont un nom different chez le cochon donc mapping ici.
# C'est juste pour permettre a seurat de faire le cell cycle score
rename_cc_features <- function(seurat_obj, rename_map) {
  stopifnot(!anyDuplicated(c(names(rename_map), unname(rename_map))))
  feature_names <- rownames(seurat_obj)
  present <- rename_map[names(rename_map) %in% feature_names]
  collides <- unname(present) %in% setdiff(feature_names, names(present))
  if (any(collides)) {
    stop(
      "Rename target(s) already in the matrix: ",
      paste(unname(present)[collides], collapse = ", ")
    )
  }
  if (length(present) < length(rename_map)) {
    message(
      "  Rename map, absent from panel (no-op): ",
      paste(unname(rename_map)[!names(rename_map) %in% feature_names], collapse = ", ")
    )
  }
  feature_names[match(names(present), feature_names)] <- unname(present)
  rownames(seurat_obj) <- feature_names
  seurat_obj
}

# --- On tag les outliers du cycle cellulaires
label_cc_tails <- function(seurat_obj, cutoffs) {
  cutoff_rows <- match(as.character(seurat_obj$sample), cutoffs$sample)
  in_s <- seurat_obj$S.Score > cutoffs$cut_S[cutoff_rows]
  in_g2m <- seurat_obj$G2M.Score > cutoffs$cut_G2M[cutoff_rows]
  seurat_obj$cc_tail <- in_s | in_g2m
  seurat_obj$cc_tail_class <- factor(
    ifelse(in_s, ifelse(in_g2m, "both", "S"), ifelse(in_g2m, "G2M", "none")),
    levels = c("none", "S", "G2M", "both")
  )
  seurat_obj
}

# --- On supprime les outliers du cycle cellulaires
remove_cc_tail_cells <- function(seurat_obj) {
  n_before <- ncol(seurat_obj)
  removed <- seurat_obj$cc_tail
  samples_before <- unique(as.character(seurat_obj$sample))
  obj <- subset(seurat_obj, cells = colnames(seurat_obj)[!removed])

  stopifnot(ncol(obj) == n_before - sum(removed), !any(obj$cc_tail))
  if (!setequal(samples_before, unique(as.character(obj$sample)))) {
    stop("Cell-cycle tail removal emptied at least one sample.")
  }

  per_sample_removed <- table(sample = as.character(seurat_obj$sample)[removed])
  obj@misc$cell_cycle$removal <- list(
    removed = sum(removed),
    retained = ncol(obj),
    pct_removed = round(100 * sum(removed) / n_before, 2),
    per_sample_removed = per_sample_removed
  )
  message(
    "  Cell-cycle tail removal: ",
    sum(removed),
    "/",
    n_before,
    " cells removed (",
    sprintf("%.2f%%", 100 * sum(removed) / n_before),
    "); ",
    ncol(obj),
    " retained"
  )
  obj
}
