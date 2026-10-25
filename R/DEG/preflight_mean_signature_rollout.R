# Explicit owner-authorized rollout scope; source+call from project root.
# Read-only sequential preflight: does NOT score cells, test contrasts or plot.
source("R/DEG/signature_conditions.R")

mean_signature_rollout_manifest <- function() {
  populations <- c("Endothelial", "Epithelial", "Granulocytes", "Macrophage_Monocytes", "DC", "NK_cells", "Stromal")
  jobs <- c("44831031", "44831032", "44831033", "44831034", "44831038", "44831035", "44831036")
  resolutions <- c("0.3", "0.3", "0.2", "0.2", "0.2", "0.2", "0.3")
  clusters <- list(1:5, 1:5, 1L, c(1:4, 6L), 1:3, 1:3, 1:6)
  do.call(rbind, lapply(seq_along(populations), function(i) {
    data.frame(population = populations[i], job = jobs[i], cluster_col = paste0("harmony_snn_res.", resolutions[i]),
      cluster = as.character(clusters[[i]]),
      input = file.path("data/job-44542861/subsets", populations[i], paste0("job-", jobs[i]), "subclustered.rds"),
      output_root = file.path("results/job-44542861/clusters_plots/subset", populations[i], paste0("job-", jobs[i]),
        "mean_signature_scores", sprintf("cluster_%02d", clusters[[i]])),
      action = ifelse(populations[i] == "Endothelial" & clusters[[i]] == 1L, "preserve_completed", "new_after_approval"),
      stringsAsFactors = FALSE)
  }))
}

preflight_mean_signature_rollout <- function(output_root,
    universe_file = "results/misc/Ensembl_ID/signature-scoring/common_count_universe.rds") {
  if (file.exists(output_root)) stop("Use a NEW preflight directory; preserve previous evidence.")
  if (!requireNamespace("SeuratObject", quietly = TRUE) || !requireNamespace("Matrix", quietly = TRUE)) stop("SeuratObject/Matrix required.")
  manifest <- mean_signature_rollout_manifest()
  stopifnot(nrow(manifest) == 28L, sum(manifest$action == "preserve_completed") == 1L)
  dir.create(output_root, recursive = TRUE)
  write.csv(manifest, file.path(output_root, "manifest.csv"), row.names = FALSE)
  universe <- readRDS(universe_file); stopifnot(length(universe$features) == 15902L, length(universe$samples) == 23L)
  protected <- c(unique(manifest$input), universe_file,
    "results/misc/Ensembl_ID/signature-scoring-v2/pig_mapping.csv",
    "results/misc/Ensembl_ID/signature-scoring-v2/link_adjudication.csv", "R/DEG/signatures.R",
    list.files(manifest$output_root[manifest$action == "preserve_completed"], recursive = TRUE, full.names = TRUE, all.files = TRUE))
  protected <- protected[!file.info(protected)$isdir]
  before <- tools::md5sum(protected)
  if (anyNA(before)) stop("Missing original source/frozen/protected file.")
  saveRDS(before, file.path(output_root, "protected_md5.rds"))
  for (i in seq_len(nrow(manifest))) {
    unit <- manifest[i, ]; dest <- file.path(output_root, paste0(unit$population, "_cluster_", unit$cluster, ".INFO"))
    writeLines(c("# Mean signature rollout preflight ONLY", "Status: WAITING_SOURCE_PREFLIGHT",
      paste("Original source:", unit$input), paste("Selection:", unit$cluster_col, unit$cluster),
      paste("Future target:", unit$output_root), paste("Action:", unit$action),
      "No score/test/plot execution. No DEG sample thresholds or pair deletion applied."), dest)
  }
  source_rows <- list(); cluster_rows <- list(); sample_rows <- list(); pair_rows <- list()
  first <- c("P4", "P10", "P4", "N4"); second <- c("N4", "N10", "Control", "Control")
  groups <- c("Control", "P4", "N4", "P10", "N10")
  for (population in unique(manifest$population)) {
    units <- manifest[manifest$population == population, ]; input <- units$input[1]
    message("PREFLIGHT_SOURCE_START ", population, " ", input)
    x <- readRDS(input); md <- x[[]]; required <- c("sample", "pig", "time_point", "pressure", units$cluster_col[1])
    stopifnot(all(required %in% names(md)), !anyNA(md[, required, drop = FALSE]),
      identical(rownames(md), colnames(x)), !anyDuplicated(rownames(md)),
      all(nzchar(as.character(md[[units$cluster_col[1]]]))), all(units$cluster %in% as.character(md[[units$cluster_col[1]]])))
    md$condition <- signature_five_conditions(md)
    panel <- unique(md[, c("sample", "pig", "time_point", "pressure", "condition"), drop = FALSE])
    panel <- panel[order(as.character(panel$sample), method = "radix"), , drop = FALSE]
    stopifnot(nrow(panel) == 23L, !anyDuplicated(panel$sample), setequal(as.character(panel$sample), universe$samples),
      !anyDuplicated(paste(panel$pig, panel$condition)))
    norm <- x@commands[["NormalizeData.RNA"]]
    stopifnot(!is.null(norm), identical(norm@assay.used, "RNA"),
      identical(norm@params$normalization.method, "LogNormalize"), norm@params$scale.factor == 10000)
    dl <- sort(SeuratObject::Layers(x[["RNA"]], search = "^data[.]"), method = "radix")
    cl <- sort(SeuratObject::Layers(x[["RNA"]], search = "^counts[.]"), method = "radix")
    stopifnot(length(dl) == 23L, length(cl) == 23L,
      setequal(sub("^data[.]", "", dl), universe$samples), setequal(sub("^counts[.]", "", cl), universe$samples))
    dc <- lapply(dl, function(l) SeuratObject::Cells(x[["RNA"]], layer = l))
    common_data <- Reduce(intersect, lapply(dl, function(l) SeuratObject::Features(x[["RNA"]], layer = l)))
    common_count <- Reduce(intersect, lapply(cl, function(l) SeuratObject::Features(x[["RNA"]], layer = l)))
    stopifnot(setequal(common_data, universe$features), setequal(common_count, universe$features),
      !anyDuplicated(unlist(dc)), setequal(unlist(dc), rownames(md)))
    max_norm_error <- 0; audited_cells <- 0L
    for (j in seq_along(dl)) {
      sample <- sub("^data[.]", "", dl[j]); ids <- dc[[j]]
      stopifnot(identical(unique(as.character(md[ids, "sample"])), sample),
        setequal(ids, SeuratObject::Cells(x[["RNA"]], layer = paste0("counts.", sample))))
      # Audit one AUTHORIZED cluster cell per layer, never score/retransform input.
      scoped <- ids[as.character(md[ids, units$cluster_col[1]]) %in% units$cluster]
      if (length(scoped)) {
        counts <- SeuratObject::LayerData(x[["RNA"]], layer = paste0("counts.", sample), cells = scoped[1])
        data <- SeuratObject::LayerData(x[["RNA"]], layer = dl[j], features = rownames(counts), cells = scoped[1])
        data <- data[rownames(counts), scoped[1], drop = FALSE]
        expected <- log1p(as.numeric(counts) / sum(counts) * 10000)
        stopifnot(all(is.finite(data@x)), all(data@x >= 0))
        max_norm_error <- max(max_norm_error, max(abs(expected - as.numeric(data))))
        audited_cells <- audited_cells + 1L
      }
    }
    stopifnot(max_norm_error < 1e-10)
    source_rows[[population]] <- data.frame(population = population, input = input, source_cells = nrow(md),
      source_md5 = unname(before[input]), normalization = "LogNormalize", scale_factor = 10000,
      normalization_timestamp = as.character(norm@time.stamp), data_layers = length(dl), count_layers = length(cl),
      data_common_features = length(common_data), count_common_features = length(common_count),
      normalization_audited_cells = audited_cells, normalization_max_error = max_norm_error, status = "PASS")
    for (i in seq_len(nrow(units))) {
      unit <- units[i, ]; ids <- rownames(md)[as.character(md[[unit$cluster_col]]) == unit$cluster]
      stopifnot(length(ids) > 0L, !anyDuplicated(ids))
      sp <- panel; sp$n_cells <- as.integer(table(factor(as.character(md[ids, "sample"]), levels = as.character(panel$sample))))
      sp$population <- population; sp$cluster <- unit$cluster
      sample_rows[[length(sample_rows) + 1L]] <- sp
      available <- sp[sp$n_cells > 0L, ]; all_complete <- TRUE; pair_counts <- integer(4)
      for (j in seq_along(first)) {
        a <- as.character(available$pig[available$condition == first[j]])
        b <- as.character(available$pig[available$condition == second[j]])
        expected_pigs <- as.character(panel$pig[panel$condition == first[j]])
        matched <- sort(intersect(a, b), method = "radix"); pair_counts[j] <- length(matched)
        complete <- setequal(a, expected_pigs) && setequal(b, expected_pigs); all_complete <- all_complete && complete
        pair_rows[[length(pair_rows) + 1L]] <- data.frame(population = population, cluster = unit$cluster,
          contrast = paste(first[j], second[j], sep = "-"), n_first_pigs = length(a), n_second_pigs = length(b),
          n_matched = length(matched), expected_pairs = length(expected_pigs),
          matched_pigs = paste(matched, collapse = ";"), first_only_pigs = paste(setdiff(a, b), collapse = ";"),
          second_only_pigs = paste(setdiff(b, a), collapse = ";"),
          absent_first_pigs = paste(setdiff(expected_pigs, a), collapse = ";"),
          absent_second_pigs = paste(setdiff(expected_pigs, b), collapse = ";"), full_original_pairing = complete)
      }
      exists <- file.exists(unit$output_root)
      status <- if (unit$action == "preserve_completed") "PRESERVE_COMPLETED" else if (exists) "OUTPUT_COLLISION" else if (!all_complete) "HOLD_INFERENCE_PAIRING_DECISION" else "PREFLIGHT_READY_NOT_EXECUTED"
      tmp <- tempfile(); saveRDS(sort(ids, method = "radix"), tmp); ids_md5 <- unname(tools::md5sum(tmp)); unlink(tmp)
      cr <- data.frame(population = population, cluster = unit$cluster, n_cells = length(ids),
        n_samples = nrow(available), n_pigs = length(unique(available$pig)),
        min_observed_sample_cells = min(available$n_cells), missing_samples = paste(sp$sample[sp$n_cells == 0L], collapse = ";"),
        pairs_P4_N4 = pair_counts[1], pairs_P10_N10 = pair_counts[2], pairs_P4_Control = pair_counts[3], pairs_N4_Control = pair_counts[4],
        all_original_pairs_present = all_complete, output_exists = exists, cell_ids_md5 = ids_md5,
        constant_difference_status = "NOT_ASSESSED_REQUIRES_SCORES", status = status, stringsAsFactors = FALSE)
      for (g in groups) {
        cr[[paste0(g, "_cells")]] <- sum(sp$n_cells[sp$condition == g])
        cr[[paste0(g, "_samples")]] <- sum(sp$n_cells[sp$condition == g] > 0L)
      }
      cluster_rows[[length(cluster_rows) + 1L]] <- cr
      if (unit$action == "preserve_completed") {
        saved <- readRDS(file.path(unit$output_root, "signature_scores.rds"))
        stopifnot(exists, identical(saved$settings$input, input), identical(saved$settings$cluster_col, unit$cluster_col),
          identical(saved$settings$clusters, unit$cluster), setequal(colnames(saved$scores), ids),
          identical(saved$metadata, md[colnames(saved$scores), names(saved$metadata), drop = FALSE]))
        rm(saved)
      }
      dest <- file.path(output_root, paste0(population, "_cluster_", unit$cluster, ".INFO"))
      cat(paste0("\nSource normalization/universe/IDs: PASS\nStatus: ", status,
        "\nCells: ", length(ids), "; observed samples: ", nrow(available), "; min observed sample cells: ", min(available$n_cells),
        "\nAvailable pairs P4-N4/P10-N10/P4-Control/N4-Control: ", paste(pair_counts, collapse = "/"),
        "\nMissing original specimens: ", cr$missing_samples,
        "\nConstant paired differences: not assessed; needs future actual scores. No inference run.\n"), file = dest, append = TRUE)
      message("PREFLIGHT_CLUSTER ", population, " C", unit$cluster, " cells=", length(ids), " samples=", nrow(available), " pairs=", paste(pair_counts, collapse = "/"), " ", status)
    }
    rm(x, md, counts, data); gc()
  }
  result <- list(manifest = manifest, sources = do.call(rbind, source_rows), clusters = do.call(rbind, cluster_rows),
    samples = do.call(rbind, sample_rows), pairs = do.call(rbind, pair_rows))
  for (n in c("sources", "clusters", "samples", "pairs")) write.csv(result[[n]], file.path(output_root, paste0(n, ".csv")), row.names = FALSE)
  saveRDS(result, file.path(output_root, "preflight.rds"))
  stopifnot(nrow(result$clusters) == 28L, nrow(result$pairs) == 112L, nrow(result$samples) == 644L,
    identical(result, readRDS(file.path(output_root, "preflight.rds"))), identical(before, tools::md5sum(protected)))
  message("PREFLIGHT_COMPLETE_28_UNITS_7_SOURCES_NO_SCORING_OR_INFERENCE_PROTECTED_HASHES_PASS")
  invisible(result)
}
