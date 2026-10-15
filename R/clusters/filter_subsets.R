# Owner-approved removal of unexpected T/B/plasma marker patterns.
# Run with Slurm from the workspace root: Rscript R/clusters/filter_subsets.R
# Original subset RDS files are read-only. No normalization or PCA is rerun.
options(warn = 1)
suppressPackageStartupMessages(library(Seurat))
source("R/io.R") # save_final_object()

job_id <- Sys.getenv("SLURM_JOB_ID")
if (!nzchar(job_id)) stop("Submit this saved-object filtering run through Slurm.")
run_id <- paste0("job-", job_id)
parent <- "job-44493014"
folders <- c("T_cell_nk", "B_cells", "Myeloid", "Granulocytes", "Endothelial", "Epithelial", "Stromal")
inputs <- file.path("data", parent, "subsets", folders, "subclustered.rds")
input_info <- file.info(inputs)[, c("size", "mtime")]
report_dir <- file.path("results", parent, "subset_filtering", run_id)
data_dirs <- file.path("data", parent, "subsets", folders, run_id)
if (dir.exists(report_dir) || any(dir.exists(data_dirs))) {
  stop("This job's output directory already exists; refusing to overwrite it.")
}
dir.create(report_dir, recursive = TRUE)
info <- file.path(report_dir, ".INFO")
rules <- c(
  "Detection uses existing SoupX-corrected RNA counts >0; unavailable features are not biological zeros.",
  "T pattern: CD3E and CD3D, together with CD5 or CD6.",
  "B pattern: at least 3 of CD79A, CD79B, MS4A1, CD19, PAX5 and BANK1.",
  "Plasma pattern: JCHAIN, MZB1 and TNFRSF17 together.",
  "T patterns are allowed in T_cell_nk. B and plasma patterns are allowed in B_cells.",
  "Remove the union of unexpected patterns. A cell matching several patterns is removed once.",
  "The owner approved marker-based exclusions; flags are not confirmed cell identities or measured contamination rates.",
  "No SoupX, normalization, feature filtering, HVG selection, PCA, integration or clustering is rerun.",
  "Existing HVGs, scale layers, reductions, graphs and cluster labels are inherited, not refitted after removal.",
  "Mixed has no saved subset RDS and is not filtered. Granulocytes_warning_check duplicates Granulocytes and is not processed."
)
writeLines(c(
  paste("Status: running; producer", run_id),
  paste("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("Inputs:", paste(inputs, collapse = "; ")),
  "Parent producer: job-44493014; SoupX count producer: job-44409807.",
  paste("Code: R/clusters/filter_subsets.R; R", getRversion(), "; Seurat", packageVersion("Seurat")),
  "Seed: not applicable; deterministic filtering.", rules,
  paste("Log:", paste0("logs/subset-filtering-", job_id, ".log"))
), info)

mapping <- read.csv("results/misc/Ensembl_ID/gene-mapping/ensembl_to_symbol_full.csv", stringsAsFactors = FALSE)
mapping <- mapping[mapping$species == "sus_scrofa" & mapping$symbol_review_required == "no" &
  mapping$gene_symbol != "" & !mapping$symbol_evidence_incomplete %in% c("true", "yes"), ]
b_genes <- c("CD79A", "CD79B", "MS4A1", "CD19", "PAX5", "BANK1")
genes <- c("CD3D", "CD3E", "CD5", "CD6", b_genes, "JCHAIN", "MZB1", "TNFRSF17")
summary_rows <- sample_rows <- cluster_rows <- mapping_rows <- list()

tryCatch({
  for (i in seq_along(folders)) {
    folder <- folders[i]
    message("Reading original subset: ", inputs[i])
    obj <- readRDS(inputs[i])
    meta <- obj[[]]
    cells <- colnames(obj)
    assay <- obj[["RNA"]]
    features <- rownames(assay)
    resolved <- vapply(genes, function(g) {
      if (g %in% features) return(g)
      ids <- intersect(mapping$ensembl_id[mapping$gene_symbol == g], features)
      if (length(ids) != 1L) stop("Cannot resolve required marker uniquely: ", g)
      ids
    }, character(1))
    counts <- matrix(NA_real_, length(cells), length(genes), dimnames = list(cells, genes))
    for (layer in Layers(assay, search = "^counts($|\\.)")) {
      x <- LayerData(assay, layer = layer)
      use <- which(resolved %in% rownames(x))
      counts[match(colnames(x), cells), use] <- t(as.matrix(x[resolved[use], , drop = FALSE]))
      rm(x)
    }
    if (anyNA(counts)) stop("A required marker is unavailable in a sample layer; refusing to infer absence.")
    t_pattern <- counts[, "CD3E"] > 0 & counts[, "CD3D"] > 0 &
      (counts[, "CD5"] > 0 | counts[, "CD6"] > 0)
    b_detected <- rowSums(counts[, b_genes, drop = FALSE] > 0)
    b_pattern <- b_detected >= 3L
    plasma_pattern <- counts[, "JCHAIN"] > 0 & counts[, "MZB1"] > 0 & counts[, "TNFRSF17"] > 0
    foreign_t <- t_pattern & folder != "T_cell_nk"
    foreign_b <- b_pattern & folder != "B_cells"
    foreign_plasma <- plasma_pattern & folder != "B_cells"
    remove <- foreign_t | foreign_b | foreign_plasma
    parent_col <- if ("parent_harmony_snn_res.0.7" %in% names(meta)) "parent_harmony_snn_res.0.7" else "harmony_snn_res.0.7"
    reasons <- vapply(seq_along(cells), function(j) {
      paste(c("T_pattern", "B_pattern", "plasma_pattern")[c(foreign_t[j], foreign_b[j], foreign_plasma[j])], collapse = ";")
    }, character(1))
    flags <- data.frame(subset = folder, cell_id = cells, sample = meta$sample,
      parent_cluster = as.character(meta[[parent_col]]),
      T_pattern = t_pattern, B_markers_detected = b_detected, B_pattern = b_pattern,
      plasma_pattern = plasma_pattern, foreign_T = foreign_t, foreign_B = foreign_b,
      foreign_plasma = foreign_plasma, removed = remove, reason = reasons,
      nCount_RNA = meta$nCount_RNA, nFeature_RNA = meta$nFeature_RNA,
      counts, check.names = FALSE)
    write.csv(flags, file.path(report_dir, paste0(folder, "_cell_flags.csv")), row.names = FALSE)
    write.csv(flags[remove, ], file.path(report_dir, paste0(folder, "_removed_cells.csv")), row.names = FALSE)
    sample_summary <- aggregate(cbind(before = rep(1L, nrow(flags)), removed = as.integer(remove),
      T_flag = as.integer(foreign_t), B_flag = as.integer(foreign_b), plasma_flag = as.integer(foreign_plasma)),
      by = list(sample = as.character(meta$sample)), FUN = sum)
    sample_summary$subset <- folder
    sample_summary$retained <- sample_summary$before - sample_summary$removed
    sample_rows[[i]] <- sample_summary
    cluster_summary <- aggregate(cbind(before = rep(1L, nrow(flags)), removed = as.integer(remove),
      T_flag = as.integer(foreign_t), B_flag = as.integer(foreign_b), plasma_flag = as.integer(foreign_plasma)),
      by = list(parent_cluster = as.character(meta[[parent_col]])), FUN = sum)
    cluster_summary$subset <- folder
    cluster_summary$retained <- cluster_summary$before - cluster_summary$removed
    cluster_rows[[i]] <- cluster_summary
    mapping_rows[[i]] <- data.frame(subset = folder, gene = genes, feature = unname(resolved))

    obj$foreign_T_flag <- foreign_t
    obj$foreign_B_flag <- foreign_b
    obj$foreign_plasma_flag <- foreign_plasma
    obj$foreign_lineage_remove <- remove
    selected_cells <- cells[!remove]
    obj <- subset(obj, cells = selected_cells)
    obj@misc$foreign_lineage_filter <- list(
      producer = run_id, input = inputs[i], rules = rules,
      removed_cells = cells[remove], flags_csv = file.path(report_dir, paste0(folder, "_cell_flags.csv")),
      inherited_analysis = TRUE, before = length(cells), removed = sum(remove), retained = length(selected_cells))
    dir.create(data_dirs[i], recursive = TRUE)
    path <- file.path(data_dirs[i], "subclustered.rds")
    message(folder, ": remove ", sum(remove), "/", length(cells), "; saving ", path)
    save_final_object(obj, path)
    reopened <- readRDS(path)
    stopifnot(
      identical(colnames(reopened), selected_cells),
      !any(colnames(reopened) %in% cells[remove]),
      identical(reopened[[]], obj[[]]), identical(Idents(reopened), Idents(obj)),
      identical(reopened@assays, obj@assays), identical(reopened@reductions, obj@reductions),
      identical(reopened@graphs, obj@graphs), identical(reopened@neighbors, obj@neighbors),
      identical(reopened@misc$foreign_lineage_filter, obj@misc$foreign_lineage_filter),
      !any(reopened$foreign_lineage_remove))
    summary_rows[[i]] <- data.frame(subset = folder, before = length(cells),
      T_flag = sum(foreign_t), B_flag = sum(foreign_b), plasma_flag = sum(foreign_plasma),
      removed = sum(remove), retained = length(selected_cells),
      pct_removed = 100 * mean(remove), output = path, validation = "passed")
    write.csv(do.call(rbind, summary_rows), file.path(report_dir, "summary.csv"), row.names = FALSE)
    cat(paste0(folder, ": saved and reopened; ", sum(remove), " removed; ", length(selected_cells), " retained.\n"), file = info, append = TRUE)
    message(folder, ": saved RDS reopened and verified")
    rm(obj, reopened, assay, counts, meta, flags); gc()
  }
  summary <- do.call(rbind, summary_rows)
  write.csv(do.call(rbind, sample_rows), file.path(report_dir, "by_sample.csv"), row.names = FALSE)
  write.csv(do.call(rbind, cluster_rows), file.path(report_dir, "by_parent_cluster.csv"), row.names = FALSE)
  write.csv(do.call(rbind, mapping_rows), file.path(report_dir, "marker_features.csv"), row.names = FALSE)
  stopifnot(identical(file.info(inputs)[, c("size", "mtime")], input_info))
  table_rows <- vapply(seq_len(nrow(summary)), function(i) {
    r <- summary[i, ]
    sprintf("| %s | %d | %d | %d | %d | %d | %.3f%% | %d |", r$subset, r$before,
      r$T_flag, r$B_flag, r$plasma_flag, r$removed, r$pct_removed, r$retained)
  }, character(1))
  writeLines(c(
    "# Subset cell filtering report", "",
    "| Subset | Before | T flags | B flags | Plasma flags | Removed once | Removed percentage | Retained |",
    "|---|---:|---:|---:|---:|---:|---:|---:|", table_rows, "",
    sprintf("The run removed %d unique cells from %d input cells and retained %d cells.",
      sum(summary$removed), sum(summary$before), sum(summary$retained)), "",
    "## Approved selection rules", "", rules, "",
    "## Inputs and outputs", "",
    paste0("The filtering producer is `", run_id, "`. Original input RDS files were not overwritten."),
    paste0("- `", inputs, "` becomes `", summary$output, "`."), "",
    "Each subset has a cell_flags.csv file covering every input cell, with marker counts and exclusion reasons.",
    "Each subset also has a removed_cells.csv file containing only excluded cell IDs and their evidence.",
    "The by_sample.csv and by_parent_cluster.csv files report exclusion counts and retained counts.",
    "The marker_features.csv file records the gene names and actual features used, including pig Ensembl ID resolutions.", "",
    "## Validation and limitations", "",
    "Every saved RDS was reopened. Retained cell IDs, metadata, identities, assays/count layers, reductions, graphs and neighbors matched the saved in-memory object.",
    "All original input file sizes and modification times were unchanged. Required markers were available in every sample layer.",
    "Flags are marker-based exclusion criteria, not experimentally confirmed identities or a measured contamination rate.",
    "The inherited analysis and existing PCA figures still describe fits made before removal. Recompute HVGs, scaling and PCA from the new RDS files before choosing downstream PCs.",
    paste0("The original execution log is `logs/subset-filtering-", job_id, ".log`."), ""
  ), file.path(report_dir, "report.md"))
  cat(paste0("Status: completed; all seven saved RDS files reopened and verified; originals unchanged.\nFinished: ",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n"), file = info, append = TRUE)
  message("Subset filtering completed. Report: ", file.path(report_dir, "report.md"))
}, error = function(e) {
  cat(paste0("Status: failed; ", conditionMessage(e), "\n"), file = info, append = TRUE)
  stop(e)
})
