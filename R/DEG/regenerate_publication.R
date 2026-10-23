# Rscript R/DEG/regenerate_publication.R <subset index 1:5> [--diagnostics-only] [--clusters=1,2]
# Or: Rscript R/DEG/regenerate_publication.R --root=<saved analysis root> [same options]
# Owner-approved replacement of production PDF/PNG figures; scientific files read-only.
suppressPackageStartupMessages(library(ggplot2))
source("R/DEG/render_saved_pseudobulk.R")
subsets <- c(Endothelial = "44831031", Epithelial = "44831032",
  Granulocytes = "44831033", Macrophage_Monocytes = "44831034", DC = "44831038")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(all(grepl("^(--root=.+|--clusters=[0-9,]+|--diagnostics-only|[1-5])$", args)))
positionals <- args[!startsWith(args, "--")]
root_arg <- grep("^--root=", args, value = TRUE)
cluster_arg <- grep("^--clusters=", args, value = TRUE)
stopifnot(length(positionals) <= 1L, length(root_arg) <= 1L, length(cluster_arg) <= 1L)
diagnostics_only <- "--diagnostics-only" %in% args
if (length(root_arg)) {
  root <- sub("^--root=", "", root_arg)
} else {
  stopifnot(length(positionals) == 1L)
  index <- as.integer(positionals)
  root <- file.path("results/job-44542861/clusters_plots/subset", names(subsets)[index],
    paste0("job-", subsets[index]), "DEG")
}
root <- normalizePath(path.expand(root), mustWork = TRUE)
record <- readLines(file.path(root, ".INFO"))
directories <- sort(dirname(Sys.glob(file.path(root, "cluster_*", "pseudobulk_fit.rds"))))
if (!length(directories)) stop("No saved pseudobulk fits in ", root, "; no figures can be regenerated.")
if (length(cluster_arg)) {
  requested <- as.integer(strsplit(sub("^--clusters=", "", cluster_arg), ",", fixed = TRUE)[[1]])
  available <- as.integer(sub("^cluster_", "", basename(directories)))
  stopifnot(length(requested) > 0L, !anyNA(requested), !anyDuplicated(requested))
  if (!all(requested %in% available)) {
    stop("Requested clusters have no eligible saved fit: ", paste(setdiff(requested, available), collapse = ", "))
  }
  directories <- directories[available %in% requested]
}
total_pairs <- 0L
for (directory in directories) {
  cluster_record <- readLines(file.path(directory, ".INFO"))
  stopifnot(startsWith(cluster_record[1], "# "))
  label <- sub("^# +", "", cluster_record[1])
  seed_line <- grep("^Seed:", cluster_record, value = TRUE)
  stopifnot(length(seed_line) == 1L, grepl("^Seed: [0-9]+", seed_line))
  seed <- as.numeric(sub("^Seed: ([0-9]+).*", "\\1", seed_line))
  stopifnot(is.finite(seed))
  # Resolve display settings from frozen records, never mutable pipeline defaults.
  selection <- grep("^Volcano/heatmap selection:", cluster_record, value = TRUE)
  if (!length(selection)) selection <- grep("^Volcano/heatmap selection:", record, value = TRUE)
  fdr_cutoff <- logfc_cutoff <- NA_real_
  if (!diagnostics_only) {
    stopifnot(length(selection) == 1L,
      grepl("FDR<=[0-9.]+", selection), grepl("abs\\(log2FC\\)>=[0-9.]+", selection))
    fdr_cutoff <- as.numeric(sub(".*FDR<=([0-9.]+).*", "\\1", selection))
    logfc_cutoff <- as.numeric(sub(".*abs\\(log2FC\\)>=([0-9.]+).*", "\\1", selection))
    stopifnot(is.finite(fdr_cutoff), is.finite(logfc_cutoff))
  }
  pb <- readRDS(file.path(directory, "pseudobulk_fit.rds"))
  hallmark_contrasts <- character()
  if (!diagnostics_only) {
    # Missing/blocked GSEA stops a full request before this cluster's figures/record are replaced.
    gsea <- read_saved_pseudobulk_gsea(directory, names(pb$results))
    hallmark_contrasts <- names(gsea$results)
  }
  pdfs <- saved_pseudobulk_figure_paths(directory, names(pb$results), hallmark_contrasts, diagnostics_only)
  targets <- c(pdfs, sub("\\.pdf$", ".png", pdfs))
  scope <- if (diagnostics_only) {
    hallmark_state <- if (file.exists(file.path(directory, "gsea", "hallmark", "unresolved_collisions.csv"))) {
      "blocked by identifier collisions"
    } else if (file.exists(file.path(directory, "gsea", "hallmark", "gsea_results.rds"))) {
      "saved results present but not read or verified"
    } else "missing saved results"
    paste("Partial scope: diagnostics only; DE/heatmap/Hallmark figures not regenerated. Hallmark:", hallmark_state)
  } else "Full production figure scope; completed saved Hallmark results required."
  info_path <- file.path(directory, ".INFO")
  note <- c("", paste("## Publication-size figure replacement", Sys.Date()),
    "Status: running", paste("Started:", Sys.time()),
    paste("Plotting allocation:", Sys.getenv("SLURM_JOB_ID", unset = "local")),
    "Script: R/DEG/regenerate_publication.R; saved-results renderer and approved publication_style.R.",
    paste("Read-only input:", directory), paste("Recorded display selection:", paste(selection, collapse = "; ")),
    paste("Recorded label seed:", seed, "; GSEA parameters read only for full regeneration."),
    scope, paste("Production targets:", length(pdfs), "PDFs and matching PNGs; preview directories untouched."),
    "MDS display labels omit shared REVO prefix; full data IDs unchanged.",
    "No fitting, DE/GSEA tests, source-object writes, changed thresholds, sample or gene exclusions.")
  # Protect all scientific files, including any previews; only constructed targets are checked as figures.
  files <- list.files(directory, pattern = "\\.(csv|rds)$", recursive = TRUE, full.names = TRUE)
  before <- tools::md5sum(files)
  writeLines(c(cluster_record, note), info_path)
  started <- Sys.time()
  tryCatch({
    counts <- render_saved_pseudobulk(directory, directory, label, fdr_cutoff, logfc_cutoff, seed, diagnostics_only)
    target_info <- file.info(targets)
    stopifnot(identical(before, tools::md5sum(files)),
      all(!is.na(target_info$mtime)), all(target_info$size > 0), all(target_info$mtime >= started))
    note[note == "Status: running"] <- if (diagnostics_only) {
      "Status: partial regeneration; diagnostics replaced only. PDF reopening and visual review pending."
    } else "Status: full production figures regenerated; PDF reopening and visual review pending."
    writeLines(c(cluster_record, note, paste("Finished:", Sys.time()),
      paste("Input consistency PASS (requested scope):", paste(names(counts), counts, collapse = "; ")),
      paste("Scientific files unchanged:", length(files), "CSV/RDS files byte-identical by MD5.")), info_path)
    cat("PASS", directory, ":", length(pdfs), "PDF/PNG pairs;", length(files), "scientific files unchanged\n", scope, "\n")
  }, error = function(e) {
    writeLines(c(cluster_record, note, paste("Status: failed:", conditionMessage(e)),
      paste("Scientific files unchanged:", identical(before, tools::md5sum(files)))), info_path)
    stop(e)
  })
  total_pairs <- total_pairs + length(pdfs)
}
cat(if (diagnostics_only) "PARTIAL (diagnostics only)" else "COMPLETED (full production figures)",
  root, length(directories), "eligible saved clusters;", total_pairs, "PDF/PNG pairs; visual review pending\n")
