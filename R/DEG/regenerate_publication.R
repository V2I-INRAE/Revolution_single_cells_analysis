# Rscript R/DEG/regenerate_publication.R <subset index 1:5> [--diagnostics-only] [--clusters=1,2]
# Owner-approved replacement of existing PDF/PNG figures; scientific files read-only.
suppressPackageStartupMessages(library(ggplot2))
source("R/DEG/render_saved_pseudobulk.R")
subsets <- c(Endothelial = "44831031", Epithelial = "44831032",
  Granulocytes = "44831033", Macrophage_Monocytes = "44831034", DC = "44831038")
index <- as.integer(commandArgs(trailingOnly = TRUE)[1])
stopifnot(length(index) == 1L, !is.na(index), index %in% seq_along(subsets))
diagnostics_only <- "--diagnostics-only" %in% commandArgs(trailingOnly = TRUE)
subset <- names(subsets)[index]
root <- file.path("results/job-44542861/clusters_plots/subset", subset,
  paste0("job-", subsets[index]), "DEG")
record <- readLines(file.path(root, ".INFO"))
selection <- grep("^Volcano/heatmap selection:", record, value = TRUE)
stopifnot(length(selection) == 1L)
# Resolve display thresholds from this frozen run record, not current pipeline defaults.
fdr_cutoff <- as.numeric(sub(".*FDR<=([0-9.]+).*", "\\1", selection))
logfc_cutoff <- as.numeric(sub(".*abs\\(log2FC\\)>=([0-9.]+).*", "\\1", selection))
stopifnot(is.finite(fdr_cutoff), is.finite(logfc_cutoff))
directories <- dirname(Sys.glob(file.path(root, "cluster_*", "pseudobulk_fit.rds")))
expected_clusters <- c(5L, 5L, 1L, 5L, 3L)[index]
stopifnot(length(directories) == expected_clusters)
cluster_arg <- grep("^--clusters=", commandArgs(trailingOnly = TRUE), value = TRUE)
if (length(cluster_arg)) {
  requested <- as.integer(strsplit(sub("^--clusters=", "", cluster_arg), ",", fixed = TRUE)[[1]])
  available <- as.integer(sub("cluster_", "", basename(directories)))
  stopifnot(!anyNA(requested), all(requested %in% available))
  directories <- directories[available %in% requested]
}
for (directory in directories) {
  label <- paste(subset, "cluster", as.integer(sub("cluster_", "", basename(directory))))
  cluster_record <- readLines(file.path(directory, ".INFO"))
  seed_line <- grep("^Seed:", cluster_record, value = TRUE)
  stopifnot(length(seed_line) == 1L)
  seed <- as.numeric(sub("^Seed: ([0-9]+).*", "\\1", seed_line))
  stopifnot(is.finite(seed))
  info_path <- file.path(directory, ".INFO")
  note <- c("", "## Publication-size figure replacement (2026-10-22)",
    "Status: running", paste("Started:", Sys.time()),
    paste("Plotting allocation:", Sys.getenv("SLURM_JOB_ID"), "array index:", index),
    "Script: R/DEG/regenerate_publication.R; saved-results renderer and approved publication_style.R.",
    paste("Read-only input:", directory), paste("Recorded display selection:", selection),
    paste("Recorded label seed:", seed, "; GSEA parameters read from saved gsea_results.rds."),
    if (diagnostics_only) "Diagnostics-only layout correction: inner titles, axis wrapping and label repulsion; no other figures regenerated." else
      "Replaces only 12 PDFs and 12 matching PNGs in diagnostics/, plots/, gsea/hallmark/.",
    "Approved cluster06 style: 170mm wide; labels9pt, axis10pt, titles11pt, tags14pt; MDS dots3.2mm, axes/guides0.25mm; embedded-font PDF and300dpi PNG.",
    "MDS display labels omit shared REVO prefix, explicitly noted in caption; full data IDs unchanged.",
    "No fitting, DE/GSEA tests, source-object writes, changed thresholds/colours, sample or gene exclusions.",
    "Run record: logs/deg-publication-rollout.INFO",
    paste("Evidence:", sprintf("logs/deg-publication-rollout-%s_%d.log",
      Sys.getenv("SLURM_ARRAY_JOB_ID", unset = Sys.getenv("SLURM_JOB_ID", unset = "local")), index)))
  writeLines(c(cluster_record, note), info_path)
  files <- list.files(directory, pattern = "\\.(csv|rds)$", recursive = TRUE, full.names = TRUE)
  files <- files[!grepl("/publication_preview/", files, fixed = TRUE)]
  before <- tools::md5sum(files)
  pdfs <- list.files(directory, pattern = "\\.pdf$", recursive = TRUE, full.names = TRUE)
  pdfs <- pdfs[!grepl("/publication_preview/", pdfs, fixed = TRUE)]
  stopifnot(length(pdfs) == 12L)
  if (diagnostics_only) pdfs <- pdfs[grepl("/diagnostics/", pdfs, fixed = TRUE)]
  started <- Sys.time()
  tryCatch({
    counts <- render_saved_pseudobulk(directory, directory, label, fdr_cutoff, logfc_cutoff, seed, diagnostics_only)
    stopifnot(identical(before, tools::md5sum(files)),
      all(file.info(pdfs)$mtime >= started),
      all(file.info(sub("\\.pdf$", ".png", pdfs))$mtime >= started))
    note[note == "Status: running"] <- "Status: figures regenerated; PDF reopening and visual review pending."
    writeLines(c(cluster_record, note, paste("Finished:", Sys.time()),
      paste("Input consistency PASS:", paste(names(counts), counts, collapse = "; ")),
      paste("Scientific files unchanged:", length(files), "CSV/RDS files byte-identical by MD5.")), info_path)
    cat("PASS", directory, ":", length(pdfs), "PDF/PNG pairs;", length(files), "scientific files unchanged\n")
  }, error = function(e) {
    writeLines(c(cluster_record, note, paste("Status: failed:", conditionMessage(e)),
      paste("Scientific files unchanged:", identical(before, tools::md5sum(files)))), info_path)
    stop(e)
  })
}
cat("COMPLETED", subset, length(directories), "clusters;", length(directories) * if (diagnostics_only) 2L else 12L, "PDF/PNG pairs\n")
