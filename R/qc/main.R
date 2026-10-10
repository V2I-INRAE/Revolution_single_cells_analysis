# Run the same QC pipeline on Rhapsody, SoupX or CellBender counts

args <- commandArgs(trailingOnly = TRUE)
input_source <- if (length(args) == 0L) "rhapsody" else args[[1]]
if (length(args) > 1L || !input_source %in% c("rhapsody", "soupx", "cellbender")) {
  stop("Usage: Rscript R/qc/main.R [rhapsody|soupx|cellbender]")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

source("R/qc/params.R")
source("R/qc/io.R")
source("R/qc/filter.R")
source("R/qc/reduce.R")
source("R/qc/scrublet.R")
source("R/qc/plots.R")

samples <- c(
  "REVO26-C", "REVO26-N4", "REVO26-N10", "REVO26-P4", "REVO26-P10",
  "REVO27-C", "REVO27-N4", "REVO27-N10", "REVO27-P4", "REVO27-P10",
  "REVO29-C", "REVO29-N4", "REVO29-P4", "REVO30-C", "REVO30-N4",
  "REVO30-N10", "REVO30-P10", "REVO30-P4", "REVO31-C", "REVO31-N4",
  "REVO31-N8", "REVO31-P10", "REVO31-P4"
)
input_paths <- setNames(vapply(samples, qc_input_path, character(1),
  input_source = input_source), samples)
if (any(!file.exists(input_paths))) {
  stop("Missing QC inputs:\n", paste(input_paths[!file.exists(input_paths)], collapse = "\n"))
}

job_id <- Sys.getenv("SLURM_JOB_ID")

run_id <- if (nzchar(job_id)) paste0("job-", job_id) else
  paste(format(Sys.time(), "%Y%m%d-%H%M%S"), Sys.getpid(), sep = "-")

labeled_dir <- file.path("data", "qc_labeled_data", run_id)
clean_dir <- file.path("data", "clean_concatenated_data", run_id)
output_dir <- file.path("results", "qc", run_id)

for (directory in c(labeled_dir, clean_dir, output_dir)) {
  dir.create(dirname(directory), recursive = TRUE, showWarnings = FALSE)
  if (!dir.create(directory, showWarnings = FALSE)) {
    stop("Cannot create new run directory: ", directory)
  }
}

info_paths <- file.path(c(labeled_dir, clean_dir, output_dir), ".INFO")
info_lines <- c(
  "# QC run", paste0("- Run: ", run_id),
  paste0("- Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "- Status: running", paste0("- Input source: ", input_source),
  "\n## Inputs", paste0("- ", names(input_paths), ": ", input_paths),
  if (input_source == "soupx") "- Ambient correction producer: 44409807; corrected counts reused, not rerun.",
  "\n## Settings",
  "- QC and Scrublet on all input counts; filter only after labeled checkpoint and plots.",
  "- Fixed QC thresholds shared by all samples; no MAD filtering. Ribosomal percentage is diagnostic only, with no filtering.",
  "- Genes >=min_features; upper thresholds inclusive; complexity >min_log10_genes_per_umi. Any QC flag or Scrublet doublet excludes.",
  "- Temporary LogNormalize/PCA/UMAP for diagnostics; saved counts unchanged.",
  paste0("- Filtering: ", paste(names(qc_params$filtering),
    unlist(qc_params$filtering), sep = "=", collapse = "; ")),
  paste0("- Diagnostic variable genes: ", qc_params$doublets$n_variable_features,
    "; PCs: ", paste(qc_params$doublets$pcs, collapse = ","),
    "; seed: ", qc_params$doublets$seed),
  "- Scrublet: min_counts=3; n_prin_comps=30; score cutoff >0.15; seed=0.",
  paste0("- BD expected doublet rate table (cells:percent): ", paste(
    qc_params$doublets$bd_multiplet_table$cells,
    qc_params$doublets$bd_multiplet_table$rate, sep = ":", collapse = ", ")),
  "- Expected rate: linear interpolation; above 20,000 input cells, linear extrapolation from 19,000:4.5% and 20,000:4.7%, without an upper clamp (user-approved).",
  "- BD table: Instrument User Guide, Doc ID 214062 Rev. 3.0, pp. 80–81; captured cells on retrieved beads. Input cell count remains the proxy; extrapolated rates are not vendor-validated measurements.",
  "- Full applied settings and per-sample thresholds are saved in object@misc$qc.",
  if (input_source == "cellbender")
    "- FPR 0.01 filtered counts; source runs and unresolved warnings: data/raw_data/cellbender/.INFO.",
  "\n## Outputs",
  paste0("- Labeled: ", file.path(labeled_dir, "labeled_concatenated.rds")),
  paste0("- Clean: ", file.path(clean_dir, "clean_concatenated_scrublet.rds")),
  paste0("- Figures: ", output_dir),
  paste0("- Log: logs/*qc-pipeline-", if (nzchar(job_id)) job_id else Sys.getpid(), ".log (batch launcher).")
)
for (info_path in info_paths) writeLines(info_lines, info_path)

message("QC run: ", run_id, "\nInput source: ", input_source, "\nLabeled data: ", labeled_dir,
  "\nClean data: ", clean_dir, "\nResults: ", output_dir)

# Keep all input cells until the labeled checkpoint and plots are saved.
labeled_objs <- list()

for (sample in samples) {
  counts <- read_qc_counts(sample, input_source)
  obj <- build_seurat_obj_from_counts(counts, sample)
  rm(counts)
  obj <- inspect_seurat_qc(obj)

  obj <- label_qc_outliers(obj)

  # --- Prepare the temporary reduction and diagnostic UMAP coordinates
  reduced_obj <- preprocess_qc_umap(obj)
  coordinates <- compute_diagnostic_umap(reduced_obj)
  obj <- AddMetaData(obj, coordinates[colnames(obj), , drop = FALSE])

  # --- Label with Scrublet on all input cells; no removal here.
  obj <- detect_scrublet_doublets(obj)
  labeled_objs[[sample]] <- label_retained_cells(obj)
  rm(obj, reduced_obj, coordinates)
  gc()
}

labeled <- concatenate_qc_objects(labeled_objs)
rm(labeled_objs)
gc()
labeled@misc$qc$run_id <- run_id
labeled@misc$qc$input_source <- input_source
labeled@misc$qc$input_paths <- input_paths
labeled@misc$qc$threshold_reference <- "fixed_shared_across_samples"
labeled@misc$qc$feature_min_inclusive <- TRUE
write_qc_object(labeled, file.path(labeled_dir, "labeled_concatenated.rds"))

# --- These five figures compare all input cells with QC-only retained cells.
plot_qc_comparison(labeled, plot_dir = output_dir)
plot_qc_umaps(labeled, plot_dir = output_dir)

clean <- filter_labeled_cells(labeled)
write_qc_object(clean, file.path(clean_dir, "clean_concatenated_scrublet.rds"))

for (info_path in info_paths) writeLines(c(
  sub("- Status: running", "- Status: completed", info_lines, fixed = TRUE),
  "\n## Outcome",
  paste0("- Finished: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("- Cells: ", ncol(labeled), " input; ", ncol(clean), " retained."),
  "- All samples processed; labeled and clean checkpoints and diagnostic figures saved."
), info_path)
