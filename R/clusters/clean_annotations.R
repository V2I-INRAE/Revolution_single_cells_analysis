# Run from the project root: Rscript R/clusters/clean_annotations.R
# Load the parent once, select original child barcodes, remove/reassign, save.
# Source objects are read-only; no normalization, doublet detection or refitting.
# Six removals = 4444 cells; eight transfers = 21157; final total = 234183.
# Rename unconventional_T_cells to Unconventional_T_and_ILCs (27327 cells).
# Targets below are compartments, not new fine-cluster identifiers.

suppressPackageStartupMessages(library(Seurat))
source("R/clusters/annotation_helpers.R", local = TRUE)

cleanup_sources_44542861 <- read.table(text = "folder|label|job|resolution
B_cells|B cells|job-44614229|0.4
CD4_T_cells|CD4_T_cells|job-44614236|0.4
CD8_T_cells|CD8_T_cells|job-44614237|0.4
DC|DC|job-44614230|0.3
Endothelial|Endothelial|job-44614231|0.3
Epithelial|Epithelial|job-44614232|0.4
Granulocytes|Granulocytes|job-44614233|0.3
Macrophage_Monocytes|Macrophage_Monocytes|job-44614234|0.4
NK_cells|NK_cells|job-44614238|0.4
Stromal|Stromal|job-44614235|0.4
unconventional_T_cells|unconventional_T_cells|job-44614239|0.4
T_cells_unresolved|T_cells_unresolved|job-44614240|0.4",
  header = TRUE, sep = "|", colClasses = "character", stringsAsFactors = FALSE)

cleanup_rules_44542861 <- read.table(text = "folder|cluster|action|cells|target
Endothelial|6|remove|318|REMOVE
Epithelial|7|remove|603|REMOVE
Macrophage_Monocytes|9|remove|445|REMOVE
Macrophage_Monocytes|10|remove|245|REMOVE
Stromal|6|remove|319|REMOVE
T_cells_unresolved|2|remove|2514|REMOVE
CD4_T_cells|3|reassign|3748|CD8_T_cells
CD4_T_cells|4|reassign|3035|CD8_T_cells
T_cells_unresolved|1|reassign|4017|CD8_T_cells
T_cells_unresolved|3|reassign|2213|CD8_T_cells
T_cells_unresolved|4|reassign|1669|CD4_T_cells
T_cells_unresolved|5|reassign|1522|Unconventional_T_and_ILCs
NK_cells|3|reassign|2892|Unconventional_T_and_ILCs
NK_cells|4|reassign|2061|Unconventional_T_and_ILCs",
  header = TRUE, sep = "|", colClasses = "character", stringsAsFactors = FALSE)

run_cleanup_44542861 <- function(
  input = "data/job-44542861/annotated_clustered.rds",
  output = "data/job-44542861/clean_annotated_object.rds"
) {
  info_file <- paste0(tools::file_path_sans_ext(output), ".INFO")
  stopifnot(!file.exists(output), !file.exists(info_file))

  # Record the exact rules before loading any objects.
  info <- c(
    "# Annotation cleanup",
    "Status: running",
    paste("Started:", Sys.time()),
    paste("Producer job:", Sys.getenv("SLURM_JOB_ID", unset = "manual")),
    paste("Input:", input),
    paste("Output:", output),
    "Method: approved doublet removal and compartment reassignment by exact barcode.",
    "Normalized RNA retained; no refitting; no random seed needed (deterministic).",
    "Sources:",
    capture.output(print(cleanup_sources_44542861, row.names = FALSE)),
    "Rules:",
    capture.output(print(cleanup_rules_44542861, row.names = FALSE))
  )
  completed <- FALSE
  validation <- "Validation: pending"
  writeLines(info, info_file)
  on.exit({
    info[info == "Status: running"] <- paste("Status:",
      if (completed) "completed" else "failed/incomplete; see execution log")
    writeLines(c(info, paste("Finished:", Sys.time()), validation), info_file)
  }, add = TRUE)

  # Preserve the original compartment and child-cluster provenance.
  obj <- readRDS(input)
  obj$cleanup_original_subset <- as.character(obj$subset)
  obj$cleanup_source_cluster <- NA_character_
  obj$cleanup_source_job <- NA_character_
  obj$cleanup_source_resolution <- NA_character_

  # Discard old fits once; cell removal/reassignment does not refit them.
  for (name in c(Reductions(obj), Graphs(obj), Neighbors(obj))) {
    obj[[name]] <- NULL
  }
  rna <- obj[["RNA"]]
  for (layer in Layers(rna, search = "^scale\\.data($|\\.)")) {
    LayerData(rna, layer = layer) <- NULL
  }
  obj[["RNA"]] <- rna
  rm(rna)

  # Read each original child once; select cells by its job and resolution.
  for (i in seq_len(nrow(cleanup_sources_44542861))) {
    entry <- cleanup_sources_44542861[i, ]
    child_path <- file.path(dirname(input), "subsets", entry$folder,
      entry$job, "subclustered.rds")
    child <- readRDS(child_path)
    cluster_col <- paste0("harmony_snn_res.", entry$resolution)
    clusters <- as.character(child[[cluster_col, drop = TRUE]])
    at <- match(colnames(child), colnames(obj))
    obj$cleanup_source_cluster[at] <- clusters
    obj$cleanup_source_job[at] <- entry$job
    obj$cleanup_source_resolution[at] <- entry$resolution

    rules <- cleanup_rules_44542861[cleanup_rules_44542861$folder == entry$folder, ]
    for (j in seq_len(nrow(rules))) {
      rule <- rules[j, ]
      barcodes <- colnames(child)[clusters == rule$cluster]
      message(entry$folder, " C", rule$cluster, ": ", rule$action,
        " ", length(barcodes), " cells")
      stopifnot(length(barcodes) == as.integer(rule$cells))

      if (rule$action == "remove") {
        obj <- remove_doublets(obj, barcodes)
      } else {
        obj <- reassign_cells(obj, barcodes, entry$label, rule$target)
      }
    }
    rm(child)
    gc(verbose = FALSE)
  }

  # Rename the original unconventional pool; incoming cells already use this label.
  barcodes <- colnames(obj)[obj$cleanup_original_subset == "unconventional_T_cells"]
  obj <- reassign_cells(obj, barcodes,
    "unconventional_T_cells", "Unconventional_T_and_ILCs")

  # Keep qualified lineage notes separate from compartment routing.
  lineage <- c(
    "Granulocytes/4" = "mast_cells",
    "unconventional_T_cells/8" = "ILC2_like",
    "T_cells_unresolved/5" = "gamma_delta_like_T",
    "NK_cells/3" = "NK_like_T_subtype_unresolved",
    "NK_cells/4" = "NK_like_T_subtype_unresolved"
  )
  key <- paste(obj$cleanup_original_subset, obj$cleanup_source_cluster, sep = "/")
  notes <- unname(lineage[key])
  notes[is.na(notes)] <- "compartment_only"
  obj$cleanup_lineage_note <- notes
  DefaultAssay(obj) <- "RNA"
  obj@misc$cleanup_44542861 <- list(
    input = input,
    sources = cleanup_sources_44542861,
    rules = cleanup_rules_44542861,
    parent_analysis = obj@misc[c("integration", "clustering")]
  )
  obj@misc$integration <- NULL
  obj@misc$clustering <- NULL

  # Save once, release the in-memory object, then verify the saved output.
  cells <- colnames(obj)
  metadata <- obj[[]]
  counts <- table(obj$subset)
  layers <- Layers(obj[["RNA"]])
  saveRDS(obj, output)
  rm(obj)
  gc(verbose = FALSE)

  saved <- readRDS(output)
  stopifnot(
    ncol(saved) == 234183L,
    identical(colnames(saved), cells),
    identical(saved[[]], metadata),
    identical(table(saved$subset), counts),
    identical(Layers(saved[["RNA"]]), layers),
    length(Reductions(saved)) == 0L,
    length(Graphs(saved)) == 0L,
    length(Neighbors(saved)) == 0L,
    length(Layers(saved[["RNA"]], search = "^scale\\.data($|\\.)")) == 0L,
    identical(as.character(Idents(saved)), as.character(saved$subset))
  )
  validation <- c(
    "Validation: saved cell IDs, metadata, compartment counts, RNA layer names and identities verified; old fits absent.",
    capture.output(print(counts))
  )
  completed <- TRUE
  invisible(output)
}

# Sourcing defines the runner without executing cleanup; the CLI takes no arguments.
if (sys.nframe() == 0L) {
  if (length(commandArgs(trailingOnly = TRUE))) {
    stop("Run without arguments; source this file to set input/output.")
  }
  run_cleanup_44542861()
}
