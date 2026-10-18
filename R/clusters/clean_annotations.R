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

  writeLines(c
    (paste("Input:", input), paste("Output:", output), paste("Started:", Sys.time()),
    paste("Producer job:", Sys.getenv("SLURM_JOB_ID", unset = "manual")),
    "Method: approved doublet removal and compartment reassignment by exact barcode.",
    "Normalized RNA retained; no refitting; no random seed needed (deterministic).",
    "Status: running", capture.output(print(cleanup_sources_44542861, row.names = FALSE)),
    capture.output(print(cleanup_rules_44542861, row.names = FALSE))), info_file
  )

  annotated_cluster <- readRDS(input)
  annotated_cluster$cleanup_original_subset <- as.character(annotated_cluster$subset)
  annotated_cluster$cleanup_source_cluster <- NA_character_
  annotated_cluster$cleanup_source_job <- NA_character_
  annotated_cluster$cleanup_source_resolution <- NA_character_

  # Discard old fits once; they are not refitted by cell removal/reassignment.
  for (name in c(
      Reductions(annotated_cluster),
      Graphs(annotated_cluster),
      Neighbors(annotated_cluster)
    )
  ) {
    annotated_cluster[[name]] <- NULL
  }
  rna <- annotated_cluster[["RNA"]]

  for (layer in Layers(rna, search = "^scale\\.data($|\\.)")) LayerData(rna, layer = layer) <- NULL

  annotated_cluster[["RNA"]] <- rna

  rm(rna)

  for (i in seq_len(nrow(cleanup_sources_44542861))) {
    entry <- cleanup_sources_44542861[i, ]
    child <- readRDS(file.path(dirname(input), "subsets", entry$folder, entry$job, "subclustered.rds"))
    clusters <- as.character(child[[paste0("harmony_snn_res.", entry$resolution), drop = TRUE]])
    at <- match(colnames(child), colnames(annotated_cluster))
    annotated_cluster$cleanup_source_cluster[at] <- clusters
    annotated_cluster$cleanup_source_job[at] <- entry$job
    annotated_cluster$cleanup_source_resolution[at] <- entry$resolution

    rules <- cleanup_rules_44542861[cleanup_rules_44542861$folder == entry$folder, ]

    for (j in seq_len(nrow(rules))) {
      rule <- rules[j, ]
      barcodes <- colnames(child)[clusters == rule$cluster]
      message(entry$folder, " C", rule$cluster, ": ", rule$action, " ", length(barcodes), " cells")
      stopifnot(length(barcodes) == as.integer(rule$cells))

      if (rule$action == "remove") {
        annotated_cluster <- remove_doublets(annotated_cluster, barcodes)
      } else {
        annotated_cluster <- reassign_cells(annotated_cluster, barcodes, entry$label, rule$target)
      }
    }
    rm(child)
    gc(verbose = FALSE)
  }

  barcodes <- colnames(annotated_cluster)[annotated_cluster$cleanup_original_subset == "unconventional_T_cells"]
  annotated_cluster <- reassign_cells(annotated_cluster, barcodes,
    "unconventional_T_cells", "Unconventional_T_and_ILCs")
  # Retain the previously agreed qualified annotations separately from routing.
  lineage <- c("Granulocytes/4" = "mast_cells", "unconventional_T_cells/8" = "ILC2_like",
    "T_cells_unresolved/5" = "gamma_delta_like_T", "NK_cells/3" = "NK_like_T_subtype_unresolved",
    "NK_cells/4" = "NK_like_T_subtype_unresolved")
  key <- paste(annotated_cluster$cleanup_original_subset, annotated_cluster$cleanup_source_cluster, sep = "/")
  notes <- unname(lineage[key])
  notes[is.na(notes)] <- "compartment_only"
  annotated_cluster$cleanup_lineage_note <- notes
  DefaultAssay(annotated_cluster) <- "RNA"
  annotated_cluster@misc$cleanup_44542861 <- list(input = input,
    sources = cleanup_sources_44542861, rules = cleanup_rules_44542861,
    parent_analysis = annotated_cluster@misc[c("integration", "clustering")])
  annotated_cluster@misc$integration <- NULL
  annotated_cluster@misc$clustering <- NULL

  # Save once, then reopen for a small saved-output check.
  cells <- colnames(annotated_cluster)
  counts <- table(annotated_cluster$subset)
  layers <- Layers(annotated_cluster[["RNA"]])
  saveRDS(annotated_cluster, output)
  rm(annotated_cluster)
  gc(verbose = FALSE)
  saved <- readRDS(output)
  stopifnot(ncol(saved) == 234183L, identical(colnames(saved), cells),
    identical(table(saved$subset), counts), identical(Layers(saved[["RNA"]]), layers),
    length(Reductions(saved)) == 0L,
    identical(as.character(Idents(saved)), as.character(saved$subset)))
  cat("\nStatus: completed; saved cell IDs, counts, RNA layers and identities verified.\n",
    capture.output(print(counts)), sep = "\n", file = info_file, append = TRUE)
  invisible(output)
}

# Sourcing is safe. Old --dry-run/CLI arguments are rejected, not silently run.
if (sys.nframe() == 0L) {
  if (length(commandArgs(trailingOnly = TRUE))) stop("Run without arguments; source this file to set input/output.")
  run_cleanup_44542861()
}
