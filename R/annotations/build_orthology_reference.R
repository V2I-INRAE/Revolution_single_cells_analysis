# One-time snapshot of Ensembl pig -> human orthology (BioMart archive host, release in
# R/annotations/config.R), used to map pig features onto human annotation-model genes.
# Run from the project root: Rscript R/annotations/build_orthology_reference.R
# Rerunning replaces the snapshot; annotation outputs then need a rerun.
source("R/annotations/config.R")
cfg <- annotation_config$orthology

attributes <- c(
  "ensembl_gene_id",
  "external_gene_name",
  "hsapiens_homolog_ensembl_gene",
  "hsapiens_homolog_associated_gene_name",
  "hsapiens_homolog_orthology_type",
  "hsapiens_homolog_orthology_confidence"
)
# completionStamp appends "[success]"; BioMart can otherwise truncate silently.
query <- paste0(
  '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>',
  '<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" ',
  'completionStamp="1" datasetConfigVersion="0.6">',
  '<Dataset name="sscrofa_gene_ensembl" interface="default">',
  paste0('<Attribute name="', attributes, '"/>', collapse = ""),
  "</Dataset></Query>"
)
response <- httr::POST(
  cfg$biomart,
  body = list(query = query),
  encode = "form",
  httr::timeout(900)
)
httr::stop_for_status(response)
lines <- strsplit(
  httr::content(response, as = "text", encoding = "UTF-8"),
  "\n",
  fixed = TRUE
)[[1]]
if (!length(lines) || lines[[length(lines)]] != "[success]") {
  stop(
    "Incomplete BioMart response: ",
    substr(paste(lines, collapse = "\n"), 1, 300)
  )
}

orthology <- read.delim(
  text = paste(lines[-length(lines)], collapse = "\n"),
  header = FALSE,
  quote = "",
  colClasses = "character",
  na.strings = character(),
  col.names = c(
    "pig_gene_id",
    "pig_gene_name",
    "human_gene_id",
    "human_gene_name",
    "orthology_type",
    "orthology_confidence"
  )
)
if (
  !all(grepl("^ENSSSCG[0-9]+$", orthology$pig_gene_id)) ||
    anyDuplicated(orthology)
) {
  stop("Unexpected BioMart rows")
}
orthology$ensembl_release <- cfg$ensembl_release
orthology$retrieved <- format(Sys.Date())

dir.create(dirname(cfg$file), recursive = TRUE, showWarnings = FALSE)
write.table(orthology, cfg$file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(
  "Ensembl release",
  cfg$ensembl_release,
  "| pig genes",
  length(unique(orthology$pig_gene_id)),
  "| with human orthologue",
  length(unique(orthology$pig_gene_id[nzchar(orthology$human_gene_id)])),
  "| rows",
  nrow(orthology),
  "\n"
)
