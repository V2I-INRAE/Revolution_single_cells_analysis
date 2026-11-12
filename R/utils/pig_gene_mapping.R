# Display annotation only: never rename matrix rows or choose GSEA canonical IDs.
prepare_pig_gene_mapping <- function(features,
    mapping_dir = "results/misc/Ensembl_ID/gene-mapping",
    reference_dir = "results/misc/msigdb_pig_2026.1") {
  features <- unique(as.character(features))
  full <- read.csv(file.path(mapping_dir, "ensembl_to_symbol_full.csv"),
    stringsAsFactors = FALSE, check.names = FALSE)
  if (anyDuplicated(full$ensembl_id)) stop("Duplicate Ensembl annotation keys")
  clean <- function(x) { x[is.na(x)] <- ""; trimws(x) }
  for (field in c("gene_symbol", "symbol_source", "symbol_review_required", "symbol_evidence_incomplete")) {
    full[[field]] <- clean(full[[field]])
  }
  # Read pig-side fields only. db_gene_symbol/target_gene_symbol are human names.
  reference <- do.call(rbind, lapply(c("hallmark", "reactome", "gobp"), function(collection) {
    x <- readRDS(file.path(reference_dir, paste0(collection, "_pig.rds")))
    if (any(is.na(x$ortholog_taxon_id) | x$ortholog_taxon_id != 9823L)) stop("Not a pig projection")
    data.frame(ensembl_id = clean(x$ensembl_gene), symbol = clean(x$gene_symbol),
      ncbi = clean(x$ncbi_gene), stringsAsFactors = FALSE)
  }))
  reference <- unique(reference[grepl("^ENSSSCG[0-9]+$", reference$ensembl_id) &
    nzchar(reference$symbol), ])
  reference_symbols <- split(reference$symbol, reference$ensembl_id)
  reference_ncbi <- split(reference$ncbi, reference$ensembl_id)
  # Assess duplicate Ensembl symbols globally, not separately in each top-N table.
  candidates <- unique(c(features, full$ensembl_id, reference$ensembl_id))
  ens <- grepl("^ENSSSCG[0-9]+([.][0-9]+)?$", candidates)
  m <- match(candidates, full$ensembl_id)
  symbol <- ifelse(ens, "", candidates)
  source <- ifelse(ens, "", "existing_feature_name")
  status <- ifelse(ens, "unmapped", "existing_feature_name")
  for (i in which(ens)) {
    row <- m[i]
    if (!is.na(row) && (full$symbol_review_required[row] != "no" ||
        tolower(full$symbol_evidence_incomplete[row]) %in% c("yes", "true"))) {
      status[i] <- "review_required"
      next
    }
    if (!is.na(row) && full$mapping_status[row] == "multiple_symbols_review") {
      status[i] <- "ambiguous_symbol"
      next
    }
    if (!is.na(row) && nzchar(full$gene_symbol[row])) {
      symbol[i] <- full$gene_symbol[row]
      source[i] <- paste0("pig_Ensembl_annotation:", full$symbol_source[row])
      status[i] <- if (grepl("^LOC[0-9]+$", symbol[i])) "locus_only" else "mapped"
      next
    }
    alternatives <- unique(reference_symbols[[candidates[i]]])
    ncbi <- unique(reference_ncbi[[candidates[i]]])
    ncbi <- ncbi[nzchar(ncbi)]
    if (length(alternatives) > 1L || length(ncbi) > 1L) {
      status[i] <- "conflicting_pathway_annotation"
    } else if (length(alternatives) == 1L && !grepl("^LOC[0-9]+$", alternatives)) {
      symbol[i] <- alternatives
      source[i] <- "frozen_MSigDB2026.1_pig_projection"
      status[i] <- "mapped_pathway_rescue"
    }
  }
  duplicated_symbol <- nzchar(symbol) & (duplicated(symbol) | duplicated(symbol, fromLast = TRUE))
  display <- candidates
  named <- ens & nzchar(symbol) & !grepl("^LOC[0-9]+$", symbol)
  # Pathway-only rescues must be collision-free; never resolve ambiguous identities by naming.
  blocked_rescue <- named & duplicated_symbol & status == "mapped_pathway_rescue"
  status[blocked_rescue] <- "pathway_symbol_collision"
  named[blocked_rescue] <- FALSE
  display[named] <- symbol[named]
  qualified <- named & duplicated_symbol
  display[qualified] <- paste0(symbol[qualified], " (", candidates[qualified], ")")
  status[qualified] <- "mapped_symbol_collision"
  out <- data.frame(gene = candidates, gene_symbol = symbol, display_gene = display,
    mapping_source = source, mapping_status = status, stringsAsFactors = FALSE)
  out <- out[match(features, out$gene), , drop = FALSE]
  rownames(out) <- NULL
  if (anyDuplicated(out$gene) || anyDuplicated(out$display_gene)) stop("Display mapping is not unique")
  out
}

annotate_pig_genes <- function(tab, mapping) {
  if (!"gene" %in% names(tab)) stop("Result table lacks its feature key: gene")
  index <- match(as.character(tab$gene), mapping$gene)
  if (anyNA(index)) stop("Features missing from pig annotation: ", paste(unique(tab$gene[is.na(index)]), collapse = ", "))
  for (field in c("gene_symbol", "display_gene", "mapping_source", "mapping_status")) {
    tab[[field]] <- mapping[[field]][index]
  }
  tab
}
