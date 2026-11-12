# One-time snapshot of KEGG pig (ssc) pathways for GSEA; academic use of the KEGG REST API.
# Run from the project root: Rscript R/DEG/build_kegg_reference.R
# Rerunning replaces the snapshot with the current KEGG release; DEG results then need a rerun.
suppressPackageStartupMessages(library(KEGGREST))
source("R/DEG/config.R")

pathways <- keggList("pathway", "ssc")
links <- keggLink("ssc", "pathway")
# keggList() returns only the entry type; the symbol is the part of column 4 before the semicolon.
genes <- read.delim(
  "https://rest.kegg.jp/list/ssc",
  header = FALSE,
  quote = "",
  col.names = c("id", "type", "position", "definition")
)
has_symbol <- grepl(";", genes$definition, fixed = TRUE) # otherwise free text, not a symbol
genes <- setNames(
  ifelse(has_symbol, sub(";.*$", "", genes$definition), NA_character_),
  genes$id
)

kegg <- data.frame(
  gs_id = sub("^path:", "", names(links)),
  ncbi_gene = sub("^ssc:", "", unname(links))
)
# Non-numeric ssc identifiers (mitochondrial genes) have no NCBI Gene ID.
kegg <- kegg[grepl("^[0-9]+$", kegg$ncbi_gene), ]
kegg$gene_symbol <- unname(genes[paste0("ssc:", kegg$ncbi_gene)])
title <- sub(" - Sus scrofa \\(pig\\)$", "", unname(pathways[kegg$gs_id]))
kegg$gs_description <- title
kegg$gs_name <- paste0(
  "KEGG_",
  gsub("^_|_$", "", gsub("[^A-Z0-9]+", "_", toupper(title)))
)
kegg$ensembl_gene <- NA_character_
info <- strsplit(keggInfo("kegg"), "\n")[[1]]
kegg$db_version <- paste(
  "KEGG pathway database",
  sub("^.*\\s", "", info[grepl("^\\s+pathway\\s", info)])
)
kegg$retrieved <- format(Sys.Date())
kegg <- unique(kegg)

dir.create(
  dirname(deg_config$kegg_reference),
  recursive = TRUE,
  showWarnings = FALSE
)
saveRDS(kegg, deg_config$kegg_reference)
cat(
  "KEGG release",
  unique(kegg$db_version),
  "| pathways",
  length(unique(kegg$gs_id)),
  "| genes",
  length(unique(kegg$ncbi_gene)),
  "| rows",
  nrow(kegg),
  "\n"
)
