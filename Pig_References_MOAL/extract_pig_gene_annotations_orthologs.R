
install.packages(
  "moalannotgene",
  repos = c(
    "https://fdumbioinfo.r-universe.dev",
    "https://cloud.r-project.org"
  )
)

pig <- getExportedValue("moalannotgene", "genedbss")
orthologs <- getExportedValue("moalannotgene", "orthogenedb")

class(pig)
class(orthologs)

# Check number of rows and columns
dim(pig)
dim(orthologs)

# Display column names
colnames(pig)
colnames(orthologs)

# Preview the first 6 rows
head(pig)
head(orthologs)

# Check memory size
object.size(pig)
object.size(orthologs)

# Create output directory
dir.create("Pig_References", showWarnings = FALSE)

# Use character IDs for reliable matching
pig$GeneID <- as.character(pig$GeneID)

id_cols <- c("tax_id", "GeneID",
             "Other_tax_id", "Other_GeneID")

orthologs[id_cols] <- lapply(
  orthologs[id_cols], as.character
)

# FILE 1: Pig gene annotations
write.table(
  pig,
  "Pig_References/Sus_scrofa_gene_annotations.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)

# FILE 2: All ortholog relationships involving pigs
pig_orthologs <- subset(
  orthologs,
  tax_id == "9823" | Other_tax_id == "9823"
)

write.table(
  unique(pig_orthologs),
  "Pig_References/Sus_scrofa_all_orthologs.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)

# FILE 3: Only human-pig orthologs
human_pig <- subset(
  orthologs,
  (tax_id == "9823" & Other_tax_id == "9606") |
    (tax_id == "9606" & Other_tax_id == "9823")
)

mapping <- data.frame(
  Pig_GeneID = ifelse(
    human_pig$tax_id == "9823",
    human_pig$GeneID,
    human_pig$Other_GeneID
  ),
  Human_GeneID = ifelse(
    human_pig$tax_id == "9606",
    human_pig$GeneID,
    human_pig$Other_GeneID
  ),
  Relationship = human_pig$relationship
)

mapping <- unique(mapping)

# Add pig gene symbols
mapping$Pig_Symbol <- pig$Symbol[
  match(mapping$Pig_GeneID, pig$GeneID)
]

# Add human gene symbols
human <- getExportedValue("moalannotgene", "genedbhs")
human$GeneID <- as.character(human$GeneID)

mapping$Human_Symbol <- human$Symbol[
  match(mapping$Human_GeneID, human$GeneID)
]

write.table(
  mapping,
  "Pig_References/Human_Pig_ortholog_mapping.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)

# Verify results
cat("Pig genes:", nrow(pig), "\n")
cat("All pig ortholog rows:", nrow(pig_orthologs), "\n")
cat("Unique human-pig mappings:", nrow(mapping), "\n")
cat("Mappings with pig symbols:",
    sum(!is.na(mapping$Pig_Symbol)), "\n")

# Display example mappings
head(mapping, 10)

# Show saved files and their location
list.files("Pig_References")
normalizePath("Pig_References")

