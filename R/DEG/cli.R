# Explicit inputs for one subset run; no references are downloaded or inferred
# from mutable historical run directories.
run_deg_cli <- function(args) {
  usage <- paste("Usage: Rscript R/DEG/main.R --input=<subclustered.rds>",
    "--subset=<label> --cluster-col=<column> --clusters=all",
    "--output-root=<subset-job/res0.X/DEG> [--references=<frozen reference directory>]",
    "(default: results/misc/msigdb_pig_2026.1)",
    "[--mapping=<pig mapping CSV>] [--pigs=<pig,pig>] [--min-matched-pigs=3]",
    "[--gsea-exclusions=<CSV with cluster,gene,optional contrast (* = all)>]")
  if (identical(args, "--help")) return(cat(usage, "\n"))
  if (any(!grepl("^--[a-z-]+=.+$", args))) stop(usage)
  keys <- sub("^--([^=]+)=.*$", "\\1", args)
  values <- sub("^--[^=]+=", "", args)
  allowed <- c("input", "subset", "cluster-col", "clusters", "output-root", "references", "mapping", "pigs", "min-matched-pigs", "gsea-exclusions")
  required <- c("input", "subset", "cluster-col", "clusters", "output-root")
  if (anyDuplicated(keys) || any(!keys %in% allowed) || !all(required %in% keys)) stop(usage)
  opt <- as.list(stats::setNames(values, keys))
  if (is.null(opt$references)) opt$references <- "results/misc/msigdb_pig_2026.1"
  refs <- file.path(opt$references, c("hallmark_pig.rds", "reactome_pig.rds", "gobp_pig.rds"))
  mapping <- if (is.null(opt$mapping)) "results/misc/Ensembl_ID/gene-mapping/pig_to_human_orthologues_annotated.csv" else opt$mapping
  if (!all(file.exists(c(opt$input, mapping, refs)))) stop("Input object, mapping or frozen reference missing; no analysis started")
  exclusions <- list()
  if (!is.null(opt$`gsea-exclusions`)) {
    x <- read.csv(opt$`gsea-exclusions`, stringsAsFactors = FALSE, colClasses = "character")
    if (!"contrast" %in% names(x)) x$contrast <- "*"
    if (!all(c("cluster", "gene") %in% names(x)) || anyNA(x[, c("cluster", "gene", "contrast")]) ||
        any(!nzchar(x$cluster)) || any(!nzchar(x$gene)) ||
        any(!x$contrast %in% c("*", pseudobulk_comparisons()$contrast)) ||
        anyDuplicated(x[, c("cluster", "gene", "contrast")])) stop("Invalid GSEA exclusions CSV")
    exclusion_keys <- ifelse(x$contrast == "*", x$cluster, paste(x$cluster, x$contrast, sep = "/"))
    exclusions <- split(x$gene, exclusion_keys)
  }
  minimum <- if (is.null(opt$`min-matched-pigs`)) 3L else as.integer(opt$`min-matched-pigs`)
  pigs <- if (is.null(opt$pigs)) NULL else strsplit(opt$pigs, ",", fixed = TRUE)[[1L]]
  clusters <- strsplit(opt$clusters, ",", fixed = TRUE)[[1L]]
  if (!identical(clusters, "all") && any(!grepl("^[0-9]+$", clusters))) stop("Use all or non-negative integer cluster IDs")
  dir.create("logs", showWarnings = FALSE)
  stem <- file.path("logs", paste0("deg-", format(Sys.time(), "%Y%m%d-%H%M%S"), "-", Sys.getpid()))
  log_path <- paste0(stem, ".log"); info_path <- paste0(stem, ".INFO")
  info <- c("# Pseudobulk DEG and GSEA", "Status: running", paste("Started:", Sys.time()),
    paste("Input:", opt$input), paste("Output:", opt$`output-root`), paste("References:", opt$references),
    paste("Cluster column:", opt$`cluster-col`, "; clusters:", paste(clusters, collapse = ", ")),
    paste("Pigs:", if (is.null(pigs)) "all" else paste(pigs, collapse = ", "), "; minimum matched:", minimum),
    paste("Additional GSEA exclusions:", if (is.null(opt$`gsea-exclusions`)) "none" else opt$`gsea-exclusions`),
    "Uniform GSEA retry: initial1000; whole-contrast retry10000 with seed1234 if any unavailable estimates; stop if still unavailable.",
    paste("Evidence:", log_path), "Method details and per-cluster statuses are recorded in the output .INFO.")
  writeLines(info, info_path)
  cat("Run record:", info_path, "\nLog:", log_path, "\n")
  log <- file(log_path, open = "wt")
  sink(log); sink(log, type = "message")
  finished <- FALSE
  on.exit({
    sink(type = "message"); sink(); close(log)
    info[info == "Status: running"] <- if (finished) "Status: finished; saved-output validation pending" else "Status: failed/incomplete; see log"
    writeLines(c(info, paste("Finished:", Sys.time())), info_path)
  }, add = TRUE)
  run_pseudobulk(input = opt$input, subset_label = opt$subset, cluster_col = opt$`cluster-col`,
    clusters = clusters, output_root = opt$`output-root`, reference_dir = opt$references,
    mapping_file = mapping, gsea_exclusions = exclusions, pigs = pigs, min_matched_pigs = minimum,
    log_path = log_path, collections = c("hallmark", "reactome", "gobp"))
  finished <- TRUE
}
