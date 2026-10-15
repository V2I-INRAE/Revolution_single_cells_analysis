# Command-line parsing and option validation for subclustering.

validate_subcluster_options <- function(only_pca, n_hvg, dims) {
  if (only_pca && (!is.null(n_hvg) || !is.null(dims))) {
    stop("--only-pca cannot be combined with --hvg or --dims.")
  }
  counts <- list(hvg = n_hvg, dims = dims)
  for (name in names(counts)) {
    value <- counts[[name]]
    if (!is.null(value) &&
        (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
          value < 1 || value > .Machine$integer.max || value %% 1 != 0)) {
      stop("--", name, " must be a positive integer.")
    }
  }
  invisible(TRUE)
}

parse_subcluster_args <- function(args) {
  usage <- paste(
    "Usage: Rscript R/clusters/subcluster.R <input.rds> \"<subset label>\" <folder>",
    "[--only-pca | --hvg=<number> --dims=<number>]"
  )
  flags <- args[startsWith(args, "--")]
  positional <- args[!startsWith(args, "--")]
  if (length(positional) != 3L ||
      any(!grepl("^(--only-pca|--(hvg|dims)=[0-9]+)$", flags)) ||
      anyDuplicated(sub("=.*$", "", flags))) {
    stop(usage)
  }
  value <- function(name) {
    flag <- grep(paste0("^--", name, "="), flags, value = TRUE)
    if (length(flag)) as.numeric(sub("^[^=]+=", "", flag)) else NULL
  }
  options <- list(only_pca = "--only-pca" %in% flags,
    n_hvg = value("hvg"), dims = value("dims"))
  do.call(validate_subcluster_options, options)
  c(list(input = positional[[1L]], subset_label = positional[[2L]],
    folder = positional[[3L]]), options)
}
