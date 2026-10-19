# Command-line parsing and option validation for subclustering.

validate_subcluster_options <- function(only_pca, n_hvg, dims,
  k_param = NULL, resolutions = NULL) {
  if (only_pca && (!is.null(n_hvg) || !is.null(dims) ||
      !is.null(k_param) || !is.null(resolutions))) {
    stop("--only-pca cannot be combined with --hvg, --dims, --k-param or --resolutions.")
  }
  counts <- list(hvg = n_hvg, dims = dims, "k-param" = k_param)
  for (name in names(counts)) {
    value <- counts[[name]]
    if (!is.null(value) &&
        (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
          value < 1 || value > .Machine$integer.max || value %% 1 != 0)) {
      stop("--", name, " must be a positive integer.")
    }
  }
  if (!is.null(resolutions) &&
      (!is.numeric(resolutions) || !length(resolutions) ||
        any(!is.finite(resolutions) | resolutions <= 0) || anyDuplicated(resolutions))) {
    stop("--resolutions must contain unique positive numbers separated by commas.")
  }
  invisible(TRUE)
}

parse_subcluster_args <- function(args) {
  usage <- paste(
    "Usage: Rscript R/clusters/subcluster.R <input.rds> \"<subset label>\" <folder>",
    "[--only-pca | --hvg=<number> --dims=<number> --k-param=<number> --resolutions=0.1,0.2,...]"
  )
  flags <- args[startsWith(args, "--")]
  positional <- args[!startsWith(args, "--")]
  if (length(positional) != 3L ||
      any(!grepl(paste0("^(--only-pca|--(hvg|dims|k-param)=[0-9]+|",
        "--resolutions=[0-9]+(\\.[0-9]+)?(,[0-9]+(\\.[0-9]+)?)*)$"), flags)) ||
      anyDuplicated(sub("=.*$", "", flags))) {
    stop(usage)
  }
  value <- function(name) {
    flag <- grep(paste0("^--", name, "="), flags, value = TRUE)
    if (length(flag)) as.numeric(sub("^[^=]+=", "", flag)) else NULL
  }
  resolution_flag <- grep("^--resolutions=", flags, value = TRUE)
  resolutions <- if (length(resolution_flag)) {
    as.numeric(strsplit(sub("^--resolutions=", "", resolution_flag), ",", fixed = TRUE)[[1L]])
  } else NULL
  options <- list(only_pca = "--only-pca" %in% flags,
    n_hvg = value("hvg"), dims = value("dims"),
    k_param = value("k-param"), resolutions = resolutions)
  do.call(validate_subcluster_options, options)
  c(list(input = positional[[1L]], subset_label = positional[[2L]],
    folder = positional[[3L]]), options)
}
