# One result tree per saved subset resolution; source objects remain read-only.
subset_resolution_dir <- function(input, resolution) {
  parts <- strsplit(normalizePath(input, mustWork = TRUE), .Platform$file.sep)[[1L]]
  i <- match("subsets", parts)
  if (is.na(i) || length(parts) < i + 3L || !grepl("^job-[0-9]+$", parts[i - 1L]) ||
      !grepl("^job-[0-9]+$", parts[i + 2L])) {
    stop("Cannot derive subset publication directory from input: ", input)
  }
  if (length(resolution) != 1L || is.na(resolution) || !grepl("^[0-9]+([.][0-9]+)?$", resolution)) {
    stop("Supply the saved clustering resolution explicitly")
  }
  file.path("results", parts[i - 1L], "clusters_plots", "subset", parts[i + 1L],
    parts[i + 2L], paste0("res", resolution))
}

subset_publication_dir <- function(input, resolution) {
  file.path(subset_resolution_dir(input, resolution), "publication_figures")
}
