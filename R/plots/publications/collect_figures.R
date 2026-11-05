# Collect the saved figure files for a single new DEG run under its subset job.
# Numerical tables and fitted objects remain in the analysis directory.
collect_deg_publication_figures <- function(input_dir, output_dir = file.path(dirname(input_dir), "publications_figures")) {
  dirs <- list.dirs(input_dir, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl("^cluster_[0-9]+$", basename(dirs))]
  figures <- unlist(lapply(dirs, function(d) list.files(d, pattern = "[.](pdf|png)$", recursive = TRUE, full.names = TRUE)), use.names = FALSE)
  if (!length(figures)) return(invisible(character()))
  root <- file.path(output_dir, basename(input_dir))
  relative <- substring(figures, nchar(input_dir) + 2L)
  targets <- file.path(root, relative)
  if (any(file.exists(targets))) stop("Publication figure target already exists; refusing to overwrite")
  for (i in seq_along(figures)) {
    dir.create(dirname(targets[i]), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(figures[i], targets[i])) stop("Could not copy publication figure: ", figures[i])
  }
  if (!identical(unname(tools::md5sum(figures)), unname(tools::md5sum(targets)))) stop("Publication copies failed checksum validation")
  invisible(targets)
}
