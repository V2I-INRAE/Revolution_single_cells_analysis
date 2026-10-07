#!/usr/bin/env Rscript
# Install all R packages required by the scrnaseq-seurat-core-analysis workflow.
# Runs inside the renv project so packages land in the project library on /work,
# then records exact versions in renv.lock (snapshot type = "all").
# Optional package names restrict installation, e.g. Rscript scripts/install_packages.R glmGamPoi.

ncpus <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
options(Ncpus = ncpus, repos = c(CRAN = "https://cloud.r-project.org"))
Sys.setenv(MAKEFLAGS = paste0("-j", ncpus))
message("Installing with Ncpus = ", ncpus)

# name = R package name used for post-install verification
pkgs <- list(
  cran = c(
    renv = "renv",
    remotes = "remotes",
    SoupX = "SoupX",
    Seurat = "Seurat",
    SeuratObject = "SeuratObject",
    scCustomize = "scCustomize",
    ggplot2 = "ggplot2",
    ggprism = "ggprism",
    patchwork = "patchwork",
    dplyr = "dplyr",
    harmony = "harmony",
    ashr = "ashr"
  ),
  bioc = c(
    glmGamPoi = "glmGamPoi",
    SingleR = "SingleR",
    celldex = "celldex",
    SingleCellExperiment = "SingleCellExperiment",
    muscat = "muscat",
    DESeq2 = "DESeq2",
    ComplexHeatmap = "ComplexHeatmap"
  ),
  github = c(
    scplotter = "pwwang/scplotter",
    DoubletFinder = "chris-mcginnis-ucsf/DoubletFinder",
    scrubletR = "Moonerss/scrubletR@1c08e58c7a551406819263603161d49e3055effc",
    lisi = "immunogenomics/lisi",
    SeuratData = "satijalab/seurat-data",
    Azimuth = "satijalab/azimuth"
  )
)

requested <- commandArgs(trailingOnly = TRUE)
if (length(requested) > 0) {
  unknown <- setdiff(requested, unlist(lapply(pkgs, names)))
  if (length(unknown) > 0) stop("Unknown packages: ", paste(unknown, collapse = ", "))
  pkgs <- lapply(pkgs, function(x) x[names(x) %in% requested])
}

status <- c()

for (repo in names(pkgs$cran)) {
  message("\n=== [CRAN] ", repo, " ===")
  try(renv::install(repo), silent = TRUE)
  status[repo] <- requireNamespace(repo, quietly = TRUE)
}

for (pkg in names(pkgs$bioc)) {
  message("\n=== [Bioconductor] ", pkg, " ===")
  try(renv::install(paste0("bioc::", pkg)), silent = TRUE)
  status[pkg] <- requireNamespace(pkg, quietly = TRUE)
}

for (pkg in names(pkgs$github)) {
  message("\n=== [GitHub] ", pkgs$github[[pkg]], " ===")
  try(
    remotes::install_github(
      pkgs$github[[pkg]],
      upgrade = "never",      # never rebuild the already-installed dependency stack
      # scplotter's optional Giotto stack is not needed for Seurat QC figures.
      dependencies = if (pkg == "scplotter") NA else TRUE,
      quiet = TRUE
    ),
    silent = TRUE
  )
  status[pkg] <- requireNamespace(pkg, quietly = TRUE)
}

message("\n\n==== INSTALL REPORT ====")
for (pkg in names(status)) {
  message(if (status[[pkg]]) "  OK      " else "  FAILED  ", pkg)
}

failed <- names(status)[!status]
if (length(failed) > 0) {
  message("\nFAILED: ", paste(failed, collapse = ", "))
}

message("\nWriting renv.lock ...")
options(repos = BiocManager::repositories())
renv::snapshot(type = "all", prompt = FALSE)
message("Done. Lockfile: ", file.path(getwd(), "renv.lock"))

quit(status = if (length(failed) > 0) 1L else 0L)
