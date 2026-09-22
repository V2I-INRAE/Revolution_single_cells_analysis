#!/usr/bin/env Rscript
# Install all R packages required by the scrnaseq-seurat-core-analysis workflow.
# Runs inside the renv project so packages land in the project library on /work,
# then records exact versions in renv.lock (snapshot type = "all").

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
    ggplot2 = "ggplot2",
    ggprism = "ggprism",
    patchwork = "patchwork",
    dplyr = "dplyr",
    harmony = "harmony",
    ashr = "ashr",
    svglite = "svglite"
  ),
  bioc = c(
    SingleR = "SingleR",
    celldex = "celldex",
    SingleCellExperiment = "SingleCellExperiment",
    muscat = "muscat",
    DESeq2 = "DESeq2",
    ComplexHeatmap = "ComplexHeatmap"
  ),
  github = c(
    DoubletFinder = "chris-mcginnis-ucsf/DoubletFinder",
    lisi = "immunogenomics/lisi",
    SeuratData = "satijalab/seurat-data",
    Azimuth = "satijalab/azimuth"
  )
)

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
      dependencies = TRUE,
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
renv::snapshot(type = "all", prompt = FALSE)
message("Done. Lockfile: ", file.path(getwd(), "renv.lock"))

quit(status = if (length(failed) > 0) 1L else 0L)
