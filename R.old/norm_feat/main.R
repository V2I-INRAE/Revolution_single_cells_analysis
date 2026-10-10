# LogNormalize the Scrublet-clean cells and plot feature/PCA diagnostics.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript R/norm_feat/main.R <clean_input.rds>")
}
input_file <- args[[1]]
if (!file.exists(input_file)) stop("Input checkpoint not found: ", input_file)
input_file <- normalizePath(input_file)

suppressPackageStartupMessages({
  library(Seurat)
})

source("R/norm_feat/params.R")
source("R/norm_feat/io.R")
source("R/norm_feat/lognorm.R")
source("R/norm_feat/pca.R")
source("R/norm_feat/plots.R")

job_id <- Sys.getenv("SLURM_JOB_ID")
run_id <- if (nzchar(job_id)) paste0("job-", job_id) else
  paste(format(Sys.time(), "%Y%m%d-%H%M%S"), Sys.getpid(), sep = "-")
data_dir <- file.path("data", "norm_feat", run_id)
output_dir <- file.path("results", "norm_feat", "lognorm", run_id)
for (directory in c(data_dir, output_dir)) {
  dir.create(dirname(directory), recursive = TRUE, showWarnings = FALSE)
  if (!dir.create(directory, showWarnings = FALSE)) {
    stop("Cannot create new run directory: ", directory)
  }
}
message("Normalization run: ", run_id, "\nInput: ", input_file,
  "\nData: ", data_dir, "\nResults: ", output_dir)

obj <- read_normalization_input(input_file)

set.seed(norm_feat_params$pca$seed)
obj_log <- run_lognormalize(
  obj,
  normalization_method = norm_feat_params$normalization$lognorm_method,
  scale_factor = norm_feat_params$normalization$lognorm_scale_factor
)
obj_log <- find_hvgs(
  obj_log,
  selection_method = norm_feat_params$normalization$hvg_selection_method,
  n_features = norm_feat_params$normalization$n_variable_features
)
obj_log <- scale_data(
  obj_log,
  features = VariableFeatures(obj_log),
  vars_to_regress = norm_feat_params$normalization$vars_to_regress
)
obj_log <- run_pca_analysis(
  obj_log,
  n_pcs = norm_feat_params$pca$n_pcs
)

obj_log@misc$norm_feat <- list(
  input_file = input_file, run_id = run_id,
  params = norm_feat_params, session_info = sessionInfo()
)
saveRDS(obj_log, file.path(data_dir, "lognorm.rds"))
plot_variable_features(
  obj_log,
  output_dir = output_dir, filename = "variable_features_lognorm.png"
)
plot_hvg_overlap(
  obj_log,
  output_dir = output_dir, filename = "hvg_overlap_lognorm.png"
)

plot_pca_elbow(
  obj_log, output_dir = output_dir, filename = "pca_elbow_lognorm.png"
)
plot_pca_loadings(
  obj_log, output_dir = output_dir, filename = "pca_loadings_lognorm.png"
)
