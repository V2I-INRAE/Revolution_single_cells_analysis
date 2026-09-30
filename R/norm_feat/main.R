# LogNormalize the Scrublet-clean cells and plot feature/PCA diagnostics.

suppressPackageStartupMessages({
  library(Seurat)
})

source("R/norm_feat/params.R")
source("R/norm_feat/io.R")
source("R/norm_feat/lognorm.R")
source("R/norm_feat/pca.R")
source("R/norm_feat/plots.R")

input_file <- file.path(
  "data", "clean_concatenated_data", "clean_concatenated_scrublet.rds"
)
data_dir <- file.path("data", "norm_feat")
output_dir <- file.path("results", "norm_feat", "lognorm")

dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

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
for (group_by in c("sample", "pressure", "time", "pressure_time")) {
  plot_pca_grouping(
    obj_log,
    group_by = group_by,
    output_dir = output_dir,
    seed = norm_feat_params$pca$seed,
    filename = paste0("pca_", group_by, "_lognorm.png")
  )
}
