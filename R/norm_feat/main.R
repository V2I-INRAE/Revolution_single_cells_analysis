# Run both normalization routes independently, then compare variable features.

suppressPackageStartupMessages({
  library(Seurat)
})

source("R/norm_feat/params.R")
source("R/norm_feat/io.R")
source("R/norm_feat/sct.R")
source("R/norm_feat/lognorm.R")
source("R/norm_feat/pca.R")
source("R/norm_feat/plots.R")

input_file <- file.path(
  "data", "clean_concatenated_data", "clean_concatenated.rds"
)
data_dir <- file.path("data", "norm_feat")
output_dir <- file.path("results", "norm_feat")

dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

obj <- read_normalization_input(input_file)

# --- SCT: normalization, variable feature selection and scaling in one step.
set.seed(norm_feat_params$pca$seed)
obj_sct <- run_sctransform(
  obj,
  vars_to_regress = norm_feat_params$normalization$vars_to_regress,
  n_genes = norm_feat_params$normalization$n_variable_features,
  feature_buffer = norm_feat_params$normalization$sct_feature_buffer
)
obj_sct <- run_pca_analysis(
  obj_sct,
  n_pcs = norm_feat_params$pca$n_pcs
)

# --- LogNormalize: restart from the same input, not from the SCT result.
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
saveRDS(obj_sct, file.path(data_dir, "sct.rds"))

# --- Just to check overlap and distinct between the two methods
sct_features <- VariableFeatures(obj_sct, assay = "SCT")
log_features <- VariableFeatures(obj_log, assay = "RNA")
all_features <- union(sct_features, log_features)
shared_features <- intersect(sct_features, log_features)

feature_membership <- data.frame(
  gene = all_features,
  selected_sct = all_features %in% sct_features,
  selected_lognorm = all_features %in% log_features
)
write.csv(
  feature_membership,
  file.path(output_dir, "variable_feature_membership.csv"),
  row.names = FALSE
)

comparison <- data.frame(
  n_sct = length(sct_features),
  n_lognorm = length(log_features),
  n_shared = length(shared_features),
  n_sct_only = length(setdiff(sct_features, log_features)),
  n_lognorm_only = length(setdiff(log_features, sct_features)),
  n_union = length(all_features),
  jaccard = length(shared_features) / length(all_features)
)

write.csv(
  comparison,
  file.path(output_dir, "variable_feature_comparison.csv"),
  row.names = FALSE
)

plot_variable_features(
  obj_sct,
  output_dir = output_dir, filename = "variable_features_sct.png"
)
plot_variable_features(
  obj_log,
  output_dir = output_dir, filename = "variable_features_lognorm.png"
)
plot_hvg_overlap(
  obj_sct,
  output_dir = output_dir, filename = "hvg_overlap_sct.png"
)
plot_hvg_overlap(
  obj_log,
  output_dir = output_dir, filename = "hvg_overlap_lognorm.png"
)

route_objects <- list(sct = obj_sct, lognorm = obj_log)
for (route in names(route_objects)) {
  route_obj <- route_objects[[route]]
  plot_pca_elbow(
    route_obj,
    output_dir = output_dir,
    filename = paste0("pca_elbow_", route, ".png")
  )
  plot_pca_loadings(
    route_obj,
    output_dir = output_dir,
    filename = paste0("pca_loadings_", route, ".png")
  )
  for (group_by in c("sample", "pressure", "time", "pressure_time")) {
    plot_pca_grouping(
      route_obj,
      group_by = group_by,
      output_dir = output_dir,
      seed = norm_feat_params$pca$seed,
      filename = paste0("pca_", group_by, "_", route, ".png")
    )
  }
}
