# Run both normalization routes independently, then compare variable features.

suppressPackageStartupMessages({
  library(Seurat)
})

source("R/norm_feat/io.R")
source("R/norm_feat/sct.R")
source("R/norm_feat/lognorm.R")
source("R/norm_feat/plots.R")

input_file <- file.path(
  "data", "clean_concatenated_data", "clean_concatenated.rds"
)
data_dir <- file.path("data", "norm_feat")
output_dir <- file.path("results", "norm_feat")
n_features <- 3000
seed <- 1234

obj <- read_normalization_input(input_file)

# --- SCT: normalization, variable feature selection and scaling in one step.
set.seed(seed)
obj_sct <- run_sctransform(
  obj,
  vars_to_regress = NULL,
  n_genes = n_features
)
saveRDS(obj_sct, file.path(data_dir, "sct.rds"))

# --- LogNormalize: restart from the same input, not from the SCT result.
set.seed(seed)
obj_log <- run_lognormalize(obj)
obj_log <- find_hvgs(obj_log, n_features = n_features)
obj_log <- scale_data(
  obj_log,
  features = VariableFeatures(obj_log),
  vars_to_regress = NULL
)
saveRDS(obj_log, file.path(data_dir, "lognorm.rds"))

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
