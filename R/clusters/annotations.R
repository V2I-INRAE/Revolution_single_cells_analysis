# Owner-approved coarse subsets for job-44493014, harmony_snn_res.0.7.
# Named vector: cluster ID -> subset label, for annotate_clusters_manual().
coarse_subsets <- c(
  `1` = "T cells & NK", `2` = "T cells & NK",
  `4` = "T cells & NK", `5` = "T cells & NK",
  `7` = "B cells", `8` = "B cells",
  `3` = "Myeloid", `9` = "Myeloid", `13` = "Myeloid",
  `14` = "Myeloid", `18` = "Myeloid",
  `6` = "Granulocytes", `15` = "Granulocytes",
  `11` = "Endothelial", `12` = "Endothelial", `20` = "Endothelial",
  `16` = "Epithelial", `17` = "Epithelial", `21` = "Epithelial",
  `23` = "Epithelial", `24` = "Epithelial",
  `10` = "Stromal", `19` = "Stromal",
  `22` = "Mixed"
)

# Owner-approved coarse subsets for job-44542861, harmony_snn_res.0.7.
# Cluster numbers are this run's own labels; verified 1-24 each exactly once
# against results/job-44542861/clusters_plots/markers/0.7/cluster_identities_res0.7.md.
coarse_subsets_44542861 <- c(
  `2` = "CD4_T_cells", `3` = "NK_cells", `4` = "unconventional_T_cells",
  `5` = "CD8_T_cells", `8` = "T_cells_unresolved",
  `7` = "B cells", `9` = "B cells",
  `1` = "Macrophage_Monocytes", `11` = "Macrophage_Monocytes",
  `14` = "Macrophage_Monocytes", `21` = "Macrophage_Monocytes",
  `6` = "Granulocytes", `15` = "Granulocytes",
  `16` = "DC", `19` = "DC",
  `12` = "Endothelial", `13` = "Endothelial",
  `17` = "Epithelial", `18` = "Epithelial", `22` = "Epithelial",
  `23` = "Epithelial", `24` = "Epithelial",
  `10` = "Stromal", `20` = "Stromal"
)

# Owner-approved coarse subsets for job-44890936, harmony_snn_res.0.4.
# Cluster numbers are this run's own labels (22 clusters); verified 1-22 each
# exactly once against
# results/job-44890936/clusters_plots/markers_res0.4/cluster_identities_res0.4.md.
# Working pre-subclustering groups; not final cell-type annotations.
coarse_subsets_44890936 <- c(
  `1` = "T_NK_ILC", `3` = "T_NK_ILC", `5` = "T_NK_ILC",
  `6` = "B_plasma", `7` = "B_plasma",
  `2` = "Myeloid", `10` = "Myeloid", `12` = "Myeloid", `20` = "Myeloid",
  `4` = "Granulocytes", `13` = "Granulocytes",
  `14` = "DC", `17` = "DC",
  `9` = "Endothelial", `11` = "Endothelial",
  `15` = "Epithelial", `16` = "Epithelial", `19` = "Epithelial",
  `21` = "Epithelial", `22` = "Epithelial",
  `8` = "Stromal", `18` = "Stromal"
)
