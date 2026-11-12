# Fixed settings of the paired pseudobulk DEG and GSEA pipeline (R/DEG/main.R).
deg_config <- list(
  # Pseudobulk eligibility
  min_cells = 20L, # cells per specimen within a cluster
  min_matched_pigs = 3L, # complete pig pairs required per contrast

  # edgeR::filterByExpr
  filter_min_count = 10,
  filter_min_total_count = 15,
  filter_large_n = 10,
  filter_min_prop = 0.7,

  # Differential expression (BH within cluster/contrast)
  deg_fdr = 0.05,
  deg_logfc = 0.5, # |log2FC| used only to colour the publication figure

  # GSEA (fgseaMultilevel; BH within collection/cluster/contrast)
  gsea_min_size = 15, # a gene set is tested only if 15-500 of its genes are ranked
  gsea_max_size = 500,
  gsea_nperm = 1000L,
  gsea_nperm_retry = 10000L, # whole contrast rerun once if any estimate is unavailable
  pathway_fdr = 0.05,
  seed = 1234,
  collections = c("hallmark", "reactome", "gobp"),

  # Frozen references
  reference_dir = "results/misc/msigdb_pig_2026.1",
  mapping_file = "results/misc/Ensembl_ID/gene-mapping/pig_to_human_orthologues_annotated.csv"
)
