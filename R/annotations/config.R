# Settings for cell-level annotation models. Annotation outputs write no .INFO
# files; they are reproducible by rerunning with these settings.

annotation_config <- list(
  azimuth = list(
    # uv environment with panhumanpy: cd R/annotations && uv sync
    python = "R/annotations/.venv/bin/python",
    model_version = "v1",
    # Every model writes its own columns, so other annotation models can coexist.
    column_prefix = "azimuth_ph_v1_",
    output_folder = "azimuth_ph_v1"
  ),
  # Pig features enter a human model only through a verified one-to-one
  # Ensembl orthologue; symbol-named features must also keep the same symbol.
  # R/annotations/azimuth_summary_subset.R: one figure per subset/resolution and one CSV.
  # Resolutions are those retained for Markers/DEG/publication figures (Nov 6).
  summary = list(
    parent_job = "44542861",
    runs = data.frame(
      subset = c("DC", "Endothelial", "Epithelial", "Stromal",
        "Macrophage_Monocytes", "Macrophage_Monocytes"),
      job = c("44831038", "44831031", "44831032", "44831036", "44831034", "44831034"),
      resolution = c("0.2", "0.3", "0.3", "0.3", "0.2", "0.4")
    ),
    notes = "R/annotations/azimuth_notes.csv",
    # Figures only: labels below this fraction of a subset's cells are drawn as Other.
    min_label_fraction = 0.01,
    output_dir = "results/job-44542861/clusters_plots/subset/azimuth_ph_v1_summary"
  ),
  # Parent summary: retain the original parent UMAP and 24 clusters at res0.7,
  # restricted by barcode to the cleaned cells; never refit an embedding.
  parent_summary = list(
    input = "data/job-44542861/clean_annotated_object.rds",
    reference = "data/job-44542861/annotated_clustered.rds",
    parent_job = "44542861",
    input_producer_job = "44827588",
    resolution = "0.7",
    output_dir = "results/job-44542861/clusters_plots/res0.7/annotations/azimuth_ph_v1"
  ),
  orthology = list(
    file = "results/misc/Ensembl_ID/gene-mapping/sscrofa_hsapiens_orthology_ensembl116.tsv",
    biomart = "https://jun2026.archive.ensembl.org/biomart/martservice",
    ensembl_release = 116L
  )
)
