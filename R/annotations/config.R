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
  orthology = list(
    file = "results/misc/Ensembl_ID/gene-mapping/sscrofa_hsapiens_orthology_ensembl116.tsv",
    biomart = "https://jun2026.archive.ensembl.org/biomart/martservice",
    ensembl_release = 116L
  )
)
