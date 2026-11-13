# REVOLUTION

Single-cell RNA-seq analysis of paired lung lobes (BD Rhapsody, 23 pig samples).
One pipeline run goes from raw counts to a clustered Seurat object: QC (including
scDblFinder doublet detection), normalization, integration and clustering in a
single session, configured by `R/config.R` and tracked in per-run `.INFO`
records.

## Setup

```bash
renv::restore()
```

Dependencies, including scDblFinder (Bioconductor) for doublet detection, are
locked in `renv.lock`; the pipeline needs no Python environment. Only the
Azimuth annotation step below uses Python, through its own uv environment.

## Usage

```bash
# 1. Set the run in R/config.R (input source, samples, integration method, thresholds)
# 2. Launch
sbatch scripts/sbatch_pipeline.sh
```

Outputs land in `data/job-<id>/` (`.INFO` + `clustered.rds`) and
`results/job-<id>/` (`.INFO` + `qc_plots/`, `norm_plots/`,
`integration_plots/`, `clusters_plots/`).

Ambient-correction and CellBender tooling live in `scripts/preprocessing/`.
Run history is recorded in [RUN_HISTORY.md](RUN_HISTORY.md).

## Cell explorer on the cluster

In a terminal on your allocated compute node, start the app:

```bash
bash /work/project/revo-pig-sc/analysis/app/run.sh
```

The script loads GCC 12.2.0 and R 4.6.1, selects the project directory and
prints the compute node name. The current data cache is already prepared.
Keep this terminal open. Press Ctrl+C to stop the app. If the app is already
running on that node, connect to it instead of starting another copy.

On your local computer, open a second terminal and create an SSH tunnel:

```bash
ssh -N -J fbraza@LOGIN_HOST \
  -L 3838:127.0.0.1:3838 fbraza@COMPUTE_NODE
```

Replace `LOGIN_HOST` with your usual Genotoul login hostname and `COMPUTE_NODE`
with the node printed by the launcher, for example `n005`. Leave this tunnel
open, then visit **http://127.0.0.1:3838** in your local browser.
The app listens only on the compute node's loopback address, not publicly.
See [app/README.md](app/README.md) for navigation and data preparation.

## Cell annotation with Pan-human Azimuth

Local Pan-human Azimuth (v1) predictions for any saved object, written to their own `azimuth_ph_v1_*` columns so other annotation models can be added alongside:

```bash
(cd R/annotations && uv sync)                    # once: Python 3.12 + panhumanpy
Rscript R/annotations/build_orthology_reference.R # once: Ensembl 116 pig-human orthology
Rscript R/annotations/azimuth.R --input=data/job-SOURCE_JOB/subsets/SUBSET/job-SUBSET_JOB/subclustered.rds
```

The model is human: pig features enter only through verified one-to-one Ensembl orthologues; unmapped model genes are zero-filled. Settings are in `R/annotations/config.R`. Output goes to `<input folder>/annotations/azimuth_ph_v1/`: `predictions.rds`/`.csv` (one row per cell barcode) and `gene_mapping.csv`; the input object is not modified. `False` in the broad/medium/fine columns marks an inconsistent hierarchy. Join the predictions by barcode to colour an existing UMAP, e.g. `obj$azimuth_ph_v1_medium <- p$azimuth_ph_v1_medium[match(colnames(obj), p$cell)]`.

For the full cleaned parent (234,183 cells), use the unchanged prediction script through Slurm, then the parent summary:

```bash
prediction_job=$(sbatch --parsable scripts/sbatch_azimuth_parent.sh)
sbatch --dependency=afterok:"$prediction_job" scripts/sbatch_azimuth_summary_parent.sh
```

`R/annotations/azimuth_summary_parent.R` produces a resolution-0.7 table and a two-panel figure (24 original clusters | Azimuth final labels), using the original parent UMAP restricted to cleaned cells. Outputs: `results/job-44542861/clusters_plots/res0.7/annotations/azimuth_ph_v1/`. The former subset summary is retained as `R/annotations/azimuth_summary_subset.R`; it is not run by the parent jobs.

## Pseudobulk DEG and pathway enrichment

Run paired differential expression and pathway enrichment for a saved subset from the project root. Replace the job IDs, subset label, and resolution with values from the input. Each run replaces the matching `res0.X/DEG/` directory and its `DEG_<contrast>.pdf` publication figures.

```bash
source_job=SOURCE_JOB
subset=SUBSET
subset_job=SUBSET_JOB
resolution=RESOLUTION

Rscript R/DEG/main.R \
  --input=data/job-${source_job}/subsets/${subset}/job-${subset_job}/subclustered.rds \
  --subset=${subset} --cluster-col=harmony_snn_res.${resolution} --clusters=all
```

Settings are in `R/DEG/config.R`. The workflow excludes specimens with fewer than 20 cells and skips contrasts with fewer than three complete pig pairs. Each contrast uses its own paired edgeR model. It writes DEG tables, sample and model audits, and Hallmark, Reactome, GO:BP, and KEGG GSEA results under `res0.X/DEG/`. Each `_gsea.csv` also records `collapsed_into` and `main_pathway` from `fgsea::collapsePathways` (significant terms explained by a more significant pathway).

Publication figures are made afterwards from the saved results, without refitting:

```bash
Rscript R/plots/publications/main.R \
  --input=data/job-${source_job}/subsets/${subset}/job-${subset_job}/subclustered.rds \
  --subset=${subset} --cluster-col=harmony_snn_res.${resolution}
```

It replaces, in `res0.X/publication_figures/`, the UMAP, the marker dot plot (`Markers/markers_top5_plus_canonical.csv`, colour capped at the 90th percentile), the top-25 marker heatmap (`Markers/markers_top25.csv`), one DEG plot per contrast, and one GSEA figure per collection and contrast with clusters on the x-axis (`gsea_hallmark/`, `gsea_gobp/`, `gsea_reactome/`). GSEA labels favour the terms in `R/plots/publications/gsea_label_whitelist.csv`. The workflow uses frozen pig pathway references in `results/misc/msigdb_pig_2026.1/`. KEGG uses a frozen snapshot of the native pig (`ssc`) pathways in `results/misc/kegg_pig/`, built once with `Rscript R/DEG/build_kegg_reference.R` (KEGG REST API, academic use only); rebuilding it takes the current KEGG release and requires a DEG rerun.

> Optional settings: `--clusters=1,2` runs selected clusters. `--output-root=<dir>/DEG` writes elsewhere, for example for a test run. `--gsea-exclusions=<CSV>` applies approved gene exclusions to GSEA.

The pipeline writes no `.INFO` or log files; redirect the console output to `logs/` if needed. The source object and reference files are read-only.
