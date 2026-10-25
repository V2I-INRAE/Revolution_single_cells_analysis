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
locked in `renv.lock`; the pipeline needs no Python environment.

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

## Pseudobulk pathway enrichment

`R/DEG/run_pseudobulk.R` runs **Hallmark and Reactome GSEA** for each eligible
cluster/contrast, writing separate `gsea/hallmark/` and `gsea/reactome/` outputs.
Supply the existing frozen `hallmark_pig.rds` and `reactome_pig.rds` directory;
no reference download or database update is performed. Both collections use
the same pig aliases, approved collision exclusions and signed edgeR rankings.
Reactome is MSigDB `C2:CP:REACTOME`; database version and reference hashes are
saved with results. BH correction is separate **within each collection,
cluster and contrast**, not pooled across collections. Reporting remains
`padj < 0.01`; unstarred terms on plots may be nonsignificant.

Each collection exports full NES/FDR/leading-edge tables, pathway mapping
coverage, ranked genes, dotplots/running-score plots and `gsea_results.rds`.
Historical Hallmark-only analyses remain readable; saved-figure regeneration
also includes Reactome when that collection's output directory exists.
Existing results are not backfilled automatically. For enrichment-only work,
`prepare_pseudobulk_gsea(..., collection = "reactome")` and
`run_pseudobulk_gsea()` accept the existing frozen references and saved
`pseudobulk_fit.rds`, respectively; use a new output directory and the original
run's approved exclusions. Do not rerun differential expression just to add a
collection.
