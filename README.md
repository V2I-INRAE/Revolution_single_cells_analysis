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
