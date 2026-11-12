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

## Pseudobulk DEG and pathway enrichment

Run paired differential expression and pathway enrichment for a saved subset from the project root. Replace the job IDs, subset label, and resolution with values from the input. The matching `DEG/` output directory must not already exist.

```bash
source_job=SOURCE_JOB
subset=SUBSET
subset_job=SUBSET_JOB
resolution=RESOLUTION

Rscript R/DEG/main.R \
  --input=data/job-${source_job}/subsets/${subset}/job-${subset_job}/subclustered.rds \
  --subset=${subset} --cluster-col=harmony_snn_res.${resolution} --clusters=all \
  --output-root=results/job-${source_job}/clusters_plots/subset/${subset}/job-${subset_job}/res${resolution}/DEG
```

The workflow excludes specimens with fewer than 25 cells and skips contrasts with fewer than three complete pig pairs. Each contrast uses its own paired edgeR model. It writes DEG tables, sample and model audits, and Hallmark, Reactome, and GO:BP GSEA results under `res0.X/DEG/`. Combined DEG plots go in the sibling `publication_figures/` directory. The workflow uses frozen pig pathway references in `results/misc/msigdb_pig_2026.1/` by default.

> Optional settings: `--pigs=<pig,pig>` limits the pig set. `--min-matched-pigs=<n>` raises the minimum pair count, which cannot be lower than three. `--references=<dir>` selects another frozen reference directory. `--gsea-exclusions=<CSV>` applies approved gene exclusions to GSEA. Run `Rscript R/DEG/main.R --help` for all options.

Each run writes its settings and status to `.INFO` files and its execution log under `logs/`. The source object and reference files are read-only. Existing result directories are not overwritten.
