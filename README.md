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

Run one saved subset locally from the project root, choosing a new output directory:

```bash
Rscript R/DEG/main.R \
  --input=data/job-44542861/subsets/DC/job-44831038/subclustered.rds \
  --subset=DC --cluster-col=harmony_snn_res.0.2 --clusters=1,2,3 \
  --output-root=results/job-44542861/clusters_plots/subset/DC/job-44831038/DEG_NEW
```

The driver retains the paired edgeR pseudobulk design, existing count filtering and GSEA ranking. It runs Hallmark, Reactome and GO:BP GSEA and ORA. The default frozen references are `results/misc/msigdb_pig_2026.1/` (MSigDB 2026.1.Hs pig projection); supply `--references=<directory>` to select an equivalent frozen reference explicitly. No database is downloaded by the analysis. Historical analyses may require their approved extra GSEA exclusions: supply `--gsea-exclusions=<CSV>` with `cluster,gene` columns. Unresolved GSEA identifier collisions block enrichment rather than being silently resolved.

ORA uses all genes with DEG FDR <= 0.05, without a fold-change cutoff, as a combined non-directional query. The background is the tested, unambiguous genes, including unannotated nonhits. Multiple tested features mapping to the same pig gene collapse to one background gene; any significant feature selects that gene. Ambiguous mappings and collapsed features are recorded. Terms with 15–500 background genes are tested with the hypergeometric upper tail, with BH adjustment per collection, cluster and contrast; report FDR <= 0.05.

Optional `--pigs=REVO26,REVO30` restricts animals; `--min-matched-pigs=2` permits an exploratory two-pig contrast (default: 3). Existing output directories are never overwritten. Detailed methods, exclusions, status and logs are recorded for each run. Source objects and references are read-only. Publication figures are also collected under the subset job's `publications_figures/` directory, with the analysis name and cluster preserved to avoid collisions. Existing results are not automatically regenerated.
