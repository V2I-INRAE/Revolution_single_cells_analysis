# REVOLUTION

Single-cell analysis of paired lung lobes. Run commands from
`/work/project/revo-pig-sc/analysis` on the cluster.

## Setup

Restore R dependencies with `renv::restore()`. Batch scripts load the required
modules; for interactive R, load `statistics/R/4.6.1` inside a compute allocation.

Create the separate Python environments once:

```bash
# QC / Scrublet
uv python install 3.11.16
uv venv --python 3.11.16 .venv-scrublet
uv pip sync --python .venv-scrublet/bin/python requirements-scrublet.txt

# scVI (CPU)
uv venv --python 3.12.13 .venv-scvi
uv pip sync --python .venv-scvi/bin/python --torch-backend=cpu requirements-scvi.txt
```

For interactive use, set `RETICULATE_PYTHON` to the chosen environment's Python
before starting a fresh R session. QC activates Scrublet automatically; scVI uses
the separate environment. See `renv.lock` and the requirements files for versions.

## Run the pipeline

These examples use existing checkpoints. For a new analysis, substitute the
upstream job IDs you intend to use; no stage automatically selects the latest run.

```bash
# QC: one entry point, sequential, Scrublet only; Rhapsody is the default
sbatch scripts/sbatch_qc.sh

# Select corrected counts instead
sbatch scripts/sbatch_qc.sh soupx
sbatch scripts/sbatch_qc.sh cellbender

# Normalize the selected Scrublet-clean checkpoint
sbatch scripts/sbatch_norm_feat.sh \
  data/clean_concatenated_data/job-44194877/clean_concatenated_scrublet.rds

# Integration: CCA is the working route
sbatch scripts/sbatch_integration_cca.sh data/norm_feat/job-44216449/lognorm.rds

# Optional reviewer/sensitivity runs: Harmony and scVI remain available
sbatch scripts/sbatch_integration.sh harmony data/norm_feat/job-44212906/lognorm.rds
sbatch scripts/sbatch_integration.sh scvi data/norm_feat/job-44212906/lognorm.rds

# Cluster a selected representation
sbatch scripts/sbatch_clustering.sh unintegrated data/integration/harmony/job-44215822/lognorm.rds
sbatch scripts/sbatch_clustering.sh harmony data/integration/harmony/job-44215822/lognorm.rds
sbatch scripts/sbatch_clustering.sh cca data/integration/cca/job-<integration_job_id>/lognorm.rds
```

For dependent submissions, use Slurm `--dependency=afterok:<job_id>[:<job_id>...]`.
Integration jobs also write six UMAP PNGs into their method/job results folder,
using the same in-memory object after saving the final RDS. Figure generation
lives in `R/integration/plots.R`; no separate plotting job is needed. Figures
show unintegrated PCA beside the selected method, grouped by sample, pressure,
time and pressure–time, plus sample facets. Diagnostic sampling/mixing metrics
are not run automatically. Concise `.INFO` records track saving and plotting;
a plotting failure is reported without losing the saved integration object.
See [clustering instructions](R/clustering/README.md) for clustering plots.
See [CCA marker instructions](R/find_markers/README.md) for all-cell marker
discovery and separate bar-plot, dot-plot and heatmap jobs.

## Outputs and tracking

`<run_id>` is `job-<SLURM_JOB_ID>`, or a timestamp plus PID outside Slurm.
`R/qc/main.R [rhapsody|soupx|cellbender]` and `scripts/sbatch_qc.sh` use one QC
workflow, defaulting to Rhapsody. Input loading in `R/qc/io.R` selects:

- **Rhapsody:** filtered MEX ZIPs under `data/raw_data/rhapsody/<sample>/`.
- **SoupX:** sparse corrected count matrices at
  `data/raw_data/soupx/<sample>/job-44409807/corrected_counts.rds`.
- **CellBender:** `<sample>_cellbender_FPR_0.01_filtered.h5` under
  `data/raw_data/cellbender/<sample>/`. Only the H5 expression matrix is read,
  not CellBender's latent groups. Source runs and unresolved convergence/cell-call
  warnings are recorded in `data/raw_data/cellbender/.INFO`.

Each run uses one source for all 23 samples, then concatenates, plots and saves
the Scrublet-clean checkpoint in the same job. Input source and paths are saved
in the objects' QC metadata and `results/qc/<run_id>/.INFO`; no source or run-ID
columns are added to cell metadata. No normalization or integration job is launched.
`scripts/sbatch_soupx.sh` remains the ambient-correction launcher.
Existing run directories are rejected, including on requeue. Flat RDS paths and
compatibility links have been removed.

| Stage | Checkpoints | Results |
|---|---|---|
| QC labeled | `data/qc_labeled_data/<run_id>/labeled_concatenated.rds` | `results/qc/<run_id>/` |
| QC clean | `data/clean_concatenated_data/<run_id>/clean_concatenated_scrublet.rds` | Same QC directory |
| Normalization | `data/norm_feat/<run_id>/lognorm.rds` | `results/norm_feat/lognorm/<run_id>/` |
| Integration | `data/integration/<method>/<run_id>/lognorm.rds` | `results/integration/<method>/<run_id>/` |
| Clustering | `data/clustering/<method>/<run_id>/lognorm.rds` | `results/clustering/<method>/<run_id>/` |

Logs are in `logs/`, with timestamps and job IDs for new runs. Check Slurm status
and logs: a saved checkpoint may precede a later failure. [RUN_HISTORY.md](RUN_HISTORY.md)
records historical jobs, inputs, outputs and outcomes.

## Methods

- **QC:** Scrublet is the sole doublet detector, run independently per sample
  on all input integer counts before any cell removal. A temporary normalized
  PCA/UMAP is used only for diagnostic figures; it does not alter saved counts.
  Scrublet uses scores >0.15 for final calls; its printed
  automatic threshold is not the applied cutoff. QC before/after plots compare
  all input cells with QC-only retained cells. Fixed QC thresholds are shared
  by all samples: 300–4,000 detected genes (inclusive), mitochondrial counts
  ≤10%, and log10(genes) / log10(UMIs) >0.8. No MAD filtering is used.
  Ribosomal percentage is calculated and plotted for diagnostics only; it does
  not control cell retention. Settings: `R/qc/params.R`.
- **Normalization:** LogNormalize, per-layer VST selection, 3,000 consensus
  variable genes and 50-PC PCA. SCT checkpoints are historical only.
  Settings: `R/norm_feat/params.R`.
- **Integration:** corrects sample batches, not pressure/time. CCA is the
  working route; Harmony and scVI remain available for reviewer/sensitivity
  runs. Unintegrated UMAP and Harmony/CCA use dimensions 1–20; scVI uses original RNA counts for 3,000 selected genes and
  20 latent dimensions. Settings: `R/integration/params.R`.
- **Clustering:** Leiden (`algorithm = 4`, `leidenbase` backend, modularity
  objective) on preserved PCA (unintegrated), Harmony or CCA, using dimensions
  1–20. Use newly generated 20-component integration checkpoints; the example job
  IDs above refer to historical runs. Resolution UMAPs and clustree plots are
  retained; numerical clustering diagnostics are no longer generated.
  Settings: `R/clustering/params.R`. Clusters are not validated cell identities.

Figures are PNG; saved objects retain run metadata. Historical shared-path
figures remain in place and may contain products from multiple jobs.

## Download raw reads

Use `scripts/download_reads.sh` from your local computer after checking its
source/destination paths. Sequential downloads reuse one SSH connection;
parallel downloads are best used with SSH keys.
