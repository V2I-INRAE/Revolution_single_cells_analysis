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
# QC: labels all input cells and saves both clean doublet-filtering branches
sbatch scripts/sbatch_qc.sh

# Normalize the selected Scrublet-clean checkpoint
sbatch scripts/sbatch_norm_feat.sh \
  data/clean_concatenated_data/job-44194877/clean_concatenated_scrublet.rds

# Integration: choose harmony, scvi or cca
sbatch scripts/sbatch_integration.sh harmony data/norm_feat/job-44212906/lognorm.rds
sbatch scripts/sbatch_integration.sh scvi data/norm_feat/job-44212906/lognorm.rds

# Cluster a selected representation
sbatch scripts/sbatch_clustering.sh unintegrated data/integration/harmony/job-44215822/lognorm.rds
sbatch scripts/sbatch_clustering.sh harmony data/integration/harmony/job-44215822/lognorm.rds
sbatch scripts/sbatch_clustering.sh scvi data/integration/scvi/job-44227281/lognorm.rds

# Compare saved scVI, Harmony and CCA embeddings without refitting
# First add t-SNE to older checkpoints; new integration runs save both embeddings.
sbatch scripts/sbatch_integration_plots.sh \
  data/integration-tsne/scvi/job-44281712/lognorm.rds \
  data/integration-tsne/harmony/job-44281709/lognorm.rds \
  data/integration-tsne/cca/job-44282203/lognorm.rds \
  results/integration/comparison/my-comparison
```

For dependent submissions, use Slurm `--dependency=afterok:<job_id>[:<job_id>...]`.
See [integration plotting instructions](R/integration/README.md) for adding t-SNE
to completed runs and selecting `--embedding=umap|tsne|both` (default: both).
See [clustering instructions](R/clustering/README.md) for clustering plots.

## Outputs and tracking

`<run_id>` is `job-<SLURM_JOB_ID>`, or a timestamp plus PID outside Slurm.
Existing run directories are rejected, including on requeue. Flat RDS paths and
compatibility links have been removed.

| Stage | Checkpoints | Results |
|---|---|---|
| QC labeled | `data/qc_labeled_data/<run_id>/labeled_concatenated.rds` | `results/qc/<run_id>/` |
| QC clean | `data/clean_concatenated_data/<run_id>/clean_concatenated_{doubletfinder,scrublet}.rds` | Same QC directory |
| Normalization | `data/norm_feat/<run_id>/lognorm.rds` | `results/norm_feat/lognorm/<run_id>/` |
| Integration | `data/integration/<method>/<run_id>/lognorm.rds` | `results/integration/<method>/<run_id>/` |
| Clustering | `data/clustering/<method>/<run_id>/lognorm.rds` | `results/clustering/<method>/<run_id>/` |

Logs are in `logs/`, with timestamps and job IDs for new runs. Check Slurm status
and logs: a saved checkpoint may precede a later failure. [RUN_HISTORY.md](RUN_HISTORY.md)
records historical jobs, inputs, outputs and outcomes.

## Methods

- **QC:** DoubletFinder and Scrublet are labeled before filtering and produce
  separate clean branches. Scrublet uses scores >0.15 for final calls; its printed
  automatic threshold is not the applied cutoff. QC before/after plots compare
  all input cells with QC-only retained cells. Cells with ribosomal counts ≤2.5%
  are excluded. Settings: `R/qc/params.R`.
- **Normalization:** LogNormalize, per-layer VST selection, 3,000 consensus
  variable genes and 50-PC PCA. SCT checkpoints are historical only.
  Settings: `R/norm_feat/params.R`.
- **Integration:** corrects sample batches, not pressure/time. Unintegrated UMAP
  and Harmony/CCA use dimensions 1–20; scVI uses original RNA counts for 3,000 selected genes and
  20 latent dimensions. Settings: `R/integration/params.R`.
- **Clustering:** unintegrated uses preserved PCA; other routes use Harmony,
  CCA or scVI, selecting dimensions 1–20. Use newly generated 20-component integration
  checkpoints; the example job IDs above refer to historical runs.
  Settings: `R/clustering/params.R`. Mixing metrics do not establish
  biological preservation, and clusters are not validated cell identities.

Figures are PNG/SVG; saved objects retain run metadata. Historical shared-path
figures remain in place and may contain products from multiple jobs.

## Download raw reads

Use `scripts/download_reads.sh` from your local computer after checking its
source/destination paths. Sequential downloads reuse one SSH connection;
parallel downloads are best used with SSH keys.
