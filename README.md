# REVOLUTION

Analysis of the paired lungs lobes single cells experiment.

## Setting up

Our data and code is on the cluster @ `/work/project/revo-pig-sc`. A set of library have been installed. The script used for that is on `scripts` folder. When you `ssh` in the cluster:

```bash
srun -c 8 --mem=64G --time=08:00:00 --cpu=6 --pty bash -i
```

Move to `/work/project/revo-pig-sc` and run:

```bash
module load statistics/R/4.6.1
```

### Scrublet environment

From the analysis directory, create and activate the Python environment:

```bash
uv python install 3.11.16
uv venv --python 3.11.16 .venv-scrublet
source .venv-scrublet/bin/activate
uv pip sync requirements-scrublet.txt
export RETICULATE_PYTHON="$VIRTUAL_ENV/bin/python"
```

Set this before starting R. `scripts/sbatch_qc.sh` activates the same environment.
The R dependency is Moonerss/scrubletR, pinned in `renv.lock`; `renv::restore()`
restores it. `R/qc/scrublet.R` calls `scrubletR::scrublet_R()` with the sample-specific
BD expected rate and `threshold = NULL`, then classifies scores **> 0.15** as
doublets in R. This avoids the package's manual-threshold argument bug without
patching it. Automatic-threshold statistics printed by Python are not the final
0.15 classifications. The Python bridge uses seed 0.
Both methods are plotted before filtering. The pipeline applies the same QC flags
and filters doublets independently with each method, saving two merged objects in
`data/clean_concatenated_data/`:

- `clean_concatenated_doubletfinder.rds`
- `clean_concatenated_scrublet.rds`

Both retain the two methods' scores and calls; `doublet_filter_method` identifies
which method controlled removal. `doublet_class` always means DoubletFinder.
The minimum-three-cells gene filter runs separately in each branch, so gene sets
may differ. Any existing `clean_concatenated.rds` is left untouched.

## Integration comparison

From the analysis directory, select one integration method:

```bash
sbatch scripts/sbatch_integration.sh harmony
sbatch scripts/sbatch_integration.sh cca
sbatch scripts/sbatch_integration.sh scvi
# Direct invocation inside a suitable compute allocation:
Rscript R/integration/main.R harmony
```

The method is required and case-insensitive; extra or unknown arguments fail.
Full-data CCA/scVI memory and runtime need to be assessed before submitting a
production job; a small pilot does not establish those requirements.

All downstream integration methods read only `data/norm_feat/lognorm.rds`,
without renormalizing or changing sample layers. The `norm_feat` pipeline remains
unchanged, and existing SCT files/results are retained but are not processed by
the integration driver. Each invocation compares unintegrated PCA with exactly
the selected sample-level integration method using dimensions 1–30, including
only those two reductions in UMAP plots and diagnostics.
Pressure and time are plotting variables, not correction variables. Corrected
representations are sensitivity analyses, not automatically preferred results.

### scVI

Restore the separate CPU Python environment from the analysis directory:

```bash
uv venv --python 3.12.13 .venv-scvi
uv pip sync --python .venv-scvi/bin/python --torch-backend=cpu requirements-scvi.txt
export RETICULATE_PYTHON="$PWD/.venv-scvi/bin/python"
```

Use a fresh R session, not one already bound to `.venv-scrublet`. The R wrapper
dependency is SeuratWrappers, recorded in `renv.lock`. Select `scvi` as the method
to run the implemented helper. No GPU support is configured in this environment.

For the LogNormalize input, `run_scvi_integration()` uses its exactly 3,000 existing
RNA HVGs. The model receives **original RNA counts**,
never log-normalized values, corrected SCT counts or residuals. Missing selected
genes in any RNA count layer stop the run. Sample batches come from those layers;
the wrapper joins counts internally without modifying the saved RNA assay.
The output reduction is `integrated_scvi`, with UMAP `umap_scvi`.

scVI uses seed 1234 (including Python), 30 latent dimensions, two hidden
layers and negative-binomial likelihood. `max_epochs = NULL` uses scvi-tools'
automatic epoch limit; the helper accepts an explicit limit for pilots. The
wrapper returns an embedding, not a saved model or training history. Selected
genes and settings are retained in `obj@misc$scvi`. A short pilot checks
compatibility, not convergence or biological preservation.

Objects are saved as `data/integration/<method>/<run_id>/lognorm.rds`.
PNGs and diagnostics go to `results/integration/<method>/<run_id>/`, with
`lognorm` in their filenames and a run-specific `diagnostic_cells.csv`.
The run ID is `job-<SLURM_JOB_ID>` under Slurm, otherwise a timestamp plus PID.
Existing run directories are never reused: reruns/requeues with the same method
and job ID stop. A failed reservation can leave an empty directory; a checkpoint
can survive a later diagnostic failure. Neither indicates a completed run.
Logs include the method and job ID (PID outside Slurm), plus a timestamp;
Slurm's bootstrap log remains `logs/integration-pipeline-<job_id>.out`.

Raw sample LISI
and silhouette scores use the same 20,000 cells across selected reductions,
sampled proportionally by sample with seed 1234 and LISI perplexity 30. The sampled
IDs, per-cell scores, summaries and run settings are retained. When comparing
separate method runs, check the ordered cell/sample pairs in `diagnostic_cells.csv`
and input/version provenance, not just the seed. These metrics describe sample
mixing; without independent cell-type labels they do not establish preservation
of biology. Original assays, counts, metadata and PCA are retained.

## Getting data locally

If you want to work locally you can download the data using the scirpt `scripts/download_reads.sh`

>Use the sequential downloads if you download sequentially so that rsync opens just one ssh connection. For parallel download it is better to setup a ssh key with the cluster, otehrwise it will ask for inputing password each time it opens a ssh connection
 
