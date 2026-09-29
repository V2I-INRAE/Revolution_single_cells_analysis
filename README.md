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

From the analysis directory, run `Rscript R/integration/main.R` inside a suitable
compute allocation. Full-data CCA memory and runtime need to be assessed before
submitting a production job; a small pilot does not establish those requirements.

The driver reads the existing `data/norm_feat/lognorm.rds` and `sct.rds`
sequentially, without renormalizing or changing sample layers/models. Each route
compares unintegrated PCA with the enabled sample-level integration methods using
dimensions 1–30. In `main.R`, `does_harmony` and `does_cca` control integration,
UMAP computation and diagnostics. Currently Harmony is enabled and CCA disabled.
Pressure and time are plotting variables, not correction variables. Corrected
representations are sensitivity analyses, not automatically preferred results.

Objects containing the selected UMAPs are saved to `data/integration/`. Comparison
plots display only UMAPs present in each object, separately for each route.
PNGs and diagnostic results are saved to `results/integration/`. Raw sample LISI
and silhouette scores use the same 20,000 cells sampled proportionally by sample
across both routes, with seed 1234 and LISI perplexity 30. The sampled IDs, per-cell
scores, summaries and run settings are retained. These metrics describe sample
mixing; without independent cell-type labels they do not establish preservation
of biology. Original assays, counts, metadata and PCA are retained.

## Getting data locally

If you want to work locally you can download the data using the scirpt `scripts/download_reads.sh`

>Use the sequential downloads if you download sequentially so that rsync opens just one ssh connection. For parallel download it is better to setup a ssh key with the cluster, otehrwise it will ask for inputing password each time it opens a ssh connection
 
