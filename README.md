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

## Getting data locally

If you want to work locally you can download the data using the scirpt `scripts/download_reads.sh`

>Use the sequential downloads if you download sequentially so that rsync opens just one ssh connection. For parallel download it is better to setup a ssh key with the cluster, otehrwise it will ask for inputing password each time it opens a ssh connection
 
