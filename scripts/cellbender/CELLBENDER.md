# CellBender CPU workflow

Run from `/work/project/revo-pig-sc/analysis`. Originals in `data/raw_data` are
read-only. Only unfiltered MEX archives are correction inputs. Filtered archives
are used for comparison of barcode calls.

Python scripts, this guide and the report patch live in `scripts/cellbender/`.
Slurm launchers remain directly under `scripts/`.

## Environment

CellBender 0.4.0 is added to the existing Python 3.12 `.venv-scvi`, retaining its
CPU-only PyTorch and all existing version pins:

```bash
uv pip install --python .venv-scvi/bin/python --torch-backend=cpu \
  --constraint requirements-scvi.txt --requirement requirements-cellbender.txt
uv pip check --python .venv-scvi/bin/python
```

The before/after package lists are in
`results/cellbender-preparation/_environment/`. Restoring with
`uv pip sync ... requirements-scvi.txt` alone removes CellBender; reinstall the
addition afterwards. Never mutate this environment while analysis jobs use it.

CellBender 0.4.0's HTML report indexes a string-indexed pandas Series with integer
labels where positions are intended. With pandas 3 this raises `KeyError` in
the latent-PCA plot. The report notebook also needs explicit inline plotting;
otherwise the batch job's `Agg` backend suppresses all HTML figures. Apply the
recorded two-line report-only fix once after a fresh installation, while no jobs
use the environment:

```bash
git apply --directory=.venv-scvi/lib/python3.12/site-packages \
  scripts/cellbender/cellbender-pandas3.patch
```

This does not change training, cell probabilities, or corrected counts. New runs
record pandas version and both report-source SHA-256 hashes in `run.json`.

## Preparation and parameter review

```bash
sbatch scripts/sbatch_cellbender_preparation.sh
```

`prepare_cellbender.py` writes `preparation.json`, `barcode_counts.tsv.gz`, and
`barcode_rank.png` in `results/cellbender-preparation/<sample>/`, and a root-level
`manifest.tsv` with absolute paths to the JSON files. Existing sample folders
are rejected rather than overwritten. The JSON records original file paths,
sizes, modification times, matrix checks, expected cells, and molecule counts at
selected ranks. Numeric-looking barcodes remain strings.

Expected cells comes from the metrics mRNA row. The plotted approximately 2x
candidate for total droplets is **not an approved parameter**. Inspect every
curve, choose a point into the background tail, and record the selection and
rationale in the JSON's `selected_total_droplets_included` and `parameter_review`
fields before running. For the reference sample REVO30-P4, 25,000 corresponds
to 303 molecules, with 13,422 BD putative cells.

## Pilot and bounded parallel jobs

```bash
.venv-scvi/bin/python scripts/cellbender/prepare_cellbender.py --samples REVO30-P4 \
  --review-total-droplets 25000 \
  --review-note 'Reference curve reviewed: 303 molecules at rank 25000'
sbatch scripts/sbatch_cellbender.sh REVO30-P4
```

The launcher requests 8 CPUs, 64 GB and 24 hours on `workq`; BLAS/OpenMP and
PyTorch threads are capped at the allocation. No CUDA is used. Every sample/run
has its own working and checkpoint directory. Inputs are extracted into
`data/cellbender/<sample>/<run_id>/<sample>_unfiltered_MEX/`; outputs go into
`results/cellbender/<sample>/<run_id>/`. `run.json` records parameters, versions,
command, elapsed time, and status. Existing run directories are rejected.

Inspect pilot QC and `sacct -j <job> --format=JobID,State,Elapsed,AllocCPUS,MaxRSS`
before submitting other samples. A reviewed array TSV has one `sample` column
with a header; each sample's reviewed parameters come from its preparation JSON.
Exclude the successful pilot.
After pilot QC acceptance, the approved starting concurrency is eight jobs;
consider twelve only after reviewing resources across samples:

```bash
sbatch --array=0-21%8 scripts/sbatch_cellbender.sh \
  results/cellbender-preparation/reviewed-array.tsv
```

The first REVO30-P4 run completed count estimation but failed report validation.
Its ELBO diagnostic recommends a half-learning-rate comparison. Keep it distinct
from the comparison run, with all other parameters unchanged:

```bash
sbatch scripts/sbatch_cellbender.sh REVO30-P4 --learning-rate 0.00005
```

## Validation and acceptance

`validate_cellbender.py` runs automatically after successful correction. It can
also be rerun inside an appropriate compute allocation with an explicit output
directory. It independently reads the raw input, checks all FPR H5s, reopens
filtered outputs with Scanpy, checks nonnegative integer counts and that no
counts were added, verifies cell probabilities/calls, and checks reports for
notebook errors. Per-gene and called-cell before/after tables and
`validation.json` are saved alongside outputs.

`barcode_group_validation.tsv` and `ambient_gene_by_barcode_group.tsv` separate
the original Rhapsody calls from additional CellBender calls at every FPR.
`called_cell_counts.tsv.gz` includes molecules and detected genes before/after
correction, with barcode-group membership. These are diagnostic comparisons,
not cell-type annotations or automatic acceptance of additional barcodes.

After applying the report fix, regenerate an existing successful CellBender run's
reports without retraining, then rerun validation inside Slurm:

```bash
.venv-scvi/bin/python scripts/cellbender/validate_cellbender.py \
  results/cellbender/REVO30-P4/job-44400520 --regenerate-reports
```

Original reports are archived under `report-regeneration-<job>/`; all six FPR H5
files are SHA-256 checked before and after regeneration. The repair and original
validation error are recorded in `run.json`. Numerical validation does not remove
ELBO warnings or establish that the extra cell calls are genuine.

The three FPRs produce separate suffixed files. After numerical validation,
unsuffixed H5/metrics/report symlinks select **FPR 0.01**. These links do not mean
QC acceptance: runs remain `validated_pending_qc` until diagnostics are reviewed.
Keep all checkpoints until QC and the FPR comparison are complete.

Inspect PDF/HTML reports for sensible priors, convergence, polarized cell
probabilities, and latent PCA structure. Record the review in `qc_review.json`
with `elbo_converged` (boolean), `qc_status`, and a review note; do not label a
failed run accepted. Reruns use new job directories, with explicit parameter
changes. The default learning-rate schedule requires retraining, not extending a
finished 150-epoch checkpoint, to increase epochs.

```bash
.venv-scvi/bin/python scripts/cellbender/summarize_cellbender.py
```

This writes `results/cellbender/summary.{tsv,json}` and includes every attempt,
pending sample, and failed run; it does not pick a latest run automatically.
Overlap uses three explicit denominators: BD calls retained, CellBender calls
also found by BD, and intersection/union (Jaccard). Ambient-associated gene
comparisons use exact names and record absent genes. Aggregate removal alone
does not establish preservation of cell-type-specific expression; that needs
validated labels on matching barcodes.
