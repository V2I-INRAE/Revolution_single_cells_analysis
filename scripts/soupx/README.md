# SoupX automatic baselines

Run from `/work/project/revo-pig-sc/analysis`:

```bash
sbatch scripts/sbatch_soupx.sh REVO30-P4

# All 23 samples, at most four running together:
sbatch --array=1-23%4 --output=logs/soupx-%A_%a.log \
  scripts/sbatch_soupx.sh scripts/soupx/samples.txt
```

Each task requests 4 CPUs, 16 GiB and 4 hours. It uses the project's existing
R/renv environment (SoupX 1.6.2) and requires no Python or GPU. Run IDs are
`job-<SLURM_JOB_ID>` for individual jobs, `job-<SLURM_ARRAY_JOB_ID>` shared across
samples in an array, or a timestamp/PID for direct execution. Array task IDs and
the sample manifest checksum are recorded per sample. Existing run directories
are rejected, not overwritten. Every sample runs independently: a failed sample
does not block other tasks and is never retried with silently changed parameters.

## Inputs and method

Inputs are the sample's unfiltered and filtered `RSEC_MolsPerCell` MEX ZIPs and
`Metrics_Summary.csv` under `data/raw_data/<sample>/`. Numeric BD barcode strings
are preserved. Shared gene IDs, names and called-cell counts must agree exactly;
genes absent from the filtered matrix must be zero in the raw called-cell columns.
The full raw feature universe is retained. Symbol repair is applied once on that
universe after unique-ID matching; duplicate feature IDs cause failure.

The default background range is strictly `0 < molecules < 100`. Its gate requires
the molecules at the Rhapsody called-cell rank to exceed 1,000, at least 2,000
non-called barcodes in range, and no called barcodes in range. On failure, stop
for review; do not silently change the range. These count checks support, but do
not prove, the background-only assumption.

All original called cells are clustered independently per sample: LogNormalize
(scale factor 10,000), 2,000 VST variable genes, scaling without regression,
30 computed PCs, PCs 1–20 for neighbors with k=20, and Louvain resolution 0.6.
Seed 1234. There is no integration and no QC exclusion. Only the preliminary
clustering copy is normalized; official SoupX receives original integer counts.

`SoupChannel` estimates the default background once. `setClusters` supplies the
complete named partition, then `autoEstCont` uses its package defaults and
`forceAccept=FALSE`. `adjustCounts(method="subtraction", roundToInt=TRUE)` uses
seed 1234 immediately before stochastic rounding.

## Outputs

- Processed data: `data/soupx/<sample>/<run_id>/`
  - `corrected_counts.rds`: sparse integer genes × original called cells.
  - `soupx_channel.rds`: fitted channel, including original called-cell counts.
  - `preliminary_clustering.rds`: labels, PCs, variable genes and Seurat commands.
  - `feature_map.tsv`: original feature IDs/symbols and shared SoupX names.
- Diagnostics: `results/soupx/<run_id>/<sample>/`
  - `run.json`, `sessionInfo.txt`, `background_check.json`.
  - Native SoupX `contamination_estimation.png`, marker/estimate/profile tables.
  - Cluster assignments, per-cell counts and gene counts by preliminary cluster.
  - Matched `markers_Raw.png` and `markers_SoupX.png`, plus
    `marker_expression_by_cluster.csv` and `missing_panel_genes.txt`.
- Log: `logs/soupx-<job_id>.log`, or `logs/soupx-<array_id>_<task_id>.log`.

Every successful baseline includes before/after marker plots using the same
cells, candidate genes, numeric cluster order and color limits. The 26-gene
candidate panel is shared with the pilot; each sample has its own preliminary
partition, not copied pilot cluster identities. Present/missing genes are recorded,
and counts, detection fractions and mean log1p expression are exported per cluster.
Plot-only LogNormalize uses each method's own cell totals; neither color nor
removal alone establishes biological accuracy. These plots support sample-level
review, not automatic cell-type annotation or acceptance. No sensitivity settings,
cell filtering, mitochondrial renaming or downstream metadata changes are applied.

Serialization validation reopens both saved count sources, checks exact IDs,
integer/nonnegative/finite values, and verifies that correction never adds counts.
Input checksums, settings, scripts' checksums, versions, warnings and elapsed time
are recorded. Errors flag `failed_requires_review` and cause a nonzero job exit.
An estimation failure retains the pre-estimation channel and available diagnostics.
`forceAccept` only overrides high-rho safeguards; it cannot repair missing markers.

Successful numerical checks yield `validated_pending_qc`, never automatic
acceptance. Inspect the native fit and marker/count tables for estimation support,
off-population reduction and preservation of genuine expression. Preliminary
clusters are not validated cell types. A rho outside 1–10% prompts review, not an
automatic parameter change. A target removal rate or correlation is not sufficient
validation. Definitive CellBender comparisons require an explicitly accepted run
and matched features/cells. A user-selected, numerically validated but unaccepted
run may be compared exploratorily; no CellBender attempt is selected automatically.

Review each sample before accepting it or choosing sensitivity settings.
Original raw archives and existing downstream checkpoints remain untouched.

## Pilot biological review and matched comparison

```bash
sbatch scripts/sbatch_soupx_review.sh data/soupx/REVO30-P4/job-44404963

# User-selected, numerically validated comparator; its QC status is preserved:
sbatch scripts/sbatch_soupx_review.sh data/soupx/REVO30-P4/job-44404963 \
  results/cellbender/REVO30-P4/job-44400520

# Completed lower-learning-rate comparator:
sbatch scripts/sbatch_soupx_review.sh data/soupx/REVO30-P4/job-44404963 \
  results/cellbender/REVO30-P4/job-44403564
```

`review_soupx.R` is specific to this pilot's fixed 21-cluster partition. It does
not refit correction, recluster, filter cells, or alter processed matrices.
Qualitative raw-marker patterns define provisional positive/negative reference
groups. Six unresolved groups and overlapping lineages are excluded from negative
references, not from the data. The mappings are post-hoc QC aids, not final cell
annotations or independent experimental validation. All 21 clusters remain visible.

The BD `cell_type_experimental.csv` contains predicted immune labels. The v2.2
[manufacturer's guide](https://www.bdbiosciences.com/content/dam/bdb/marketing-documents/products-pdf-folder/software-informatics/rhapsody-sequence-analysis-pipeline/Rhapsody-Sequence-Analysis-Pipeline-UG.pdf)
describes a classifier trained on human PBMCs. Those predictions are tabulated
but not used as a pig-lung biological reference.

Review outputs use `results/soupx/<review-run-id>/REVO30-P4/`. Tables report
per-gene counts, detection fractions and retention in provisional references.
Raw/SoupX marker dot plots use `scplotter::FeatureStatPlot`, a shared color scale,
and Seurat LogNormalize only on plotting copies. Dot size is the fraction of
cells expressing the gene (0–1). CellBender mode requires an
explicit numerically validated run; unaccepted comparators are labeled exploratory
in metadata and plots, not promoted to accepted. It uses FPR 0.01 on all Rhapsody calls,
matches original feature IDs, and reports correlations of both corrected totals
and removed totals. Counts come from the full H5, independently of CellBender's
cell-call decision: original cells it no longer calls are retained and explicitly
flagged in the cell table and agreement metadata. Extra CellBender calls are
never added to the comparison. No cell-quality filtering is applied here.

CellBender H5 input is read with its official Python `load_data` API via
reticulate and the existing `.venv-scvi` environment. Seurat's standard H5 reader
mistakes this format's latent-data groups for expression matrices. No source H5
is converted or modified; Python is needed only for this comparison reader.

Source and output checks do not replace inspection of the rendered figures or
method diagnostics. No other sample is automatically launched by this review.
