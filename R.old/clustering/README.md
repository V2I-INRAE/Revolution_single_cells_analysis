# Clustering

Run from the analysis root. Each job reads an explicit checkpoint and clusters
all cells using one representation: `pca` (unintegrated), `harmony`, or
`integrated_cca`. It does not renormalize or rerun integration.
Neighbor graphs and UMAPs use components 1–20. Clustering uses Leiden
(`algorithm = 4`) with the R `leidenbase` backend and modularity objective;
`leidenbase` is recorded in `renv.lock`. Settings are in `R/clustering/params.R`.

```bash
sbatch scripts/sbatch_clustering.sh unintegrated data/integration/harmony/job-44215822/lognorm.rds
sbatch scripts/sbatch_clustering.sh harmony data/integration/harmony/job-44215822/lognorm.rds
sbatch scripts/sbatch_clustering.sh cca data/integration/cca/job-<integration_job_id>/lognorm.rds
```

These are historical examples, not automatic latest-file selection. Replace
the integration job IDs with new runs using the 20-component settings; old
30-component UMAPs are rejected rather than silently reused or overwritten.
The unintegrated route uses only the preserved PCA from the Harmony checkpoint;
it does not use corrected coordinates. `data/norm_feat/job-44212906/lognorm.rds` is also a
valid unintegrated input. If it has no UMAP, the stage computes one once.
Saved UMAPs must document matching input reduction/dimensions and are reused.

Outputs are isolated by method and run ID (`job-<SLURM_JOB_ID>` or timestamp/PID):

- `data/clustering/<method>/<run_id>/lognorm.rds`: full object, original assays,
  saved embeddings, method-specific NN/SNN graphs and all resolution columns.
  `obj@misc$clustering` retains settings, partitions, provenance and session information.
- `results/clustering/<method>/<run_id>/cluster_assignments.csv`: assignments
  exported for tabular use. No separate analysis bundle is written.

Existing run directories are rejected. Use the full `lognorm.rds` directly for
clustering plots and marker analysis. Wait for successful job completion before
loading it: a failed or running job can leave a partially written checkpoint.
The last resolution becomes Seurat's active identity; it is **not** a selected
optimal clustering. Use the explicit resolution metadata columns.

Jobs run resolutions sequentially. The job still requests 128 GB; inspect
measured memory before changing that request.

## Plot completed runs without recomputing

```bash
sbatch scripts/sbatch_clustering_plots.sh \
  data/clustering/cca/job-<id>/lognorm.rds \
  --plot-resolutions=0.1,0.2,0.3,0.4
```

Replace `<id>` with the clustering job ID; use the corresponding `harmony` or
`unintegrated` path for those methods. Each plotting job reads one full checkpoint.
No comparison directory or cross-run comparison files are created.
The optional last argument limits **figures only** to saved resolutions; omitting
it plots all resolutions. The input path, selected resolutions and output directory
are recorded in the job log. The source object retains all saved resolutions,
even when their UMAP panels and clustree rows are omitted.

`umap_resolutions.png` and `clustree.png` are written to
`results/clustering/<method>/<run_id>/`, using the method and run ID stored in each
object, overwriting those figures if they already exist. Resolution plots use the
saved UMAP; plotting does not fit missing embeddings. Clustree shows transitions
between resolutions, not resampling stability.

## Explore CCA UMAP parameters with fixed clusters

```bash
# Pilot the last configuration, inspect it, then run the remaining fourteen.
sbatch --array=14 scripts/sbatch_clustering_umap_explore.sh \
  data/clustering/cca/job-44452483/lognorm.rds
sbatch --array=0-13%2 scripts/sbatch_clustering_umap_explore.sh \
  data/clustering/cca/job-44452483/lognorm.rds
```

Submit the remaining array only after the pilot succeeds, its four PNGs are
inspected, and its memory use fits the allocation (2 CPUs/48 GB per task).
An explicit test name can replace array selection; without either, the default
is `min.dist_0.5`. Defaults, ordered tests and the direct `RunUMAP` helper are
shared with integration exploration in `R/utils/umap_explore.R`.

This CCA-only entry point requires saved resolutions 0.1/0.2/0.3/0.4 and a baseline
matching the exploration reference settings. It recomputes only UMAP on the saved
CCA coordinates, never integration, neighbour graphs or Leiden. Each variant saves
four baseline-left/test-right PNGs under
`results/clustering/cca/<run_id>/umap_explore/<test_name>/`, alongside `.INFO`.
Existing variant directories are rejected; neither the source object nor the
ordinary clustering plots are overwritten. No object or tested coordinates are
saved, so later reuse of a tested embedding requires recomputation.

The baseline settings come from the checkpoint. Automatic baseline epochs are
recorded as automatic; consult the original producer log for the effective count.
For integration 44282203, `logs/20261001-234739-integration-cca-44282203.log` records
200 epochs, matching the explicit 200 used for exploration. Each task checks
finite/aligned coordinates and hashes the original CCA coordinates, baseline UMAP
and fixed label columns before and after fitting/rendering. Clustree is not repeated
because the assignments do not change. UMAP neighbour count is separate from the
clustering graph's `k`; visual separation does not establish biological improvement.
