# Clustering comparison

Run from the analysis root. Each job reads an explicit checkpoint and clusters
all cells using one representation: `pca` (unintegrated), `harmony`, `integrated_cca`, or
`integrated_scvi`. It does not renormalize or rerun integration.
Neighbor graphs, silhouette diagnostics and UMAPs use components 1–20:
PCA/Harmony/CCA components for those routes and all 20 latent components for scVI.
Settings are in `R/clustering/params.R`.

```bash
sbatch scripts/sbatch_clustering.sh unintegrated data/integration/harmony/job-44215822/lognorm.rds
sbatch scripts/sbatch_clustering.sh harmony data/integration/harmony/job-44215822/lognorm.rds
sbatch scripts/sbatch_clustering.sh scvi data/integration/scvi/job-44227281/lognorm.rds
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
- `results/clustering/<method>/<run_id>/`: assignments, diagnostic cell IDs,
  per-cell silhouettes, cluster sizes, sample composition and summaries as CSV;
  `analysis.rds` also holds settings, source provenance and session information.

Existing run directories are rejected. `analysis.rds` is written last and is the
input to the plotting stage. A failed run may leave a directory but no completed
bundle. The last resolution becomes Seurat's active identity; it is **not** a
selected optimal clustering. Use the explicit resolution metadata columns.

Silhouettes use exactly 50,000 cells, selected proportionally by sample with
canonical barcode ordering and a fixed seed. Smaller inputs stop rather than
silently changing that number. One Euclidean distance vector is reused across
resolutions. This vector alone is about 10 GB, with additional copies during
silhouette calculation. The job requests 128 GB; inspect measured memory before
changing that request. Jobs run resolutions sequentially.

Singleton diagnostic clusters retain the standard silhouette value zero.
Absent clusters have no score; invalid partitions (one cluster or one cluster
per sampled cell) have NA scores and an explicit status. The summaries include
cell-weighted and cluster-weighted means, and report absent/singleton clusters.
Clusters below 10 full-data cells are flagged, never removed.

## Plot completed runs without recomputing

```bash
sbatch scripts/sbatch_clustering_plots.sh \
  results/clustering/comparison/my-comparison \
  results/clustering/unintegrated/job-<id>/analysis.rds \
  results/clustering/harmony/job-<id>/analysis.rds \
  results/clustering/scvi/job-<id>/analysis.rds \
  --plot-resolutions=0.4,0.6,0.8
```

Replace `<id>` with each clustering job ID. The comparison output directory must
be new. One or two completed methods can also be plotted. Baseline PCA/HVG and
metadata fingerprints, cell sets, settings and diagnostic IDs must match.
The optional last argument limits **figures only** to saved resolutions; omitting
it plots all resolutions. `plotted_resolutions.csv` records the figure scope.
Comparison CSVs and the source analysis bundles retain all resolutions, including
1.0, even when its UMAP panel, clustree row and comparison figure points are omitted.

New runs also save t-SNE using the same representation/components, seed 1234 and
perplexity 30. Resolution plots are exported for both embeddings with separate
`umap_` and `tsne_` prefixes; clustree and diagnostics are generated only once.
For older completed runs, first [add t-SNE to new checkpoint copies](../integration/README.md#add-t-sne-to-completed-runs)
without rerunning integration or clustering, then supply the copied `analysis.rds`
bundles to the plotting command above. Plotting does not fit missing embeddings.
