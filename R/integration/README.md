# Integration embeddings

New integration runs save UMAP and t-SNE for unintegrated PCA and the selected
method (scVI, Harmony or CCA). t-SNE uses the same input components as UMAP,
perplexity 30 and the saved run seed (1234 for current runs). Clustering remains
based on the high-dimensional representation, never on UMAP or t-SNE.

`setup_integration()` validates the input; it does not transform it. It checks
cell/metadata/PCA alignment, consecutive available PCs, sample labels, variable
features in PCA loadings and scaled expression, presence of RNA count layers, and one
normalized layer (or SCT model) per sample covering all cells exactly once.
CCA already uses this same validation before `run_cca_integration()`.

## Add t-SNE to completed runs

Run from the analysis root. This fits only t-SNE and writes new copies; original
checkpoints, UMAPs, integration representations, graphs, cluster assignments
and diagnostic results are not overwritten or recomputed.

```bash
sbatch scripts/sbatch_add_tsne.sh integration \
  data/integration/scvi/job-44281712/lognorm.rds \
  data/integration-tsne/scvi/job-44281712

sbatch scripts/sbatch_add_tsne.sh integration \
  data/integration/harmony/job-44281709/lognorm.rds \
  data/integration-tsne/harmony/job-44281709

sbatch scripts/sbatch_add_tsne.sh integration \
  data/integration/cca/job-44282203/lognorm.rds \
  data/integration-tsne/cca/job-44282203

sbatch scripts/sbatch_add_tsne.sh clustering \
  results/clustering/scvi/job-44282549/analysis.rds \
  data/clustering-tsne/scvi/job-44282549
```

Select explicit source runs. Output directories must be new. Integration
copies keep the original run-directory basename so the comparison reader's
provenance checks still apply. Clustering takes the completed `analysis.rds`
bundle, saves an augmented `lognorm.rds`, and copies the bundle with its
checkpoint path updated; assignments, settings and diagnostics are retained.
Repeat the clustering command for each completed method you want to plot,
including CCA and the unintegrated PCA baseline.

## Plot saved embeddings

After the embedding jobs finish, the existing R plotting entry points produce
both UMAP and t-SNE by default, without fitting embeddings. Submit explicit
checkpoint paths and a new output directory:

```bash
sbatch scripts/sbatch_integration_plots.sh \
  data/integration-tsne/scvi/job-44281712/lognorm.rds \
  data/integration-tsne/harmony/job-44281709/lognorm.rds \
  data/integration-tsne/cca/job-44282203/lognorm.rds \
  results/integration/comparison/tsne-scvi-harmony-cca

sbatch scripts/sbatch_clustering_plots.sh \
  results/clustering/comparison/tsne-scvi \
  data/clustering-tsne/scvi/job-44282549/analysis.rds \
  --plot-resolutions=0.4,0.6,0.8
```

For integration, optionally append `--embedding=umap`, `--embedding=tsne`, or
`--embedding=both`. Only the selected embeddings are loaded, checked and
plotted, so `--embedding=umap` also works with older UMAP-only checkpoints.
The same arguments work with `Rscript R/integration/plot_main.R`. The batch
script now takes paths instead of job IDs. Existing output directories are
rejected rather than overwriting figures. A plotting job can depend on all
three embedding jobs with `--dependency=afterok:<scvi_tsne_job>:<harmony_tsne_job>:<cca_tsne_job>`.

Integration plots retain all existing grouping/faceting options, with `umap_`
or `tsne_` filename prefixes and matching axis labels. Clustering saves
`umap_resolutions` and `tsne_resolutions` (PNG/SVG), plus one clustree.
The resolution filter remains figures-only. Clustering embedding figures are
written beside each supplied analysis bundle, as in the existing workflow.

Existing embeddings are reused only when their saved input dimensions and
settings match; incompatible embeddings stop rather than being overwritten.
Integration comparison still requires the same baseline data and coordinates
across the three methods. Its existing corrected-input provenance restriction
is unchanged. Mixing and silhouette diagnostics remain in the original
high-dimensional representations; they are not calculated on t-SNE.
