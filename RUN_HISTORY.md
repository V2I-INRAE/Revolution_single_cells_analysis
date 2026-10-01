# Pipeline run history

Audited 2026-10-01. Times are CEST; dates are in 2026. Status comes from Slurm
accounting; settings/inputs from logs and saved diagnostics. This is an inventory,
not a recommendation of methods or resolutions. Submission commands: [README](README.md).

## QC and normalization: surviving files

| Job | Run interval | Checkpoints | Status / log |
|---|---|---|---|
| 44194877 | Sep 30 00:06–02:18 | [Labeled](data/qc_labeled_data/job-44194877/labeled_concatenated.rds), [DoubletFinder clean](data/clean_concatenated_data/job-44194877/clean_concatenated_doubletfinder.rds), [Scrublet clean](data/clean_concatenated_data/job-44194877/clean_concatenated_scrublet.rds) | COMPLETED; [log](logs/20260930-0006-qc-pipeline.log) |
| 44212906 | Sep 30 13:25–13:46 | [LogNormalize](data/norm_feat/job-44212906/lognorm.rds) | COMPLETED; [log](logs/20260930-1325-norm-feat-pipeline.log) |
| 44176755 | Sep 29 01:43–02:55 | [Legacy SCT](data/norm_feat/job-44176755/sct.rds) | COMPLETED; [log](logs/20260929-0143-norm-feat-pipeline.log) |

Files were moved unchanged into job folders on Oct 1. Attribution uses logs,
accounting and modification times, not internal job metadata. Flat-path links
were removed. Earlier overwritten checkpoints, including 44176755's LogNormalize
output, are unavailable; historical shared-path figures may mix jobs.

## Integration

Job-specific checkpoints are `data/integration/<method>/job-<ID>/lognorm.rds`;
diagnostics are `results/integration/<method>/job-<ID>/`. All inputs were read from
the former `data/norm_feat/` flat paths, whose contents changed between runs.

| Job | Run interval | Method / input cells | Dimensions | Outcome / log |
|---|---|---|---|---|
| 44184472 | Sep 29 13:58–16:48 | Legacy Harmony; 262,936 | PCA/Harmony 30 | CANCELLED; no uniquely attributable output; [log](logs/20260929-1358-integration-pipeline.log) |
| 44190509 | Sep 29 17:17–18:30 | Legacy Harmony; LogNormalize + SCT, 262,936 each | PCA/Harmony 30 | COMPLETED; shared [LogNormalize](data/integration/lognorm.rds)/[SCT](data/integration/sct.rds) outputs; [log](logs/20260929-1717-integration-pipeline.log) |
| 44193502 | Sep 29 21:35–22:08 | [Harmony](data/integration/harmony/job-44193502/lognorm.rds); 262,936 | PCA/Harmony 30 | COMPLETED; [log](logs/20260929-213551-integration-harmony-44193502.log) |
| 44193503 | Sep 29 21:35–22:22 | [scVI](data/integration/scvi/job-44193503/lognorm.rds); 262,936 | PCA/scVI 30 | COMPLETED; [log](logs/20260929-213551-integration-scvi-44193503.log) |
| 44215817 | Sep 30 15:47–16:36 | [scVI](data/integration/scvi/job-44215817/lognorm.rds); 262,303 | PCA/scVI 30 | COMPLETED; [log](logs/20260930-154712-integration-scvi-44215817.log) |
| 44215822 | Sep 30 15:47–16:18 | [Harmony](data/integration/harmony/job-44215822/lognorm.rds); 262,303 | PCA/Harmony 30 | COMPLETED; [log](logs/20260930-154738-integration-harmony-44215822.log) |
| 44227281 | Sep 30 22:37–23:20 | [scVI](data/integration/scvi/job-44227281/lognorm.rds); 262,303 | PCA 30 / scVI 20 | COMPLETED; [log](logs/20260930-223724-integration-scvi-44227281.log) |

All runs used 23 samples and seed 1234. Harmony: sample correction, maximum 10
iterations. scVI: 3,000 selected RNA genes, original counts, two layers,
negative-binomial likelihood; logs show 30 epochs. Diagnostic sampling was
20,000 cells for legacy runs, 50,000 thereafter; LISI perplexity 30.
No completed CCA run is recorded.

### Integration plots

Outputs: `results/integration/comparison/scvi-job-<S>_harmony-job-<H>/`.
Each directory contains `comparison_sources.csv` and PNG/SVG figures.

| Plot job | Run interval | Inputs: scVI / Harmony | Outcome / log |
|---|---|---|---|
| 44215826 | Sep 30 16:36–16:45 | 44215817 / 44215822 | COMPLETED; output directory reused by 44220807; [log](logs/20260930-163630-integration-plots-44215826.log) |
| 44220807 | Sep 30 18:34–18:53 | 44215817 / 44215822 | COMPLETED; same output directory; [log](logs/20260930-183438-integration-plots-44220807.log) |
| 44227466 | Sep 30 23:20–23:46 | 44227281 / 44215822 | COMPLETED; [log](logs/20260930-232055-integration-plots-44227466.log) |

## Clustering

All nine jobs completed with 262,303 cells. U = unintegrated PCA, H = Harmony,
S = scVI. U/H inputs came from Harmony **44215822**; S inputs from scVI
**44215817** (A/B) or **44227281** (C).

Checkpoints: `data/clustering/<method>/job-<ID>/lognorm.rds`.
Results: `results/clustering/<method>/job-<ID>/`, including CSVs and `analysis.rds`
with settings and upstream provenance.

| Profile | Dimensions: PCA / Harmony / scVI | Neighbors | Resolutions |
|---|---|---|---|
| A | 30 / 30 / 30 | 20 | 0.4, 0.6, 0.8, 1.0 |
| B | 30 / 30 / 30 | 60 | 0.1, 0.2, 0.3, 0.4 |
| C | 30 / 30 / 20 | 100 | 0.05, 0.1, 0.15, 0.2, 0.25, 0.3 |

Common settings: Annoy Euclidean, 50 trees, SNN pruning 1/15, Louvain,
10 starts and 10 iterations, seed 1234, 50,000 diagnostic cells. Counts below follow
each profile's resolution order; they are not selected optimal cluster numbers.

| Job / checkpoint | Start | Method / profile | Cluster counts | Log |
|---|---|---|---|---|
| [44221464](data/clustering/unintegrated/job-44221464/lognorm.rds) | Sep 30 19:27 | U / A | 34/40/43/51 | [log](logs/20260930-192711-clustering-unintegrated-44221464.log) |
| [44221465](data/clustering/harmony/job-44221465/lognorm.rds) | Sep 30 19:27 | H / A | 33/37/44/52 | [log](logs/20260930-192741-clustering-harmony-44221465.log) |
| [44221467](data/clustering/scvi/job-44221467/lognorm.rds) | Sep 30 19:27 | S / A | 38/49/55/60 | [log](logs/20260930-192742-clustering-scvi-44221467.log) |
| [44222697](data/clustering/harmony/job-44222697/lognorm.rds) | Sep 30 21:16 | H / B | 19/24/28/30 | [log](logs/20260930-211633-clustering-harmony-44222697.log) |
| [44222698](data/clustering/unintegrated/job-44222698/lognorm.rds) | Sep 30 21:16 | U / B | 19/24/29/32 | [log](logs/20260930-211635-clustering-unintegrated-44222698.log) |
| [44222701](data/clustering/scvi/job-44222701/lognorm.rds) | Sep 30 21:16 | S / B | 26/33/34/37 | [log](logs/20260930-211654-clustering-scvi-44222701.log) |
| [44228330](data/clustering/scvi/job-44228330/lognorm.rds) | Sep 30 23:50 | S / C | 17/23/27/30/31/31 | [log](logs/20260930-235041-clustering-scvi-44228330.log) |
| [44228331](data/clustering/unintegrated/job-44228331/lognorm.rds) | Sep 30 23:50 | U / C | 15/19/21/24/25/25 | [log](logs/20260930-235042-clustering-unintegrated-44228331.log) |
| [44228335](data/clustering/harmony/job-44228335/lognorm.rds) | Sep 30 23:50 | H / C | 15/18/20/22/22/24 | [log](logs/20260930-235047-clustering-harmony-44228335.log) |

### Clustering plots

Outputs: `results/clustering/comparison/unintegrated-job-<U>_harmony-job-<H>_scvi-job-<S>/`.
Successful jobs also save method UMAP/clustree figures in each method's results folder.

| Plot job | Run interval | Input jobs: U / H / S | Outcome / log |
|---|---|---|---|
| 44221470 | Sep 30 19:54–19:56 | 44221464 / 44221465 / 44221467 | FAILED, 1:0: >45-colour guard; partial comparisons, no method UMAP/clustree figures; [log](logs/20260930-195423-clustering-plots-44221470.log) |
| 44222706 | Sep 30 21:57–22:04 | 44222698 / 44222697 / 44222701 | COMPLETED; [log](logs/20260930-215708-clustering-plots-44222706.log) |
| 44228340 | Oct 1 01:09–01:17 | 44228331 / 44228335 / 44228330 | COMPLETED; [log](logs/20261001-010945-clustering-plots-44228340.log) |

## Limits

- Completed jobs exited 0:0. File presence alone does not prove a complete run.
- Older jobs did not capture immutable input hashes or per-job source revisions;
  current settings must not be applied retroactively.
- Most multi-GB checkpoints were not deserialized for this audit; evidence came
  from accounting, logs, saved diagnostics/bundles and file checks.
- Sample-mixing metrics and exploratory clusters do not establish cell identities
  or biological preservation.
