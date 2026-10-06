# Pipeline run history

Audited through **2026-10-02 22:53 CEST**. Times are CEST; dates are in 2026.
Status comes from Slurm accounting; settings and inputs from job logs, saved
diagnostics and clustering `analysis.rds` bundles. This is an inventory, not a
recommendation of methods or resolutions.
Submission commands: [README](README.md).

## QC and normalization: surviving files

| Job | Run interval | Checkpoints | Status / log |
|---|---|---|---|
| 44176755 | Sep 29 01:43–02:55 | [Legacy SCT](data/norm_feat/job-44176755/sct.rds) | COMPLETED; [log](logs/20260929-0143-norm-feat-pipeline.log) |
| 44194877 | Sep 30 00:06–02:18 | [Labeled](data/qc_labeled_data/job-44194877/labeled_concatenated.rds), [DoubletFinder clean](data/clean_concatenated_data/job-44194877/clean_concatenated_doubletfinder.rds), [Scrublet clean](data/clean_concatenated_data/job-44194877/clean_concatenated_scrublet.rds) | COMPLETED; [log](logs/20260930-0006-qc-pipeline.log) |
| 44212906 | Sep 30 13:25–13:46 | [LogNormalize](data/norm_feat/job-44212906/lognorm.rds) | COMPLETED; [log](logs/20260930-1325-norm-feat-pipeline.log) |
| 44276132 | Oct 1 17:25–19:32 | [Labeled](data/qc_labeled_data/job-44276132/labeled_concatenated.rds), [DoubletFinder clean](data/clean_concatenated_data/job-44276132/clean_concatenated_doubletfinder.rds), [Scrublet clean](data/clean_concatenated_data/job-44276132/clean_concatenated_scrublet.rds); [18 QC PNG/SVG files](results/qc/job-44276132/) | COMPLETED 0:0; [log](logs/20261001-172502-qc-pipeline-44276132.log) |
| 44280920 | Oct 1 19:53–20:24 | [LogNormalize RDS](data/norm_feat/job-44280920/lognorm.rds), [60 PNG/SVG files](results/norm_feat/lognorm/job-44280920/) | COMPLETED 0:0 **but not the intended two-covariate analysis**: `percent_mito` absent; only `nFeature_RNA` regressed. Preserved, not used below. [Log](logs/20261001-195329-norm-feat-pipeline-44280920.log) |
| 44281196 | Oct 1 20:36–21:14 | [Corrected LogNormalize RDS](data/norm_feat/job-44281196/lognorm.rds), [60 PNG/SVG files](results/norm_feat/lognorm/job-44281196/) | COMPLETED 0:0; **both** `percent.mt` and `nFeature_RNA` regressed. This is the input to the new integration and unintegrated clustering runs. [Log](logs/20261001-203622-norm-feat-pipeline-44281196.log) |

The **older** files were moved unchanged into job folders on Oct 1. Their
attribution uses logs, accounting and modification times, not internal job
metadata. Flat-path links were removed. Earlier overwritten checkpoints,
including 44176755's LogNormalize output, are unavailable; historical
shared-path figures may mix jobs. The new jobs above wrote run-specific paths.

**New upstream lineage and parameters.** QC 44276132 read 23 sample-specific
filtered MEX count archives under `data/raw_data/` (not a previous Seurat RDS).
It labeled 299,393 cells; 231,993 passed QC. DoubletFinder and Scrublet clean
Seurat RDS branches retain 228,077 and **226,738** cells, respectively, across
all 23 samples. Saved QC provenance and the [job log](logs/20261001-172502-qc-pipeline-44276132.log)
record: `nFeature_RNA > 200` and `<= min(4500, median + 4×MAD)` per sample;
`log10GenesPerUMI > 0.8`; `percent.mt <= min(20, median + 4×MAD)` with the
mitochondrial reference cells selected on feature-count and complexity criteria;
`percent.ribo > 2.5`; genes detected in **at least 10 retained cells per sample**.
ND4L is included once among 13 mitochondrial genes. Both clean cell sets were
checked against the labeled flags. The BD multiplet-rate lookup clamped values
outside its 0–4% table range; this is a modeling limitation, not a failed job.

Both normalization jobs took **the Scrublet clean Seurat RDS from 44276132** as
input, not the DoubletFinder branch. They used LogNormalize (scale factor
10,000), VST-selected 3,000 variable genes, scaling/regression, 50-PC PCA and
seed 1234. Job 44280920 requested `percent_mito` plus `nFeature_RNA`, but its
[log](logs/20261001-195329-norm-feat-pipeline-44280920.log) explicitly warns the
first name was missing from QC metadata and confirms only the second was
regressed. Job 44281196 corrected the name to `percent.mt`; its
[log](logs/20261001-203622-norm-feat-pipeline-44281196.log) confirms both
regressions, and missing covariates now stop the stage. Both outputs are Seurat
RDS files; only **44281196** feeds the subsequent runs. Historical input
objects are not identified solely by a mutable flat path.

## Integration

Job-specific checkpoints are `data/integration/<method>/job-<ID>/lognorm.rds`;
diagnostics are `results/integration/<method>/job-<ID>/`. **Earlier** inputs
were read from former `data/norm_feat/` flat paths, whose contents changed
between runs; the new runs pass the explicit `job-44281196` path.

| Job | Run interval | Method / input cells | Dimensions | Outcome / log |
|---|---|---|---|---|
| 44184472 | Sep 29 13:58–16:48 | Legacy Harmony; 262,936 | PCA/Harmony 30 | CANCELLED; no uniquely attributable output; [log](logs/20260929-1358-integration-pipeline.log) |
| 44190509 | Sep 29 17:17–18:30 | Legacy Harmony; LogNormalize + SCT, 262,936 each | PCA/Harmony 30 | COMPLETED; shared [LogNormalize](data/integration/lognorm.rds)/[SCT](data/integration/sct.rds) outputs; [log](logs/20260929-1717-integration-pipeline.log) |
| 44193502 | Sep 29 21:35–22:08 | [Harmony](data/integration/harmony/job-44193502/lognorm.rds); 262,936 | PCA/Harmony 30 | COMPLETED; [log](logs/20260929-213551-integration-harmony-44193502.log) |
| 44193503 | Sep 29 21:35–22:22 | [scVI](data/integration/scvi/job-44193503/lognorm.rds); 262,936 | PCA/scVI 30 | COMPLETED; [log](logs/20260929-213551-integration-scvi-44193503.log) |
| 44215817 | Sep 30 15:47–16:36 | [scVI](data/integration/scvi/job-44215817/lognorm.rds); 262,303 | PCA/scVI 30 | COMPLETED; [log](logs/20260930-154712-integration-scvi-44215817.log) |
| 44215822 | Sep 30 15:47–16:18 | [Harmony](data/integration/harmony/job-44215822/lognorm.rds); 262,303 | PCA/Harmony 30 | COMPLETED; [log](logs/20260930-154738-integration-harmony-44215822.log) |
| 44227281 | Sep 30 22:37–23:20 | [scVI](data/integration/scvi/job-44227281/lognorm.rds); 262,303 | PCA 30 / scVI 20 | COMPLETED; [log](logs/20260930-223724-integration-scvi-44227281.log) |
| 44281709 | Oct 1 22:05–22:32 | [Harmony](data/integration/harmony/job-44281709/lognorm.rds); **226,738** | PCA/Harmony **20** | COMPLETED 0:0; [diagnostics](results/integration/harmony/job-44281709/diagnostics_lognorm.rds), [log](logs/20261001-220536-integration-harmony-44281709.log) |
| 44281712 | Oct 1 22:05–23:06 | [scVI](data/integration/scvi/job-44281712/lognorm.rds); **226,738** | PCA/scVI **20** | COMPLETED 0:0; [diagnostics](results/integration/scvi/job-44281712/diagnostics_lognorm.rds), [log](logs/20261001-220558-integration-scvi-44281712.log) |
| 44282203 | Oct 1 23:47–Oct 2 22:51 (23:03:57) | CCA on **226,738** cells from [norm_feat 44281196](data/norm_feat/job-44281196/lognorm.rds); [checkpoint](data/integration/cca/job-44282203/lognorm.rds) | PCA/CCA **20** | COMPLETED 0:0; 8 CPUs, 512 GB, **72-hour** limit; [diagnostic RDS](results/integration/cca/job-44282203/diagnostics_lognorm.rds), [summary CSV](results/integration/cca/job-44282203/diagnostics_lognorm.csv), [metrics PNG](results/integration/cca/job-44282203/integration_metrics_lognorm.png), [log](logs/20261001-234739-integration-cca-44282203.log) |

All runs used 23 samples and seed 1234. Harmony: sample correction, maximum 10
iterations. scVI: 3,000 selected RNA genes, **original RNA counts** (not
log-normalized expression), two layers, negative-binomial likelihood,
`max_epochs=NULL` (automatic limit); the new scVI 44281712 log shows **35**
epochs, not the 30 observed in earlier runs. Diagnostic sampling was 20,000
cells for legacy runs, 50,000 thereafter; iLISI perplexity 30 and Euclidean
sample-label silhouettes. New Harmony/scVI checkpoint logs both identify the
exact [corrected norm_feat input](data/norm_feat/job-44281196/lognorm.rds):
PCA and each integrated UMAP/diagnostic use axes 1–20. The CCA script
[`sbatch_integration_cca.sh`](scripts/sbatch_integration_cca.sh) invokes the
same driver with CCAIntegration, LogNormalize, existing variable features,
original PCA 1–20 and integrated reduction `integrated_cca`. For completed CCA
44282203 the diagnostic RDS contains 50,000 cells each for `pca` and
`integrated_cca` on axes 1–20. Its log records completion of both UMAP fits;
the driver saves their `umap`/`umap_cca` coordinates in the checkpoint, but
does **not** produce a separate CCA UMAP image. The 7.0 GB checkpoint was not
independently deserialized on the memory-limited interactive runner, so its
embedded provenance and embeddings have not been directly inspected.

### Integration plots

Outputs: `results/integration/comparison/scvi-job-<S>_harmony-job-<H>/`.
Each directory contains `comparison_sources.csv` and PNG/SVG figures.

| Plot job | Run interval | Inputs: scVI / Harmony | Outcome / log |
|---|---|---|---|
| 44215826 | Sep 30 16:36–16:45 | 44215817 / 44215822 | COMPLETED; output directory reused by 44220807; [log](logs/20260930-163630-integration-plots-44215826.log) |
| 44220807 | Sep 30 18:34–18:53 | 44215817 / 44215822 | COMPLETED; same output directory; [log](logs/20260930-183438-integration-plots-44220807.log) |
| 44227466 | Sep 30 23:20–23:46 | 44227281 / 44215822 | COMPLETED; [log](logs/20260930-232055-integration-plots-44227466.log) |
| 44282037 | Oct 1 23:08–23:30 | **44281712 / 44281709** | COMPLETED 0:0; [source manifest](results/integration/comparison/scvi-job-44281712_harmony-job-44281709/comparison_sources.csv) records both Seurat checkpoint paths, [comparison output](results/integration/comparison/scvi-job-44281712_harmony-job-44281709/) contains seven nonempty PNG/SVG pairs (sample/pressure/time/pressure-time and three sample facets); [log](logs/20261001-230833-integration-plots-44282037.log) |

### Integration t-SNE augmentation and replotting: Oct 3

Submitted Oct 3 at 20:47–20:49 CEST. These jobs fit only t-SNE from the
completed integration representations; they do not rerun integration,
normalization, clustering or diagnostics. Each augmented copy retains the
original run-directory basename under `data/integration-tsne/<method>/`.
t-SNE uses components 1–20, seed 1234 and perplexity 30; saved UMAPs are retained.

| Job | Task / source integration | State checked at Oct 3 20:49 CEST | Log |
|---|---|---|---|
| 44310707 | Add baseline and scVI t-SNE from 44281712 | RUNNING | [log](logs/20261003-204712-add-tsne-44310707.log) |
| 44310708 | Add baseline and Harmony t-SNE from 44281709 | RUNNING | [log](logs/20261003-204712-add-tsne-44310708.log) |
| 44310709 | Add baseline and CCA t-SNE from 44282203 | RUNNING | [log](logs/20261003-204712-add-tsne-44310709.log) |
| 44310723 | Plot both UMAP and t-SNE for all three augmented checkpoints | PENDING; `afterok:44310707:44310708:44310709` | [Slurm output](logs/integration-plots-44310723.out) |

The plot job writes a new [comparison directory](results/integration/comparison/umap-tsne-scvi-job-44281712_harmony-job-44281709_cca-job-44282203/).
It is expected to produce 16 PNG/SVG pairs and `comparison_sources.csv`.
Submission and running states are not completion evidence; results remain
pending until job exit codes, logs and outputs are checked.

## Clustering

The **earlier nine** jobs (profiles A–C) completed with 262,303 cells. U =
unintegrated PCA, H = Harmony, S = scVI. For A–C, U/H inputs came from Harmony
**44215822** (U clustered its *preserved PCA*, not Harmony); S inputs came from
scVI **44215817** (A/B) or **44227281** (C). The six new jobs below (D/E)
process **226,738** cells; U reads norm_feat **44281196 directly**, H reads
Harmony **44281709**, and S reads scVI **44281712**. Those two integration
checkpoints also originate from norm_feat **44281196**.

Checkpoints: `data/clustering/<method>/job-<ID>/lognorm.rds`.
Results: `results/clustering/<method>/job-<ID>/`, including CSVs and `analysis.rds`
with settings and upstream provenance.

| Profile | Dimensions: PCA / Harmony / scVI | Neighbors | Resolutions |
|---|---|---|---|
| A | 30 / 30 / 30 | 20 | 0.4, 0.6, 0.8, 1.0 |
| B | 30 / 30 / 30 | 60 | 0.1, 0.2, 0.3, 0.4 |
| C | 30 / 30 / 20 | 100 | 0.05, 0.1, 0.15, 0.2, 0.25, 0.3 |
| D | **20 / 20 / 20** | **100** | **0.05, 0.1, 0.15, 0.2, 0.25, 0.3** |
| E | **20 / 20 / 20** | **100** | **0.4, 0.6, 0.8, 1.0** |

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
| [44282104](data/clustering/harmony/job-44282104/lognorm.rds) | Oct 1 23:23 | H / D; [input Harmony 44281709](data/integration/harmony/job-44281709/lognorm.rds) | 12/15/16/16/17/20 | COMPLETED 0:0; [bundle](results/clustering/harmony/job-44282104/analysis.rds), [log](logs/20261001-232327-clustering-harmony-44282104.log) |
| [44282105](data/clustering/unintegrated/job-44282105/lognorm.rds) | Oct 1 23:23 | U / D; [input norm_feat 44281196](data/norm_feat/job-44281196/lognorm.rds) (PCA) | 11/14/16/17/18/19 | COMPLETED 0:0; [bundle](results/clustering/unintegrated/job-44282105/analysis.rds), [log](logs/20261001-232328-clustering-unintegrated-44282105.log) |
| [44282108](data/clustering/scvi/job-44282108/lognorm.rds) | Oct 1 23:23 | S / D; [input scVI 44281712](data/integration/scvi/job-44281712/lognorm.rds) | 15/20/23/28/31/32 | COMPLETED 0:0; [bundle](results/clustering/scvi/job-44282108/analysis.rds), [log](logs/20261001-232346-clustering-scvi-44282108.log) |
| [44282546](data/clustering/harmony/job-44282546/lognorm.rds) | Oct 2 01:07 | H / E; [input Harmony 44281709](data/integration/harmony/job-44281709/lognorm.rds) | 20/24/28/30 | COMPLETED 0:0; [bundle](results/clustering/harmony/job-44282546/analysis.rds), [log](logs/20261002-010723-clustering-harmony-44282546.log) |
| [44282547](data/clustering/unintegrated/job-44282547/lognorm.rds) | Oct 2 01:07 | U / E; [input norm_feat 44281196](data/norm_feat/job-44281196/lognorm.rds) (PCA) | 21/24/27/30 | COMPLETED 0:0; [bundle](results/clustering/unintegrated/job-44282547/analysis.rds), [log](logs/20261002-010726-clustering-unintegrated-44282547.log) |
| [44282549](data/clustering/scvi/job-44282549/lognorm.rds) | Oct 2 01:07 | S / E; [input scVI 44281712](data/integration/scvi/job-44281712/lognorm.rds) | 33/38/42/**47** | COMPLETED 0:0; [bundle](results/clustering/scvi/job-44282549/analysis.rds), [log](logs/20261002-010739-clustering-scvi-44282549.log) |

For D/E, each small `analysis.rds` was read: the bundle records the exact input
path (and nested norm_feat path for H/S), identical 226,738 cell IDs, shared
baseline fingerprint prefix `0c5a330481738ebb`, settings, complete assignments
and 50,000 diagnostic IDs. Profiles D and E have respectively 300,000 and
200,000 sampled-cell silhouette rows *per method* (six or four resolutions),
with available scores. All six run-specific CSVs per job exist. These are
Seurat RDS checkpoints plus CSV/RDS **diagnostic outputs**, not annotated cell
types. Most full checkpoints could not be independently deserialized on the
memory-limited login runner; scVI E's checkpoint was independently read in a
separate 64 GB validation job **44282993**, which completed 0:0 and matched
assignments, 20-dimensional latent space and saved UMAP provenance.

### Clustering plots

Outputs: `results/clustering/comparison/unintegrated-job-<U>_harmony-job-<H>_scvi-job-<S>/`.
Successful jobs also save method UMAP/clustree figures in each method's results folder.

| Plot job | Run interval | Input jobs: U / H / S | Outcome / log |
|---|---|---|---|
| 44221470 | Sep 30 19:54–19:56 | 44221464 / 44221465 / 44221467 | FAILED, 1:0: >45-colour guard; partial comparisons, no method UMAP/clustree figures; [log](logs/20260930-195423-clustering-plots-44221470.log) |
| 44222706 | Sep 30 21:57–22:04 | 44222698 / 44222697 / 44222701 | COMPLETED; [log](logs/20260930-215708-clustering-plots-44222706.log) |
| 44228340 | Oct 1 01:09–01:17 | 44228331 / 44228335 / 44228330 | COMPLETED; [log](logs/20261001-010945-clustering-plots-44228340.log) |
| 44282442 | Oct 2 00:43–00:51 | **44282105 / 44282104 / 44282108** (D) | COMPLETED 0:0; [ordered input manifest](results/clustering/comparison/unintegrated-job-44282105_harmony-job-44282104_scvi-job-44282108/comparison_sources.csv), [comparison](results/clustering/comparison/unintegrated-job-44282105_harmony-job-44282104_scvi-job-44282108/) has 18 method-resolution metric rows, 33 agreement comparisons and 5,165 contingency rows plus three PNG/SVG figure pairs. Per-method UMAP/clustree PNG/SVG are in the three D result directories. [Log](logs/20261002-004357-clustering-plots-44282442.log) |
| 44283792 | Oct 2 06:57–07:04 | **44282547 / 44282546 / 44282549** (E) | COMPLETED 0:0; [ordered source manifest](results/clustering/comparison/unintegrated-job-44282547_harmony-job-44282546_scvi-job-44282549/comparison_sources.csv), [plot subset](results/clustering/comparison/unintegrated-job-44282547_harmony-job-44282546_scvi-job-44282549/plotted_resolutions.csv), [comparison](results/clustering/comparison/unintegrated-job-44282547_harmony-job-44282546_scvi-job-44282549/) with 12 method-resolution metric rows, 21 agreement comparisons and 6,699 contingency rows (each pair totals 226,738 cells); three nonempty comparison PNG/SVG pairs and two per-method PNG/SVG pairs (UMAP/clustree) for each E bundle. [Log](logs/20261002-065719-clustering-plots-44283792.log) |

For profile **E**, all three clustering jobs completed 0:0 and their saved
bundles passed source, baseline, settings, cell-set and diagnostic compatibility
checks. The scVI E bundle has **47** clusters at resolution 1.0, exceeding the
UMAP plotting palette's **45** colours. Faouzi chose to **omit 1.0 from plots
only**: the scoped UMAP/clustree and comparison figures show 0.4/0.6/0.8;
the saved clustering objects and comparison **CSV tables retain all four
resolutions**, including 1.0. The saved `plotted_resolutions.csv` identifies
the visual subset. The 1.0 partition was not removed from the analysis.

## Limits

- Unless marked failed/cancelled, listed production jobs exited 0:0. File
  presence alone does not prove a complete run.
- Older jobs did not capture immutable input hashes or per-job source revisions;
  current settings must not be applied retroactively.
- Most multi-GB checkpoints were not independently deserialized on the
  memory-limited login runner; evidence came from accounting, logs, saved
  diagnostics/bundles and file checks. scVI E is the stated high-memory exception.
- Sample-mixing metrics and exploratory clusters do not establish cell identities
  or biological preservation.
