# Pipeline run history

Audited through **2026-10-02 22:53 CEST**. Times are CEST; dates are in 2026.
Status comes from Slurm accounting; settings and inputs from job logs, saved
diagnostics and clustering `analysis.rds` bundles. This is an inventory, not a
recommendation of methods or resolutions.

SVG copies of plot outputs were removed on Oct 7. Historical references below
to PNG/SVG pairs describe the original outputs; PNG copies remain.

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

### QC parameter record by run

These are execution-specific settings; they must not be inferred from the
current QC source files. A ribosomal rule changed between runs.

| QC job | Input / doublet route | Cell and mitochondrial QC | Ribosomal treatment | Evidence / limitation |
|---|---|---|---|---|
| 44194877 | Historical raw filtered MEX; DoubletFinder and Scrublet branches | `nFeature_RNA > 200`; the log records per-sample feature ceilings capped at 5,000; `log10GenesPerUMI > 0.8`; per-sample mitochondrial upper ceilings were applied | The log labels cells below versus at/above 2.5%, but its `qc_outlier` summary has no ribosomal flag. A 2.5% ribosomal filter is therefore not confirmed for this run. | [log](logs/20260930-0006-qc-pipeline.log); the exact historical MAD multipliers are not recoverable from retained run metadata. |
| 44276132 | 23 raw filtered MEX archives; DoubletFinder and Scrublet branches | `nFeature_RNA > 200` and `<= min(4500, median + 4×MAD)`; `log10GenesPerUMI > 0.8`; `percent.mt <= min(20, median + 4×MAD)` after feature/complexity reference-cell selection; genes in at least 10 retained cells per sample | Cells with `percent.ribo <= 2.5` were flagged and excluded. | [log](logs/20261001-172502-qc-pipeline-44276132.log); saved provenance. |
| 44415111 | 23 SoupX-corrected count matrices from `job-44409807`; Scrublet only | `nFeature_RNA > 200` and `<= min(4500, median + 4×MAD)`; `log10GenesPerUMI > 0.8`; `percent.mt <= min(20, median + 4×MAD)`; genes in at least 10 retained cells per sample | Cells strictly above each sample's `median(percent.ribo) + 3×MAD` were flagged; all input cells, including zeros, contributed. No lower cutoff or cap was used. The logged upper ceilings range from 18.36% to 26.80%. | [log](logs/20261007-171743-qc-with-soupx-44415111.log). |

For 44276132 and 44415111, Scrublet used original counts on all input cells,
with the retained-cell call threshold fixed at score >0.15 and seed 1234.

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
| 44417308 | Oct 7 21:29–21:31 (00:01:57) | CCA on **247,672** cells from [norm_feat 44416449](data/norm_feat/job-44416449/lognorm.rds) | PCA/CCA **20** | CANCELLED by user (scancel) before integration computed; no checkpoint; run directories (incl. diagnostic_cells.csv) removed at user request same day; [log](logs/20261007-212915-integration-cca-44417308.log) |

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

### Retired integration embedding augmentation and replotting: Oct 3

The temporary augmented copies and their comparison directory were removed when
the analysis was made UMAP-only. The job records and logs are retained below as
historical provenance; they are not active workflow instructions.

| Job | Task / source integration | State checked at Oct 3 20:49 CEST | Log |
|---|---|---|---|
| 44310707 | Retired augmented scVI checkpoint from 44281712 | RUNNING | [log](logs/20261003-204712-add-tsne-44310707.log) |
| 44310708 | Retired augmented Harmony checkpoint from 44281709 | RUNNING | [log](logs/20261003-204712-add-tsne-44310708.log) |
| 44310709 | Retired augmented CCA checkpoint from 44282203 | RUNNING | [log](logs/20261003-204712-add-tsne-44310709.log) |
| 44310723 | Retired comparison of the augmented checkpoints | PENDING; `afterok:44310707:44310708:44310709` | [Slurm output](logs/integration-plots-44310723.out) |

The associated augmented comparison directory was removed. Submission and
running states are not completion evidence; results remain pending until job
exit codes, logs and outputs are checked.

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

## CCA marker discovery: Oct 6

Input: [CCA clustering 44310241](results/clustering/cca/job-44310241/analysis.rds),
whose checkpoint contains 226,738 cells and 20,916 RNA features. Its seven
resolutions are 0.1/0.2/0.3/0.4/0.6/0.8/1.0 (16/21/25/25/29/33/36 clusters).
The source clustering job completed 0:0. Preflight job **44394946** independently
read the full checkpoint, confirmed 23 sample-specific LogNormalize data layers,
and completed a 50-gene, all-cell API test (0:0, 5m10s, peak RSS about 41 GiB).
That restricted-feature test is software validation, not the marker analysis.

| Job | Task | State checked Oct 6 |
|---|---|---|
| 44395117 | Real-data validation and figure checks | COMPLETED 0:0; 13m09s, peak RSS about 41 GiB |
| 44395220 | All-cell marker discovery, all seven resolutions | CANCELLED at user's request; [log](logs/20261006-154439-find-markers-44395220.log) |
| 44395226 | Bar plots, dot plots and all-cell heatmaps | CANCELLED at user's request |
| 44396189 | Production plot check | CANCELLED at user's request |

Validation checked normalized values against original counts, preserved values
across layer joining, all-cell counts, direct `FindMarkers` agreement, independently
calculated fold changes, Bonferroni adjustment, ranking, empty-marker clusters,
failed-test propagation, CSV/RDS round trips and source-file preservation.
All three figure types were rendered and inspected using the restricted-feature
software-test results. These figures are not the production marker results.

Production output: `results/find_markers/cca/job-44395220/res_*`.
Discovery uses RNA log-normalized data, Wilcoxon/Presto, positive markers,
`min.pct=0.25`, `logfc.threshold=0.25`, no detection-difference filter or cell cap,
and retains returned results before selecting Bonferroni-adjusted p-values <0.05.
Top markers are ranked by decreasing log2FC. No normalization, integration,
clustering, cell subsampling or annotation is rerun. The separate plotting job
uses saved marker tables and the original checkpoint; fresh unregressed scaling
is confined to the displayed genes. Both production jobs request 128 GB and
run through the [marker batch scripts](R/find_markers/README.md).

The production run was stopped for code review; partial outputs are not a
completed analysis. The subsequent code simplification removes marker RDS
bundles and retains per-resolution CSVs. Earlier validation applies to the
earlier implementation, not the revised code. No new analysis or plotting run
is authorized pending review.

## SoupX pilot: Oct 7

REVO30-P4 uses its original unfiltered/filtered RSEC MEX archives and metrics
under `data/raw_data/REVO30-P4/`. The pilot keeps all 13,422 Rhapsody calls;
no existing QC-cleaned checkpoint or CellBender-corrected counts are inputs.
The full raw gene universe is retained after checking that omitted filtered
rows are zero in called cells. Background starts at `0 < molecules < 100`,
conditional on a >1,000-molecule called-cell-rank proxy, at least 2,000
non-called background barcodes, and no called barcodes in that range.

Approved preliminary clustering: per-sample LogNormalize, 2,000 VST variable
genes, 30 computed PCs, PCs 1–20 for k=20 neighbors, Louvain resolution 0.6,
seed 1234, no integration or regression. SoupX 1.6.2 uses automatic estimation
without forced acceptance, subtraction, and seeded integer rounding. Existing
pipeline checkpoints remain unchanged. Numerical success is not QC acceptance.

| Job | Task | Last verified state |
|---|---|---|
| 44404957 | First pilot attempt | FAILED 1:0, 28s; input reader closed an already-closed gzip connection, before clustering/correction |
| 44404963 | Pilot after compressed-input connection fix | COMPLETED 0:0; 1m35s, peak RSS 5,432,652 KiB; numerically validated, pending biological QC |

The reader failure was reproduced and fixed by explicitly opening the gzip
connection before `Matrix::readMM`. ZIP/gzip reader, asymmetric alignment,
zero-padding rejection, integer/no-added-count validation, background-boundary
and actual metrics parsing tests pass in `scripts/soupx/test_io.R`. The launcher
passes Bash syntax and ShellCheck. Jobs request 4 CPUs, 16 GiB and 4 hours.

Scripts: `scripts/soupx/`; launcher: `scripts/sbatch_soupx.sh`.
Processed data: `data/soupx/REVO30-P4/job-44404963/`.
Diagnostics: `results/soupx/job-44404963/REVO30-P4/`.
Log: `logs/soupx-44404963.log`. The failed attempt is retained separately.
Review the pilot's fit and gene/count tables before other sample submissions;
no accepted CellBender comparator has been selected.

The background gate passed: 1,294 molecules at rank 13,422, with 134,703
non-called barcodes in range and no called barcodes in range. Preliminary
clustering yielded 21 groups. The fit used 100 marker genes and 697 independent
estimates; rho was 0.031 with half-maximum interval 0.015–0.064. The native
estimation PNG was inspected and showed a readable, unimodal posterior.
Counts changed from 84,811,621 to 82,183,607 molecules (3.10% removed), retaining
28,246 features and all 13,422 cells. Saved count matrices were reopened and
validated. Independently read cell/cluster CSV totals agree with those values.
The sole analysis warning is SoupX's deprecated Matrix `giveCsparse` argument;
no package patch was applied. At this stage marker preservation and an accepted
matched CellBender comparison remained open; see the subsequent review below.

### SoupX biological QC and user-selected CellBender comparison

The user selected `results/cellbender/REVO30-P4/job-44400520` for the first
comparison. It is numerically validated but remains unaccepted because of
convergence warnings and unresolved additional-cell calls. Its QC status is
unchanged; the comparison is exploratory, not accepted-comparator validation.

| Review job | Result |
|---|---|
| 44405241 | FAILED 1:0; plot ordering arguments were expression strings, not label vectors; removed unnecessary ordering overrides |
| 44405362 | FAILED 1:0; Seurat's H5 reader interpreted CellBender latent groups as matrices; switched to official CellBender `load_data` |
| 44405494 | FAILED 1:0; `normalizePath` resolved the venv Python symlink to the base interpreter; preserved the absolute venv path instead and verified the import |
| 44405863 | COMPLETED 0:0; 1m30s, peak RSS 5,578,044 KiB; tables independently checked and all three marker figures inspected |

Scripts: `scripts/soupx/review_soupx.R`, `scripts/sbatch_soupx_review.sh`.
Completed review: `results/soupx/job-44405863/REVO30-P4/`.
Failed review folders/logs are retained but are not completed outputs.

All 28,246 features and 13,422 original called cells match. Raw totals are
84,811,621 molecules; SoupX retains 82,183,607 (3.09865% removed), and
CellBender FPR 0.01 retains 81,406,722 (4.01466% removed). Corrected per-cell
totals correlate at r = 0.9998357856871326, whereas removed per-cell totals
correlate at r = −0.5740135066340926. High corrected-total agreement therefore
does not establish equivalent corrections. Cell/gene totals match the original
SoupX table and CellBender validation; 9,426 additional CB calls are excluded.

Pilot marker QC passed with limitations. Across 26 provisional markers,
SoupX retains 98.55–100% of positive-reference molecules per gene (99.80%
weighted) and removes 35.40% of negative-reference molecules in aggregate.
CellBender retains 99.61–100% (99.77% weighted) and removes 40.01% in negative
references. More removal is not evidence of better correction. Four markers
(C1QA, CSF3R, DCN, SFTPC) overlap rho-estimation markers; the other 22 give
similar SoupX results: 99.81% positive retention and 35.26% negative removal.
Post-hoc group definitions and six unresolved clusters limit interpretation.
BD human-PBMC-trained predictions are not pig-lung ground truth.

Raw, SoupX and exploratory CellBender dot plots show all 21 fixed clusters
and 26 genes on a common expression scale, without clipping. Strong marker
patterns remain. Quantitative marker counts and detection changes were checked
independently from CSVs; figures alone do not establish preservation.
`qc_review.json` now records the completed, qualified pilot review. Original
execution `run.json` is retained unchanged as the numerical-run record.
No raw or corrected matrices, CellBender files or downstream checkpoints were
modified. No other samples were submitted. Changes remain local/uncommitted.

### Comparison rerun with completed CellBender job 44403564

The user requested comparison with the latest completed lower-learning-rate
CellBender run (5e-05, 150 epochs). Slurm and numerical validation passed;
its existing reviewed convergence diagnostics passed without report warnings.
Cell-call/biological acceptance remains pending and its QC files were not changed.

The primary comparison retains the same 13,422 Rhapsody cells, 28,246 genes,
fixed raw clusters and post-hoc marker panels. No cell-quality filtering or
doublet removal was added. The full H5 contains all original barcodes; one
original cell, 55584626, is no longer called by CellBender and has zero corrected
molecules (raw 3,797; SoupX 3,721). It is retained and explicitly flagged rather
than silently dropped. The 5,419 additional CellBender calls are excluded.

| Review job | Result |
|---|---|
| 44408133 | FAILED 1:0 at rendering; passing only dot_size_name triggered partial list matching to dot_size upstream, treating the label as a function |
| 44408206 | COMPLETED 0:0; 1m27s, peak RSS 5,723,448 KiB; explicit package-default dot-size function avoids partial matching and preserves detection fractions |

Completed output: `results/soupx/job-44408206/REVO30-P4/`.
The failed folder/log is retained but is not a completed review. All three
figures were inspected, with the corrected `Fraction expressing` legend (0–1).
The pre-existing rotated `Features` annotation remains; no labels are clipped.

On all original cells, raw totals remain 84,811,621 molecules, SoupX retains
82,183,607 (3.09865% removed), and new CellBender retains 81,056,997 (4.42702%
removed). Corrected-total correlation is 0.999791369724995 and removed-total
correlation is −0.564400184270654. On the 13,421 jointly called cells, CellBender
removal is 4.42274%, corrected-total correlation 0.9998011364546961 and
removed-total correlation −0.5742833937479039: the qualitative conclusion is
unchanged. The new run fixes the previous convergence concern, but high
corrected-total agreement still does not establish equivalent corrections.

Across the 26 provisional markers, weighted positive-reference retention is
99.80% for SoupX and 99.72% for new CellBender; negative-reference removal is
35.40% and 40.71%, respectively. Excluding four fitting markers gives 99.81%
and 99.70% positive retention. More removal is not evidence of better correction;
marker references remain post-hoc and six clusters unresolved.

Independent checks confirmed cell/feature IDs, cell-call flags, unchanged
Raw/SoupX counts/clusters/marker statistics, gene/cell totals, marker retention,
detection changes, both correlations, and source script/H5 hashes. The jointly
called totals agree with CellBender's own validation. Source SoupX QC metadata
links this review while preserving the previous comparison. No count matrices,
CellBender files or downstream checkpoints changed; no other samples submitted.
Scripts and documentation changes remain local/uncommitted.

### All-sample automatic baselines, parallel array

The user authorized SoupX for all 23 samples and parallel processing. All
unfiltered/filtered RSEC archives and mRNA cell-count metrics were checked for
presence; counts range from 10,417 to 17,082 original calls. The fixed manifest
is `scripts/soupx/samples.txt`. Methods and background gates match the pilot;
no manual rho, sensitivity runs, cell filtering, doublet calls, MT-prefix
renaming or downstream checkpoint changes are included.

The launcher now supports array dispatch. Each task requests 4 CPUs, 16 GiB,
4 hours; at most four run concurrently. Array jobs share a run ID while sample
directories remain separate. Task IDs, manifest and script hashes are recorded.
Shared `scripts/soupx/plots.R` produces Raw/SoupX marker plots and per-cluster
marker count/detection/log-expression tables for each successful baseline.
It reuses the pilot's 26 candidate genes, not its cell-type/cluster assignments.
Missing genes are recorded and no automatic biological acceptance is performed.

| Job | Last verified result |
|---|---|
| 44409538 | Plotting preflight FAILED 1:0: Assay5 constructor is exported by SeuratObject, not Seurat; corrected namespace |
| 44409543 | Plotting preflight COMPLETED 0:0, 50s; exported marker counts matched independent original cluster totals |
| 44409555 | Initial 23-task array stopped: first four tasks FAILED 1:0 at plotting because an empty annotation label created a zero-length variable name; remaining tasks cancelled, outputs retained |
| 44409564 | Dependent verification cancelled with the stopped initial array |
| 44409746 | Exact final plotting helper preflight COMPLETED 0:0, 53s; both rendered figures inspected after restoring the working nonempty Features label |
| 44409807 | Replacement array `--array=1-23%4`: all 23 tasks COMPLETED 0:0; all samples numerically validated |
| 44409813 | Independent saved-count and marker arithmetic verification COMPLETED 0:0, 3m15s; all 23 samples passed |
| 44410208 | Plot-only array: all 23 tasks COMPLETED 0:0; widened plots, with identical numerical marker summaries |

Current data: `data/soupx/<sample>/job-44409807/`.
Current diagnostics: `results/soupx/job-44409807/<sample>/`.
Logs: `logs/soupx-44409807_<task>.log`; verification log:
`logs/soupx-verify-44409813.log`; rerender logs:
`logs/soupx-plots-44410208_<task>.log`. Independent verification reopened every
saved matrix and checked integer/nonnegative/no-added counts, marker sums,
detection fractions and cell-wise log1p normalization against exported tables.
The temporary verification/rerender scripts were removed after execution.

All 299,393 original called cells were retained across 23 libraries. Automatic
rho estimates range from 2.3% to 4.8%; all default background gates passed,
with called-cell-rank/background-upper-bound margins of 11.05–23.11.
Each sample has all 26 candidate markers and its own 21–27 preliminary clusters.
REVO30-P4 corrected counts are byte-identical to pilot job-44404963
(MD5 `fe426d2c611f7500c18ac993daa98699`), reproducing that run rather than tuning
toward a different reference removal percentage.

All 23 native contamination diagnostic plots were inspected in labeled review
contact sheets: no major secondary peaks, boundary maxima or curves tracking
the broad prior were apparent; minor shoulders/right tails do not constitute
biological acceptance, and the half-maximum spans are not confidence intervals.
Final Raw/SoupX marker pairs were inspected for representative 21-, 23- and
27-cluster samples. Wider export dimensions separate large horizontally
adjacent dots without changing expression/detection data or comparison scales.
Each `plot_rerender.json` records the new helper hash and unchanged marker
summary; original execution `run.json` hashes remain unchanged.

Every sample remains `validated_pending_qc`: numerical completion is not
biological acceptance. Sample-specific marker preservation review, downstream
cell QC and the proposed sensitivity analyses remain pending. No CellBender
files, raw data or downstream checkpoints were modified. The requested
ten-minute progress schedule was cleared after all relevant jobs finished.
Plot-helper and run-history edits remain local/uncommitted.

### SoupX input route for the existing QC pipeline

At the user's request, ten superseded `data/soupx` run folders were removed:
eight `job-44409555` folders and the REVO30-P4 pilot folders `job-44404957`
and `job-44404963`. Only `job-44409807` remains for all 23 samples; historical
`results/soupx` diagnostics were retained. The latest corrected RDS files were
not changed.

`R/qc/io.R` now loads corrected SoupX counts through `build_soupx_seurat_obj`,
sharing the original mitochondrial/barcode renaming and sample/pig/pressure/time
metadata construction with the raw-MEX loader. No job/run-ID cell metadata is
added. `R/qc/main_soupx.R` selects `job-44409807` and otherwise preserves
`main.R`'s QC, doublet, merging, plotting and checkpoint workflow.

After removing the initially added run-ID column, verification job 44411222
COMPLETED 0:0 in 3m37s. All 23 full-sample objects retained exact corrected
counts, cell IDs, mitochondrial names/percentages and expected metadata.
A 69-cell subset merge spanning every sample preserved layers and metadata;
the original raw-MEX loader regression, missing-input rejection, unchanged
source hashes and entry-point equivalence checks passed. Oracle source review
found no defects. No downstream QC was launched; `scripts/sbatch_qc.sh` still
targets the unchanged raw-data entry point. Changes remain local/uncommitted.

### Upper-only ribosomal QC: final multiplier 3 MAD

The user replaced the lower 2.5% exclusion with an upper-only per-sample rule,
initially 4 MAD and then 3 MAD. The final threshold is median + 3 × R's default
scaled MAD, computed from all input ribosomal percentages including 0%, without
log transformation, floor or cap. Only values strictly above are flagged;
equality and low percentages pass this criterion. Feature, complexity and
mitochondrial rules are unchanged. Cutoffs are stored as `ribo_max` per sample;
plots show per-sample upper segments and high-ribosomal UMAP labels.

Verification 44414023 COMPLETED 0:0 in 3m37s. All 23 SoupX samples passed
independent cutoff/flag calculations and non-ribosomal QC regression checks.
Thresholds span 18.35848–26.80025%; 1,401 of 299,393 cells (0.46795%) are flagged
by the ribosomal criterion, versus 222 (0.07415%) at 4 MAD in verification
44413950. These are criterion flags, not necessarily additional QC exclusions.
Software fixtures covered zeros, scaled MAD, strict equality and absence of
caps/floors. Initial test 44413862 failed solely on comparing Seurat's named
logical vector with an unnamed expected vector; diagnostic 44413925 confirmed
matching values, and the expectation was corrected without changing the method.

Both filtering branches passed exact retained-cell/count checks on REVO26-C and
REVO27-N4. Historical doublet labels and saved UMAP coordinates were reused only
for software/plot tests; no inference or embedding was rerun. Final before/after
and UMAP artifacts in `.amp/in/artifacts/qc-ribo-mad3/` were inspected, and plotted
threshold segments were numerically verified. Temporary verification code and
intermediate 4-MAD figures were removed. No full QC pipeline, checkpoint writes
or source-data modifications were performed. Changes remain local/uncommitted.

### Sequential QC on SoupX input, Scrublet only

Final requested layout keeps both existing input routes: `R/qc/main.R` reads
raw filtered MEX matrices and is launched by `scripts/sbatch_qc.sh`;
`R/qc/main_soupx.R` loads the corrected SoupX run `job-44409807` into Seurat
objects and is launched by `scripts/sbatch_qc_with_soupx.sh`. Both entry points
process samples sequentially and use Scrublet only. The proposed parallel array,
separate collector and their launchers were withdrawn. `main.R` was briefly
removed in error, then restored with the requested Scrublet-only changes.
`scripts/sbatch_soupx.sh` remains the unchanged ambient-correction launcher.

DoubletFinder calls, its unused clustering/sweep workflow, clean branch and
diagnostic panel are removed from QC; the BD interpolation helper now lives in
`scrublet.R`. The existing BD table, including 20,000 cells at 4.7%, is retained.
Scrublet uses original input counts/all cells and its existing score >0.15 call;
temporary normalized PCA/UMAP supplies diagnostic coordinates only. The accepted
3-MAD upper-only ribosomal rule and other QC settings remain unchanged.

Pilot array 44414192 completed two samples. Collector 44414252 and dependent
check 44414254 were cancelled when parallelization was withdrawn. Cancellation
interrupted the test clean-RDS write; check 44414264 failed to read that incomplete
file. Those disposable pilot outputs were removed, not accepted as valid QC.

Replacement sequential verification 44414269 COMPLETED 0:0 in 4m39s, peak RSS
about 8.2 GiB. It executed the actual `main_soupx.R` with only the sample list
restricted to REVO26-C/REVO30-P4. Independent rereading confirmed exact input
counts, mitochondrial/QC percentages, metadata/cell IDs, Scrublet calls,
3-MAD flags, retained-cell sets and per-sample >=10-cell gene filtering, including
genes detected in exactly 9 versus 10 retained cells. No DoubletFinder namespace
was loaded. All 16 PNG/SVG files were produced; the Scrublet UMAP was inspected.
REVO26-C retained 10,294/11,604 cells (88.71079%); REVO30-P4 retained
11,945/13,422 (88.99568%). Review plots and the numerical check table are in
`.amp/in/artifacts/qc-scrublet-test/`; disposable test RDS files were removed.

Full sequential QC job **44415111** was submitted with
`sbatch scripts/sbatch_qc_with_soupx.sh` (6 CPUs, 128 GiB, 5-hour limit) and
COMPLETED 0:0, Oct 7 17:17:43–17:59:14 CEST. Read-only output-check job
**44415118** FAILED 1:0 at the comparison of saved parameters with the mutable
current `params.R`, after the intentional feature/mitochondrial MAD change.
That check did not demonstrate a saved-data defect or complete the remaining
validation. Preserved outputs are
`data/qc_labeled_data/job-44415111/labeled_concatenated.rds`,
`data/clean_concatenated_data/job-44415111/clean_concatenated_scrublet.rds`
and `results/qc/job-44415111/`. Logs are
`logs/*-qc-with-soupx-44415111.log` and
`logs/qc-soupx-output-check-44415118.log`. Biological acceptance remains open.
No raw/SoupX source
counts or unrelated CellBender jobs were changed. Code changes remain local,
uncommitted and unpushed; no downstream normalization/integration was launched.

### SoupX QC rerun with feature/mitochondrial/ribosomal multipliers all 3

QC job **44415373** was submitted Oct 7 at 18:05:49 CEST through the existing
`scripts/sbatch_qc_with_soupx.sh` and COMPLETED 0:0, 18:06:16–18:47:12
(40m56s, peak batch RSS 36,131,152 KiB). It processed all
23 `data/soupx/<sample>/job-44409807/corrected_counts.rds` inputs sequentially,
with Scrublet only (6 CPUs, 128 GiB, 5-hour limit). No ambient correction or
downstream normalization/integration was rerun. Independent validation job
**44415439** COMPLETED 0:0, 18:47:19–18:51:08 (3m49s, peak batch RSS
20,575,340 KiB); computational checks passed, not biological acceptance.

Captured settings: `nFeature_RNA > 200` and
`<= min(4500, median + 3×MAD)` per sample; complexity >0.8;
mitochondrial percentage `<= min(20, median + 3×MAD)` using the sample's
feature/complexity-passing reference cells; ribosomal percentage
`<= median + 3×MAD` from all sample input cells, including zeros, without
cap/floor/lower cutoff. MAD is R's default scaled MAD (constant 1.4826).
Genes require detection in at least 10 retained cells per sample.
Scrublet receives all input integer counts, uses the current interpolated BD
multiplet-rate table, min_counts=3 and 30 PCs, and applies scores >0.15;
its Python random_state is **0**, not the diagnostic seed. Temporary diagnostic
LogNormalize/PCA/UMAP uses 2,000 variable features, PCs 1–20 and seed **1234**;
saved checkpoints contain counts, not those temporary normalized assays.

Validation compared the saved execution parameters with the pre-submission
settings, not live `params.R`, and decoded **eight PNGs**,
not the old 16 PNG/SVG files. All 23 exact saved input count matrices,
IDs/metadata, independently calculated QC metrics/flags/thresholds, Scrublet-only
calls/retention and sample-specific gene filtering passed, including the
9/10-cell boundary. All input SHA256 hashes and captured QC source hashes
matched at completion review. Logs:
`logs/20261007-180617-qc-with-soupx-44415373.log`,
`logs/qc-soupx-verify-44415439.log`. Actual parameters and sample thresholds
remain embedded in checkpoint `misc$qc` and documented in `.INFO`.

The [labeled checkpoint](data/qc_labeled_data/job-44415373/labeled_concatenated.rds)
retains all **299,393** input cells. **252,687** pass QC; Scrublet calls 6,115
doublets among all input cells and removes 5,015 additional QC-passing cells.
The [clean checkpoint](data/clean_concatenated_data/job-44415373/clean_concatenated_scrublet.rds)
retains **247,672 cells (82.724713%)** across all 23 samples. Per-sample retention
is 76.337665–85.953120%; sample layers retain 17,578–18,860 genes.
The original validation log records every sample's counts and thresholds.
Feature ceilings range from
3588.8547 to 4500; mitochondrial ceilings 12.348004–18.918261%; ribosomal ceilings
18.358476–26.800253%. The ribosomal criterion flags 1,401 cells, not necessarily
1,401 unique QC exclusions.

All eight final diagnostic PNGs were inspected after 19:23 CEST on Oct 7:
all sample labels/panels and legends are readable, without clipping or missing
panels. Before/after plots show QC-only retention, not doublet removal; sample
UMAPs are independent, not a joint embedding or biological annotation. Dense
point overplotting obscures violin interiors, and shared axes leave empty space
in after panels; figures were not changed. The log contains the existing
R-native UWOT default-method notice, import-replacement warning and plotting
scale-replacement messages, with no execution error or out-of-range BD clamp.

Accurate `.INFO` files now accompany the labeled-data, clean-data and results
directories. At the user's request, the extra provenance folder, duplicate
logs/source snapshots and validation CSV/RDS were removed; the results folder
contains only the eight diagnostic PNGs and `.INFO`. Original logs and data
checkpoints are unchanged. Future runs start `.INFO` with their inputs/settings
and update it with the final status/checks, without extra provenance folders.
Ten-minute monitoring ended after verified completion; the parent
thread was notified. Old run 44415111 and unrelated/CellBender edits/jobs remain
untouched. No QC code was changed for this run; records and outputs remain local,
uncommitted and unpushed. Biological acceptance/annotation remains open.

### LogNormalize on SoupX/QC 44415373

Job **44416449** was submitted through the unchanged
`scripts/sbatch_norm_feat.sh` with
`data/clean_concatenated_data/job-44415373/clean_concatenated_scrublet.rds`.
It ran Oct 7, **19:48:39–20:27:41 CEST** (39m02s, peak RSS about 44.6 GiB) and
**COMPLETED 0:0**. Preflight independently reopened the counts-only input:
**247,672 cells**, **21,084** RNA features, all **23** samples and separate
count layers, aligned cell metadata and finite regression covariates.

Configured method: LogNormalize (10,000), per-layer VST/consensus 3,000 HVGs,
ScaleData regression of `percent.mt` and `nFeature_RNA`, 50-PC PCA, seed 1234.
Resources: 6 CPUs, 256 GiB, 5-hour limit; Seurat 5.5.1/SeuratObject 5.4.0.
Checkpoint: `data/norm_feat/job-44416449/lognorm.rds` (7.6 GB);
diagnostics: `results/norm_feat/lognorm/job-44416449/` (30 PNGs).
Log: `logs/20261007-194839-norm-feat-pipeline-44416449.log`; `.INFO` files
record inputs/settings and verification.

Independent verification passed: cell/gene IDs, metadata and `misc$qc`
unchanged; all 23 counts layers exactly identical to the input; 23 per-sample
`data.*` layers plus one joined `scale.data` (3,000 × 247,672, Seurat's
standard regression output); saved settings match the configured run;
LogNormalize values reproduced with 0 difference on 3 samples × 400 cells; an
independent `lm` reproduction of regression/standardization/±10 clipping
matched stored scaled values to ~1e-15, with covariate correlations ≈0
(mean 0.002); 3,000 HVGs per sample and 3,000 consensus (each selected in ≥8
samples, 1,064 in all 23); PCA loadings use all 3,000 HVGs and 50 finite PCs
align with cells (PC1–30 explain 91.7%). All 30 diagnostic PNGs were
inspected; only cosmetic notes (label/leader overlaps on dense VST plots;
the blue-dot `VizDimLoadings` style matches prior verified run 44281196).
The joined `scale.data` means/sds deviate slightly from 0/1 solely because of
the regression and Seurat's default ±10 clipping (980,637 entries at +10).

Computational verification only; biological acceptance remains open. No
method/code changes, integration, clustering, commits or pushes; unrelated
jobs are untouched.

### Integration plotting: single-checkpoint UMAPs (CCA working route)

Per user decision, CCA is the working integration route; cross-method
comparison figures are withdrawn. `plot_integration_umaps()` now plots the
checkpoint's declared method (`misc$integration$method`) beside its
unintegrated PCA UMAP and stops, naming the reduction, if the method UMAP is
missing; layouts/dimensions scale to two panels. `R/integration/plot_main.R`
takes `<integration_checkpoint.rds> <new_output_dir>`, validates the
integration metadata, renders the existing groupings (sample, pressure, time,
pressure_time plus per-method sample facets; 6 PNGs) and writes a concise
`.INFO` per the run-record rule. `scripts/sbatch_integration_plots.sh` uses
the same two arguments. `read_integration_comparison()` (including its stale
job-44281196 pin) is removed from `R/integration/io.R`; `comparison_sources.csv`
is no longer written. `main.R` per-method runs, `cca.R`/`harmony.R`/`scvi.R`,
`plot_integration_metrics()` (per-run iLISI/silhouette violin) and
`R/clustering` comparisons are unchanged; Harmony/scVI stay runnable for
reviewers. README documents CCA-first usage.

Plot-only verification on historical `data/integration/cca/job-44282203/lognorm.rds`
(no integration computed): job 44417083 COMPLETED 0:0 in ~9 min, producing all
six expected PNGs (2-panel Unintegrated|CCA sample figure; pressure, time,
pressure_time; unintegrated and CCA sample facets) plus `.INFO`. Job 44417228
COMPLETED 0:0: existing output directory rejected before reading the 7 GB
checkpoint, and a checkpoint stripped of `umap_cca` fails with the named
reduction. Figures inspected: sample figure shows exactly two panels with the
23-sample legend; faceted CCA shows 23 labeled panels. Review copies:
`.amp/in/artifacts/integration-plot-cleanup/figures/`. A test expectation was
corrected after the first negative run (the check greps for the named
reduction; the production error message was improved accordingly). Temporary
test scripts were removed. No integration job was launched; changes remain
local/uncommitted.

### Fresh Harmony on SoupX/LogNormalize 44416449

Job **44435612** was submitted Oct 8 through the unchanged
`scripts/sbatch_integration.sh harmony data/norm_feat/job-44416449/lognorm.rds`
and **COMPLETED 0:0**, **18:47:51–19:23:20 CEST** on n038 (35m29s;
peak batch RSS 47,537,944 KiB). Existing read-only validation 44436241 passed.
Input was previously verified as **247,672 cells**, 21,084 RNA features and
23 samples. Final RDS is 7.6 GB; diagnostic CSV/RDS/PNG files exist.

Current execution settings: direct Harmony sample correction on PCA axes 1–20,
maximum 10 iterations, no projected loadings, seed 1234; PCA/Harmony UMAPs;
50,000 proportionally stratified diagnostic cells, iLISI perplexity 30 and
Euclidean sample-label silhouettes. Resources: 8 CPUs, 128 GiB, 24-hour limit.
Fresh directories: `data/integration/harmony/job-44435612/` and
`results/integration/harmony/job-44435612/`; both started concise `.INFO`
records using the `# Local provenance` template. Original log:
`logs/20261008-184751-integration-harmony-44435612.log`.

Previous Harmony **44435481** was explicitly user-cancelled and its run folders
removed before this submission. CCA **44417322** remained running and untouched.
Monitoring was changed from two hours to five minutes at the user's request;
terminal settings/data/
embedding/diagnostic checks and figure inspection will be recorded here.
No method/code changes, provenance folders, snapshots, duplicate logs, extra
validation bundles, commits or pushes were made for this submission.

Harmony 2.0.5 reached **10/10 iterations without its explicit convergence
message**: the outer convergence criterion did not pass. Both UMAP fits
finished. The user approved retaining this baseline to completion, then a
fresh **30-iteration maximum** run from the same normalization checkpoint.
Only `max_iter` in `R/integration/harmony.R` was changed (10 → 30) afterward;
the running baseline had already loaded and executed the 10-iteration call.
Other settings, convergence tolerance, CCA and unrelated work are unchanged.
R parsing and direct call-argument checks passed; rerun submission is pending
in this initial approval record; subsequent submission/cancellation is below.
No commits or pushes.

The user subsequently requested skipping integration diagnostics in future
runs while retaining their code. `R/integration/main.R` no longer sources
diagnostic/plot code, selects diagnostic cells, computes/exports mixing metrics,
or records unused diagnostic settings. Integration, both UMAP embeddings and
the single final Seurat RDS are retained. `diagnostics.R` and metric plotting
in `plots.R` remain available separately. This shared-driver change applies to
future Harmony/CCA/scVI runs; running baseline 44435612 and CCA 44417322 already
loaded their previous driver and remain untouched. Mocked driver tests passed
for all three routes, checking method/UMAP calls, actual RDS save/read, retained
provenance and absence of diagnostic calls/outputs; real-data validation of
the forthcoming run was pending at that check.

### Cancelled 30-iteration Harmony rerun 44436227

Fresh job **44436227** used the same input and launcher; started Oct 8 at
**19:24:54 CEST** on n006 (8 CPUs, 128 GiB, 24-hour limit). Only the approved
iteration ceiling (30) and omitted diagnostic orchestration differed from the
baseline. Input checks passed (247,672 cells / 23 samples), then the job was
**CANCELLED at the user's request at 19:26:42 CEST**, during Harmony 1/30
(1m48s; batch exit 0:15). Parent issued `scancel`. No final RDS or UMAP
embedding was saved; no completed convergence or output validation is available.

Both `job-44436227` run folders initially received final `# Local provenance`
`.INFO` records, then were removed at the user's explicit request on Oct 8:
`data/integration/harmony/job-44436227/` and
`results/integration/harmony/job-44436227/`. Original logs were retained. Log:
`logs/20261008-192454-integration-harmony-44436227.log`. Monitoring schedule
was cleared; **do not resubmit or start plotting jobs**. Completed baseline
44435612 and running CCA 44417322 remain untouched. No commits or pushes.

Baseline's produced metric PNG was inspected: both panels/axes/categories
are readable and unclipped; mixing distributions overlap and do not establish
biological preservation. The baseline driver saved both UMAP embeddings but
does not render UMAP figures; `R/integration/plot_main.R` is the separate
plotting entry point and was not launched here. Read-only baseline validation
job **44436241** was submitted before the stop request, with temporary code
under `.amp/in/` and original log `logs/harmony-verify-44436241.log`. It had
already **COMPLETED 0:0 at 19:30:33 CEST** when the user asked to cancel it,
so no active validation job remained to cancel. Its log reports exact preserved
IDs/metadata/RNA assay/PCA/upstream provenance, aligned finite embeddings and
saved settings, diagnostic score/summary agreement. No analysis objects or
bundles were written; temporary verification code was removed. Baseline files
were not edited. Harmony's iteration ceiling/objective history are not embedded
in misc$integration; the original execution log and .INFO record the ceiling.

**User direction:** stop submitting resource-consuming extra validation/check
jobs for successful runs unless an actual failure requires investigation or
the user explicitly requests them. Main thread informed; no replacement,
further checks or plotting jobs launched here. Monitoring remains cleared.

### UMAP plotting within each integration job

At the user's request, future integration runs now call
`plot_integration_figures()` from `R/integration/plots.R` after saving their
final RDS and before releasing the in-memory object. The plotting function
owns grouping/faceting, filenames and the six-PNG output check; `main.R`
only calls it and records saving/plotting status in concise `.INFO` files.
Plotting failures propagate while retaining the saved object. Existing plot
styles and unintegrated-plus-selected-method panels are unchanged; diagnostic
sampling/mixing metrics remain disabled.

Removed the obsolete `R/integration/plot_main.R`, its dedicated batch launcher
and `R/integration/README.md`; root README describes the combined workflow.
Local syntax and mocked control-flow checks passed for Harmony/CCA/scVI and
the plotting-failure path; no actual figures were rendered in these checks.
No cluster jobs, existing output changes, commits or pushes. Running CCA and
the completed baseline retain their already-loaded code and existing files.

### Requested Harmony series: 10, 30 and 50 iterations

The user explicitly requested removal of both `job-44435612` folders from
`data/integration/harmony/` and `results/integration/harmony/`; removed Oct 8,
with original logs retained. Historical completion/validation above describes
that removed baseline, not an available checkpoint.

The new requested series runs sequentially on the same
`data/norm_feat/job-44416449/lognorm.rds`, using the existing Harmony launcher
and unchanged sample correction, axes 1–20, seed 1234 and convergence tolerance.
`R/integration/params.R` now owns `harmony_params$max_iter`; the wrapper uses
that setting. Actual loaded parameters and log-derived convergence are saved
in `.INFO` and `misc$integration`. Convergence at the final allowed iteration
is distinguished from exhaustion without convergence. Local control-flow
checks passed; no separate validation/check/plot jobs.

| Iteration ceiling | Job | Status / original log |
|---|---|---|
| 10 | 44436475 | CANCELLED Oct 8 20:06:32 CEST after 13m37s; did not converge within 10 iterations; interrupted Harmony UMAP before final saving; `logs/20261008-195255-integration-harmony-44436475.log` |
| 10 (fresh restart) | 44437024 | COMPLETED 0:0; Oct 8 20:12:57–20:47:30 CEST, 34m33s; did not converge within 10 iterations; final RDS and six PNGs saved; `logs/20261008-201257-integration-harmony-44437024.log` |
| 30 | 44437561 | COMPLETED 0:0; Oct 8 20:50:59–21:23:59 CEST, 33m00s; converged after 11 iterations; final RDS and six PNGs saved; `logs/20261008-205059-integration-harmony-44437561.log` |
| 50 | Not submitted | Skipped at Faouzi's request after convergence at 11 iterations in the 30-ceiling run |

Faouzi explicitly approved resuming the Harmony-only series with fresh job
44437024 after the stop. After its successful completion, only
`harmony_params$max_iter` changed from 10 to 30 for job 44437561, using the
same launcher/input. No additional code/method changes or validation jobs.

All six 10-iteration PNGs were inspected (sample facets downscaled): no blank
labeled panels or visible clipping. Dense points, similar colors and low-contrast
facets limit assessment. The sample UMAP has localized color enrichment despite
broad overlap; neither numerical convergence nor correction quality follows
from this appearance. Both `.INFO` records include terminal status and limits.
No object reopening or independent checks of saved settings, cells/features,
metadata/count layers or embedding finiteness/alignment, per the resource constraint.

The 30-ceiling run completed successfully, with convergence after 11 iterations
confirmed in its log and both `.INFO` records. All six PNGs inspected (facets
downscaled): populated labeled panels without visible clipping; localized sample
color enrichment remains despite broad overlap. Dense points/similar colors and
low-contrast facets limit interpretation; numerical convergence is not evidence
of uniformly good mixing or correction quality. Both `.INFO` records finalized
with accounting and the same independent-object-check limitations above.

Each job renders the six existing UMAP PNGs in the same process after saving
the final object; diagnostics remain off. Five-minute monitoring ended after
the 30-ceiling run; no 50-iteration job submitted. Params remain at 30.
Convergence and actual outcomes are recorded per job. This resumed series
must not touch jobs owned by other threads, even if a request names them;
refer such requests to the owning thread. No commits or pushes.

## CCA job 44417322: mistaken cancellation, Oct 8

- Input: `data/norm_feat/job-44416449/lognorm.rds` (247,672 cells,
  23 samples); CCA, axes 1–20, seed 1234; original diagnostic-enabled driver.
- Terminal status: **CANCELLED**, Oct 7 21:37:44–Oct 8 20:05:53 CEST
  (22h28m09s). `sacct` records CANCELLED by 18174; batch exit 0:15.
- Cause: [Harmony thread](https://ampcode.com/threads/T-01a11c69-3a3d-7238-baf1-c845dda65c5e)
  mistakenly issued `scancel 44417322`, despite the leave-CCA-untouched
  constraint. Faouzi clarified that Harmony was intended and CCA must not
  be stopped. This was **not an explicit user request to cancel CCA**.
  The inaccurate data `.INFO` attribution has been corrected; both run
  `.INFO` records now reflect the terminal status and actual cause.
- Last log evidence: CCA anchor finding, ending with 8,696 anchors and
  `Running CCA`. No final checkpoint, UMAPs, diagnostics or metrics PNG
  saved. Data directory contains `.INFO`; results contain `.INFO` and
  `diagnostic_cells.csv` only. Existing folders, partial outputs and logs
  are preserved. Earlier references to CCA remaining running/untouched
  above describe the state before this cancellation.
- Original log: `logs/20261007-213744-integration-cca-44417322.log`.
- Old CCA monitoring cleared. No recovery or check jobs launched. Restart
  was withheld until Faouzi gave fresh explicit approval; see new run below.
  Changes remain local/uncommitted; unrelated edits and jobs were not changed.

## CCA job 44436792: approved restart, Oct 8

- Faouzi explicitly requested a rerun after correction of the cancellation
  records. Submitted the unchanged launcher with input
  `data/norm_feat/job-44416449/lognorm.rds` (247,672 cells, 23 samples).
- RUNNING on n027; Slurm start Oct 8 20:09:26 CEST; driver `.INFO` start
  20:09:40 CEST. CCA, axes 1–20, seed 1234; 8 CPUs, 512 GiB, 72-hour limit.
- Current local driver writes its own `.INFO`, integrates and saves PCA/CCA
  UMAPs and renders six UMAP PNGs. Diagnostic sampling/mixing metrics are
  disabled; no integration code or parameters changed for this submission.
- Checkpoint pending: `data/integration/cca/job-44436792/lognorm.rds`;
  results: `results/integration/cca/job-44436792/`.
- Original log: `logs/20261008-200926-integration-cca-44436792.log`.
- Two-hour monitoring resumed for this job; cancelled job 44417322's folders,
  partial outputs and logs preserved. No separate recovery/check/plot jobs.

## QC sensitivity comparison 44437522: Oct 8

- User-approved comparison of mitochondrial cap 15%, ribosomal +2 MAD, and both,
  against SoupX QC 44415373 (299,393 input cells, 247,672 baseline clean cells).
- Reused the original labeled checkpoint's QC metrics, saved execution settings,
  Scrublet calls and sample UMAP coordinates. No ambient correction, doublet
  detection, embedding, normalization or integration rerun; no clean RDS replacement.
- Job **44437522 COMPLETED 0:0**, 20:40:58–20:46:35 CEST (5m37s).
  Independent metadata/table check and UMAP legend refinement **44437549
  COMPLETED 0:0**, 20:48:29–20:49:46 (1m17s).

| Scenario | Retained cells | Input retained | Additional baseline-clean exclusions |
|---|---:|---:|---:|
| Baseline | 247,672 | 82.72% | 0 |
| Mitochondrial cap 15% only | 245,599 | 82.03% | 2,073 |
| Ribosomal +2 MAD only | 239,426 | 79.97% | 8,246 |
| Both | 237,410 | 79.30% | 10,262 |

- Overlap: 57 cells. All 23 samples remain. Both changes remove 8.94% of
  baseline-clean REVO29-N4 cells, the largest sample-relative reduction.
- [Outputs](results/qc/sensitivity-job-44437522/): concise `.INFO`,
  `per_sample_retention.csv`, and six inspected PNGs per alternative (five native
  before/after QC plots and one newly-excluded-cell UMAP). Baseline plots unchanged.
- Baseline reproduced exactly; all 92 sample/scenario rows independently checked;
  cell sets nested, combined set equals intersection, source checkpoint hashes,
  Scrublet calls and UMAP coordinates unchanged. All 18 final PNGs decoded and
  inspected; unused UMAP legend categories removed without changing counts.
- Original logs: `logs/qc-sensitivity-44437522.log` and
  `logs/qc-sensitivity-plotcheck-44437549.log`. No provenance folders, copied logs
  or extra validation bundles. Temporary run code removed after completion.
- Computational sensitivity analysis, not proof of biological improvement.
  Baseline 44415373, live QC parameters and unrelated/CellBender work unchanged.

## Independent-threshold QC sensitivity 44437686: Oct 8

- Faouzi approved estimating each QC threshold separately on all input cells
  within each sample, then combining flags before removal. The only historical
  dependency was mitochondrial MAD estimation on gene/complexity-passing cells.
- Target corrected by Faouzi to SoupX QC **44415373**, not CellBender 44436358.
  Saved metrics/settings, Scrublet calls and independent sample UMAPs reused.
- Job **44437686 COMPLETED 0:0**, 21:13:36–21:18:00 CEST (4m24s).

| Scenario | Retained cells | Input retained |
|---|---:|---:|
| Historical baseline | 247,672 | 82.72% |
| All-input design, unchanged baseline settings | 247,282 | 82.59% |
| Previous combined sensitivity (mito15/ribo2) | 237,410 | 79.30% |
| All-input design, new combined settings | 230,817 | 77.10% |

- New combined: genes >=300 and <=min(4000, median+3 MAD), mitochondrial
  <=min(15%, median+3 MAD), ribosomal <=median+2 MAD; complexity >0.8.
- Design alone: 413 newly excluded and 23 newly retained, net loss 390
  (0.16% of baseline). Mitochondrial ceilings decreased in 19 samples and
  increased in four; retention sets are therefore not necessarily nested.
- New combined versus previous combined: 6,610 newly excluded, 17 newly
  retained, net loss 6,593 (2.78%). All 23 samples remain. Largest relative
  net loss versus previous combined: REVO29-N4, 738/8,068 cells (9.15%).
- Historical flags/thresholds/retention and previous combined per-sample results
  reproduced. Sparse counts agree with saved gene/UMI metrics; IDs aligned;
  alternate MAD derivations and gene-threshold independence checks passed.
  All 92 exported rows reloaded and gain/loss identities independently checked;
  source checksum, Scrublet calls and UMAP coordinates unchanged.
- [Outputs](results/qc/sensitivity-job-44437686/): concise `.INFO`,
  `per_sample_retention.csv`, and six inspected PNGs per new scenario. All 12
  PNGs decoded; labels/legends complete with no clipping. Native violin plots
  retain dense point overplotting. Original log: `logs/qc-independent-44437686.log`.
- No replacement clean RDS, live QC code changes, ambient correction or
  downstream reruns. No provenance folders, copied logs or validation bundles;
  own temporary scripts removed. Existing runs and CellBender preserved.
- Computational sensitivity only; retention changes do not establish biological
  quality. Findings reported to Faouzi and the parent thread; monitoring cleared.

## Production SoupX QC 44438067: new combined limits, Oct 8

- Faouzi authorized adopting the new combined limits, independent all-input
  labeling, and a full QC run through `sbatch scripts/sbatch_qc.sh soupx`.
- Code: `R/qc/filter.R` now estimates mitochondrial MAD on all input cells per
  sample, independently of gene/complexity flags; gene minimum is inclusive.
  `R/qc/params.R`: genes >=300, ceiling min(4000, median+3 MAD), mitochondrial
  <=min(15%, median+3 MAD), ribosomal <=median+2 MAD; complexity >0.8 unchanged.
- All input cells receive QC and Scrublet labels before any removal. Clean cells
  have no positive QC flag and no Scrublet doublet call. Sequential 23-sample loop;
  same SoupX producer 44409807, now under `data/raw_data/soupx/`.
- Production **44438067 COMPLETED 0:0**, Oct 8 23:06:54-23:43:36 CEST on n006
  (36m42s); 6 CPUs, 128 GiB, 5-hour limit; peak batch RSS 35,221,796K.
  No normalization/integration or ambient correction rerun by this QC thread.
- Labeled: `data/qc_labeled_data/job-44438067/labeled_concatenated.rds`;
  clean: `data/clean_concatenated_data/job-44438067/clean_concatenated_scrublet.rds`;
  figures: `results/qc/job-44438067/` (eight PNGs decoded and visually inspected).
- `.INFO` initialized in all three directories with actual inputs/settings/seeds.
  Finalized with outcome and validation; no provenance folders or copied logs.
  Saved objects capture actual settings, threshold reference and inclusive minimum.
  Original log: `logs/20261008-230655-qc-pipeline-44438067.log`.
- Preflight **44438064 COMPLETED 0:0**: actual REVO29-N4 thresholds/union flags,
  gene equality boundaries, independence from gene/complexity criteria, unchanged
  Scrublet scores and sensitivity agreement passed. Earlier preflight 44438061
  failed on named-versus-unnamed vector comparison; score values proved identical
  after correcting the test. Logs: `logs/qc-preflight-44438061.log` and
  `logs/qc-preflight-44438064.log`; production code was not altered for this test fix.
- Validation **44438256 COMPLETED 0:0**, Oct 9 00:00:19-00:03:59 CEST (3m40s).
  Log: `logs/qc-combined-verify-44438256.log`. Passed saved settings, independently
  derived per-sample thresholds, source/clean counts aligned by gene/cell IDs,
  metadata, union flags, Scrublet-only labels and sample-specific gene 9/10 boundary.
- Verified **299,393 input -> 235,569 QC-passing -> 230,817 Scrublet-clean (77.10%)**;
  all 23 samples retained; exact per-sample agreement with sensitivity 44437686.
  Scrublet calls 6,115 doublets overall, excluding 4,752 additional QC-passing cells.
- Earlier validator 44438069 failed on assuming historical scores must be identical.
  Diagnostic 44438245 isolated REVO30-P10: existing BD table additions at 18,000
  and 19,000 change its expected rate from 0.040191333333333336 to
  0.040163999999999998. Maximum score difference 0.00017719301723395642;
  no doublet call changed. Scrublet prior-odds mapping explains every score change
  (maximum residual 3.7747582837255322e-15); other 22 samples' scores identical.
- Validator 44438254 then failed on assuming source/merged gene order identical;
  aligning by IDs resolved it. These were test assumptions, not demonstrated data
  defects; production data/code were not altered to make validation pass.
  Original logs: `logs/qc-combined-verify-44438069.log`,
  `logs/qc-score-check-44438245.log`, `logs/qc-combined-verify-44438254.log`.
- Parent and normalization owner receive verified completion and exact clean
  checkpoint; downstream LogNormalize belongs to the normalization thread.
  Verification scratch removed and ten-minute monitoring cleared at completion.
  Existing runs, unrelated edits/jobs and CellBender preserved.
- Changes remain local/uncommitted; no commits or pushes.

## LogNormalize on verified QC 44438067: Oct 9

- Normalization **44438396** submitted through unchanged
  `scripts/sbatch_norm_feat.sh`, **COMPLETED 0:0**;
  **00:09:20–00:44:19 CEST** (34m59s), peak RSS 43150692 KiB.
- Input: `data/clean_concatenated_data/job-44438067/clean_concatenated_scrublet.rds`,
  explicitly released after QC validator 44438256 passed. Independent focused
  preflight: **230,817 cells**, **20,775 features**, 23 counts-only sample layers,
  aligned IDs/metadata and finite regression covariates.
- Unchanged LogNormalize 10,000, 3,000 per-layer VST/consensus HVGs,
  ScaleData regression of percent.mt/nFeature_RNA, 50-PC PCA, seed 1234.
  Seurat 5.5.1/SeuratObject 5.4.0; 6 CPUs, 256 GiB, 5-hour limit.
- Outputs: `data/norm_feat/job-44438396/lognorm.rds` and
  `results/norm_feat/lognorm/job-44438396/`; both finalized with `.INFO`.
  Original log: `logs/20261009-000920-norm-feat-pipeline-44438396.log`.
- Independent saved-object verification passed: exact cell/gene IDs, all
  metadata/misc$qc and all 23 original count matrices unchanged; corresponding
  finite per-sample data layers and joined 3000 × 230817 scale.data aligned.
  Saved settings match execution; 3000 HVGs/sample and consensus, all 3000
  used in 50 finite aligned PCs. PC1–30 explain 91.52816% among computed PCs.
- Independent normalized values on 40 cells/sample (920 total) match exactly;
  three-gene lm regression/scaling/upper-cap reproduction matches to 4.44e-16,
  including SFTPC's 1116 capped cells. No PCA features omitted.
- All 30 PNGs decoded and inspected (23 VST via temporary sheets, other seven
  individually): no empty/missing panels, minor VST label/leader crowding and
  REVO31-C rightmost tick truncation. Existing blue loading-point style retained.
- Verification scratch removed; computational checks do not establish biological
  acceptance. Faouzi requested two main-thread handoffs for new Harmony and CCA
  using this same verified checkpoint; integration is owned by their threads.
- Both separate handoff messages were sent to the main thread after verification;
  ten-minute normalization monitoring cleared. Integration launches remain with
  the main/Harmony/CCA threads, not this normalization owner.
- No normalization source/method changes, commits or pushes; old runs and
  unrelated jobs/edits remain untouched.

## Harmony on new normalization 44438396: Oct 9

- Faouzi authorized this new Harmony run through the normalization/main-thread
  handoff. Submitted **44438700** Oct 9 00:52:53 CEST; **COMPLETED 0:0**,
  00:53:19–01:24:54 CEST on n002 (31m35s); peak RSS 33816596 KiB.
  Log: `logs/20261009-005319-integration-harmony-44438700.log`.
- Exact new input: `/work/project/revo-pig-sc/analysis/data/norm_feat/job-44438396/lognorm.rds`,
  derived from QC 44438067; normalization owner verified 230,817 cells,
  20,775 RNA genes, 23 samples, 3,000 consensus HVGs and 50 PCs.
  This is not the older normalization 44416449.
- Unchanged launcher/settings: sample correction, axes 1–20, max_iter 30,
  seed 1234, project.dim FALSE, default convergence tolerance; 8 CPUs,
  128 GiB, 24-hour limit. No source or method changes for this submission.
- Driver created fresh `data/integration/harmony/job-44438700/` and
  `results/integration/harmony/job-44438700/` with concise `.INFO` at 00:53:31 CEST;
  both confirm the new input and 30-iteration ceiling.
  Final `lognorm.rds` and six modular UMAP PNGs saved in the same job.
- Explicit convergence after **10 iterations**, confirmed in log and both `.INFO`;
  UMAP read 230,817 cells. All six figures inspected (sample facets downscaled):
  populated labeled panels without visible clipping. Localized sample color
  enrichment remains despite broad overlap; dense points/similar colors and
  low-contrast facets limit interpretation. Convergence is not proof of good correction.
- Both `.INFO` finalized; five-minute monitoring ended. Diagnostics/metrics off;
  no extra check jobs or object reopening. Saved settings, cells/features/metadata,
  count layers and embedding finiteness/alignment were not independently revalidated.
  No unrelated jobs/edits touched, no commits or pushes.

## CCA job 44438698: new normalization input, Oct 9

- Authorized new run after norm_feat 44438396 verification; input:
  `/work/project/revo-pig-sc/analysis/data/norm_feat/job-44438396/lognorm.rds`,
  derived from QC 44438067 (230,817 cells, 20,775 RNA genes, 23 samples,
  3,000 consensus HVGs, 50 PCs). Not the older norm_feat 44416449 input.
- Submitted with existing dedicated CCA launcher; RUNNING on n019, Slurm
  start Oct 9 00:52:49 CEST; driver `.INFO` start 00:53:14 CEST.
- CCA axes 1–20, seed 1234; 8 CPUs, 512 GiB, 72-hour limit. Unmodified local
  driver writes concise `.INFO` and renders six modular UMAP PNGs in the same
  job; diagnostic sampling and mixing metrics remain disabled.
- Checkpoint pending: `data/integration/cca/job-44438698/lognorm.rds`;
  results: `results/integration/cca/job-44438698/`.
- Original log: `logs/20261009-005249-integration-cca-44438698.log`.
- Existing CCA 44436792 remains RUNNING and untouched. Three-hour monitoring
  covers both runs; no additional validation jobs or changes to Harmony,
  integration methods, parameters or unrelated work. No commits/pushes.

## Limits

- Only jobs explicitly marked completed have verified exit status. File
  presence alone does not prove a complete run.
- Older jobs did not capture immutable input hashes or per-job source revisions;
  current settings must not be applied retroactively.
- Most multi-GB checkpoints were not independently deserialized on the
  memory-limited login runner; evidence came from accounting, logs, saved
  diagnostics/bundles and file checks. scVI E is the stated high-memory exception.
- Sample-mixing metrics and exploratory clusters do not establish cell identities
  or biological preservation.
