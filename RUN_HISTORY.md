# Pipeline run history

Compact index of recorded runs through **Oct 9, 2026**; times are CEST.
This index was condensed without refreshing scheduler states or rerunning analyses.
**Completed ≠ biologically accepted.** Nonterminal states below are last recorded observations.

## Run lineage

Selected branches; tables give other parents, comparisons and retries. Arrows mean data input, not chronological succession.

```text
Raw RSEC archives → SoupX 44409807 [299,393 cells; 23 samples]
├─ QC 44461874 → LogNormalize 44462433 [214,221 cells; fixed-threshold QC branch; verified]
│  ├─ Harmony 44464751 [completed; converged at iteration 10/30]
│  │  └─ Leiden clustering 44465100 [COMPLETED; CSV-validated; resolutions 0.1–0.4]
│  │     └─ Cluster plots 44466714 [COMPLETED; both PNGs inspected]
│  ├─ Harmony 44465095 [completed; axes 1–30; converged at iteration 7/30]
│  └─ CCA 44464754 [submitted Oct 9; three-hour monitoring]
├─ QC 44438067 → LogNormalize 44438396 [230,817 cells; latest-input branch]
│  ├─ Harmony 44438700 [completed; converged at iteration 10/30]
│  │  ├─ UMAP exploration 44442925 [14/14 completed; plots only against saved Harmony baseline]
│  │  └─ Leiden clustering 44452494 [COMPLETED; validated; resolutions 0.1–0.4]
│  │     └─ Cluster plots 44456242 [COMPLETED; both PNGs inspected]
│  └─ CCA 44438698 [last recorded RUNNING, Oct 9]
├─ QC 44415373 → LogNormalize 44416449 [247,672 cells; earlier-input branch]
│  ├─ Harmony 44437024 [completed; no convergence within 10 iterations]
│  ├─ Harmony 44437561 [completed; converged at iteration 11/30]
│  └─ CCA 44436792 [COMPLETED Oct 9; verified checkpoint + 6 PNGs]
│  QC 44415373 also supplies sensitivity analyses 44437522 and 44437686.
└─ QC 44415111 [completed; independent output check incomplete]

Raw filtered MEX → QC 44276132 → Scrublet-clean → LogNormalize 44281196
                                             [226,738 cells; historical branch]
  ├─ PCA directly → unintegrated clustering D/E
  ├─ Harmony 44281709 → Harmony clustering D/E
  │  └─ Leiden clustering 44452500 [COMPLETED; validated; resolutions 0.1–0.4]
  │     └─ Cluster plots 44456247 [COMPLETED; both PNGs inspected]
  ├─ scVI 44281712 → scVI clustering D/E
  └─ CCA 44282203
     ├─ UMAP exploration 44445961 [COMPLETED 15/15; plots only against saved CCA baseline]
     └─ Leiden clustering 44452483 [COMPLETED; validated; resolutions 0.1–0.4]
        ├─ Cluster plots 44456254 [COMPLETED; both PNGs inspected]
        └─ Cluster-coloured UMAP exploration 44458934_14 + 44459373 [COMPLETED 15/15; 60 PNGs inspected; fixed labels preserved]

CCA clustering 44310241 [integration parent not recorded here]
  └─ Marker discovery 44395220 [cancelled; partial only]
```

## Leiden clustering — Oct 9, 2026

| Job / date | Stage | Parent / input | Execution | Output / checks | Record |
|---|---|---|---|---|---|
| 44452483 / Oct 9 | Leiden, CCA | `data/integration/cca/job-44282203/lognorm.rds` | COMPLETED 0:0; 47m08s; 44.91 GiB peak | Validated 226,738 preserved cells; 15/18/18/20 clusters at 0.1/0.2/0.3/0.4; no diagnostics | [Record](results/clustering/cca/job-44452483/.INFO), [validation](logs/validate-cca-44452483-44455323.log) |
| 44452494 / Oct 9 | Leiden, Harmony | `data/integration/harmony/job-44438700/lognorm.rds` | COMPLETED 0:0; 40m09s; 42.8 GiB peak | Validated 230,817 preserved cells; 13/15/17/20 clusters at 0.1/0.2/0.3/0.4; no diagnostics | [Record](results/clustering/harmony/job-44452494/.INFO), [validation](logs/validate-clustering-44452494-44455209.log) |
| 44452500 / Oct 9 | Leiden, Harmony | `data/integration/harmony/job-44281709/lognorm.rds` | COMPLETED 0:0; 48m14s; 40.8 GiB peak | Validated 226,738 preserved cells; 12/16/17/18 clusters at 0.1/0.2/0.3/0.4; validator assertion corrected, no clustering retry | [Record](results/clustering/harmony/job-44452500/.INFO), [validation](logs/validate-clustering-44452500-44455415.log) |
| 44465100 / Oct 9 | Leiden, Harmony | `data/integration/harmony/job-44464751/lognorm.rds` | COMPLETED 0:0; 40m02s; 38.7 GiB peak | CSV check: 214,221 cells, complete labels; 14/16/17/19 clusters at 0.1/0.2/0.3/0.4; no saved-object revalidation; no figures (plotting job separate, not run) | [Record](results/clustering/harmony/job-44465100/.INFO) |

### Cluster plots

| Job / date | Stage | Parent / input | Execution | Output / checks | Record |
|---|---|---|---|---|---|
| 44456254 / Oct 9 | CCA cluster UMAP + clustree | `data/clustering/cca/job-44452483/lognorm.rds` | COMPLETED 0:0; 2m19s | Both PNGs inspected, resolutions 0.1–0.4; readable labels; minor clustree node overlaps; no fitting | [Record](results/clustering/cca/job-44452483/.INFO), [log](logs/20261009-180005-clustering-plots-44456254.log) |
| 44456242 / Oct 9 | Harmony cluster UMAP + clustree | `data/clustering/harmony/job-44452494/lognorm.rds` | COMPLETED 0:0; 2m19s | Both PNGs inspected, resolutions 0.1–0.4; minor small-island label placement limitations; no fitting | [Record](results/clustering/harmony/job-44452494/.INFO), [log](logs/20261009-180005-clustering-plots-44456242.log) |
| 44456247 / Oct 9 | Harmony cluster UMAP + clustree | `data/clustering/harmony/job-44452500/lognorm.rds` | COMPLETED 0:0; 2m15s | Both PNGs inspected, resolutions 0.1–0.4; readable/unclipped labels and edges; no fitting | [Record](results/clustering/harmony/job-44452500/.INFO), [log](logs/20261009-180005-clustering-plots-44456247.log) |
| 44466714 / Oct 10 | Harmony cluster UMAP + clustree | `data/clustering/harmony/job-44465100/lognorm.rds` | COMPLETED 0:0; 2m08s; 24.6 GiB peak | Both PNGs inspected, resolutions 0.1–0.4 (14/16/17/19 clusters); readable/unclipped labels, edges and captions; no fitting | [Record](results/clustering/harmony/job-44465100/.INFO), [log](logs/20261010-000119-clustering-plots-44466714.log) |

### Cluster-coloured UMAP exploration

| Job / date | Stage | Parent / input | Execution | Output / checks | Record |
|---|---|---|---|---|---|
| 44458934_14 / Oct 9 | CCA UMAP pilot, fixed Leiden labels | `data/clustering/cca/job-44452483/lognorm.rds` | COMPLETED 0:0; 9m57s; 27.01 GiB peak | Test 14; four paired PNGs inspected; preservation checks passed; dense point overlap remains | [Record](results/clustering/cca/job-44452483/umap_explore/min.dist_1_n.neighbors_75_repulsion.strength_2/.INFO), [log](logs/20261009-183244-cluster-umap-explore-44458934_14.log) |
| 44459373 / Oct 9 | CCA UMAP exploration, remaining 14 configurations | `data/clustering/cca/job-44452483/lognorm.rds` | COMPLETED 14/14, all 0:0; 7m29s–10m39s/task; 21.43–27.02 GiB peak | 56 paired PNGs inspected; all preservation checks passed; dense overlap remains; no biological ranking | [Record](results/clustering/cca/job-44452483/umap_explore/.INFO) |

## Latest-input branch

Parents identify the actual input artifact. “Checks” describes computational evidence only; biological acceptance remains open.

| Job / date | Stage | Parent / input | Execution | Output / checks | Record |
|---|---|---|---|---|---|
| 44409807 / Oct 7 | SoupX, 23-task array | Raw unfiltered/filtered RSEC + metrics | All 23 COMPLETED 0:0 | 299,393 cells retained; rho 2.3–4.8%; counts verified; sample-specific biological QC pending | [Record](results/soupx/job-44409807/.INFO), [results](results/soupx/job-44409807/) |
| 44438067 / Oct 8 | QC, independent all-input thresholds | 44409807: corrected counts | COMPLETED 0:0 | 235,569 QC-passing → **230,817 Scrublet-clean** (77.10%); validator 44438256 passed; 8 PNGs inspected | [Record](data/clean_concatenated_data/job-44438067/.INFO), [results](results/qc/job-44438067/) |
| 44461874 / Oct 9 | QC, fixed thresholds (no MAD) | 44409807: corrected counts | COMPLETED 0:0 | Genes 300–4,000, mitochondrial ≤10%, complexity >0.8; no ribosomal filter; 218,666 QC-passing → **214,221 Scrublet-clean** (71.55%); 7 PNGs, 3 changed figures inspected | [Record](data/clean_concatenated_data/job-44461874/.INFO), [results](results/qc/job-44461874/) |
| 44462433 / Oct 9 | LogNormalize + PCA | 44461874: Scrublet-clean | COMPLETED 0:0; 21:46:32–22:19:52 CEST (33m20s) | 214,221 cells, 20,516 genes; counts/metadata preserved; normalization/scaling/PCA verified (3,000 HVGs used in 50 PCs; PC1–30 = 91.57%); 30 PNGs inspected | [Record](data/norm_feat/job-44462433/.INFO), [results](results/norm_feat/lognorm/job-44462433/) |
| 44464751 / Oct 9 | Harmony, maximum 30 iterations | 44462433: normalized checkpoint | COMPLETED 0:0; 22:34:40–23:03:11 CEST (28m31s) | Fixed-threshold branch: 214,221 cells; converged at 10; RDS + 6 inspected PNGs; **no independent saved-object revalidation** | [Record](data/integration/harmony/job-44464751/.INFO), [results](results/integration/harmony/job-44464751/) |
| 44465095 / Oct 9 | Harmony, axes 1–30, ceiling 30 | 44462433: normalized checkpoint | COMPLETED 0:0; 23:10:42–23:41:43 CEST (31m01s) | Faouzi-requested 30-PC variant; converged at 7; RDS + 6 inspected PNGs; `dims` remains 1:30 (shared with CCA/scVI); **no independent saved-object revalidation** | [Record](data/integration/harmony/job-44465095/.INFO), [results](results/integration/harmony/job-44465095/) |
| 44464754 / Oct 9 | CCA, axes 1–20 | 44462433: normalized checkpoint | RUNNING; start 22:35:10 CEST on n020 | Fixed-threshold branch (214,221 cells); same launcher/settings, seed 1234; six same-job UMAP PNGs, diagnostics off; three-hour monitoring | [Record](data/integration/cca/job-44464754/.INFO), [log](logs/20261009-223510-integration-cca-44464754.log) |
| 44438396 / Oct 9 | LogNormalize + PCA | 44438067: Scrublet-clean | COMPLETED 0:0 | 230,817 cells, 20,775 genes; counts/metadata preserved; normalization/scaling/PCA verified; 30 PNGs inspected | [Record](data/norm_feat/job-44438396/.INFO), [results](results/norm_feat/lognorm/job-44438396/) |
| 44438700 / Oct 9 | Harmony, maximum 30 iterations | 44438396: normalized checkpoint | COMPLETED 0:0 | Converged at 10; RDS + 6 inspected PNGs; **no independent saved-object revalidation** | [Record](data/integration/harmony/job-44438700/.INFO), [results](results/integration/harmony/job-44438700/) |
| 44442925 / Oct 9 | UMAP exploration, 14-task array | 44438700: Harmony checkpoint and saved default UMAP | All 14 COMPLETED 0:0; 12:12–13:13 CEST | 14 sample comparisons present; recorded settings/file times checked; no object persistence or plot inspection | [Record](results/integration/harmony/job-44438700/.INFO), [results](results/integration/harmony/job-44438700/) |
| 44438698 / Oct 9 | CCA, axes 1–20 | 44438396: normalized checkpoint | Last recorded RUNNING, Oct 9 | Final checkpoint pending in recorded state; start 00:52:49, n019 | [Record](data/integration/cca/job-44438698/.INFO), [log](logs/20261009-005249-integration-cca-44438698.log) |

### Checkpoint locations

| Stage | Path convention |
|---|---|
| SoupX | `data/raw_data/soupx/<sample>/job-44409807/corrected_counts.rds` (moved from `data/soupx/`) |
| Labeled QC | `data/qc_labeled_data/job-<ID>/labeled_concatenated.rds` |
| Clean QC | `data/clean_concatenated_data/job-<ID>/clean_concatenated_scrublet.rds` |
| Normalization | `data/norm_feat/job-<ID>/lognorm.rds` |
| Integration | `data/integration/<method>/job-<ID>/lognorm.rds` |
| Clustering | `data/clustering/<method>/job-<ID>/lognorm.rds`; settings/assignments: `results/clustering/<method>/job-<ID>/analysis.rds` |

### Execution settings that distinguish QC runs

MAD means R's scaled MAD (constant 1.4826); thresholds are per sample. For the fully specified runs below, complexity must exceed 0.8 and genes require detection in ≥10 retained cells/sample.

| QC job | Input / doublet routes | Retained gene-count range | Mitochondrial ceiling / reference cells | Ribosomal rule |
|---|---|---|---|---|
| 44194877 | Raw MEX; DoubletFinder + Scrublet | >200; logged upper cap 5,000 | Per-sample ceilings; exact MAD settings not recoverable | 2.5% labels logged, but exclusion not confirmed |
| 44276132 | Raw MEX; DoubletFinder + Scrublet | >200; ≤min(4500, median + 4 MAD) | min(20%, median + 4 MAD); feature/complexity-passing cells | Retain >2.5% |
| 44415111 | SoupX 44409807; Scrublet | >200; ≤min(4500, median + 4 MAD) | min(20%, median + 4 MAD); feature/complexity-passing cells | ≤median + 3 MAD, all input cells |
| 44415373 | SoupX 44409807; Scrublet | >200; ≤min(4500, median + 3 MAD) | min(20%, median + 3 MAD); feature/complexity-passing cells | ≤median + 3 MAD, all input cells |
| 44438067 | SoupX 44409807; Scrublet | **≥300**; ≤min(4000, median + 3 MAD) | min(15%, median + 3 MAD); **all input cells independently** | ≤median + 2 MAD, all input cells |

- Upper-only ribosomal rules include zeros in estimation and have no floor, cap or log transform; equality passes.
- SoupX QC uses corrected input counts, not uncorrected raw counts. Scrublet calls scores >0.15 on all input cells before removal.
- For 44415373, Scrublet `random_state=0`, min_counts=3, 30 PCs; diagnostic PCA/UMAP seed=1234. Earlier history reported seed 1234 for 44276132/44415111 without distinguishing these roles; do not extrapolate it.
- Normalization 44281196/44416449/44438396/44462433: LogNormalize 10,000; 3,000 VST/consensus variable genes; regression of `percent.mt` + `nFeature_RNA`; 50 PCs; seed 1234.
- Recent integrations: axes 1–20, seed 1234; Harmony corrects sample, `project.dim=FALSE`. New driver saves RDS then six UMAP PNGs in the same job; mixing diagnostics disabled. Plot failure can leave a valid saved checkpoint but a failed job.

## Earlier SoupX branch and integration retries

| Job / date | Stage / difference | Parent / input | Execution | Output / checks | Record |
|---|---|---|---|---|---|
| 44415111 / Oct 7 | QC, feature/mitochondrial 4 MAD | 44409807: corrected counts | COMPLETED 0:0 | Labeled/clean RDS retained; validator 44415118 failed on mutable-parameter comparison; remaining checks incomplete | [Log](logs/20261007-171743-qc-with-soupx-44415111.log) |
| 44415373 / Oct 7 | QC, all three multipliers 3 MAD | 44409807: corrected counts | COMPLETED 0:0 | 299,393 input → 252,687 QC-passing → **247,672 clean**; validation 44415439 passed; 8 PNGs inspected | [Record](data/clean_concatenated_data/job-44415373/.INFO) |
| 44416449 / Oct 7 | LogNormalize + PCA | 44415373: Scrublet-clean | COMPLETED 0:0 | 247,672 cells, 21,084 genes; saved counts/metadata, normalization/scaling/PCA verified; 30 PNGs inspected | [Record](data/norm_feat/job-44416449/.INFO) |
| 44417308 / Oct 7 | CCA first attempt | 44416449 | CANCELLED by user | Before integration; no checkpoint; run directories deleted | [Log](logs/20261007-212915-integration-cca-44417308.log) |
| 44417322 / Oct 7–8 | CCA, diagnostic-enabled driver | 44416449 | CANCELLED, batch 0:15 | **Mistaken cancellation by Harmony thread, not user authorization**; no checkpoint; `.INFO`/diagnostic-cell list retained | [Record](data/integration/cca/job-44417322/.INFO), [log](logs/20261007-213744-integration-cca-44417322.log) |
| 44436792 / Oct 8–9 | CCA, approved restart | 44416449 | COMPLETED 0:0; 26h14m45s; ~48.9 GiB peak | 8.1 GB [checkpoint](data/integration/cca/job-44436792/lognorm.rds) + 6 inspected UMAP PNGs; independent reopen: IDs/metadata identical, `integrated_cca` dims 1–20 finite (3,000x20 loadings), both UMAPs 2-D aligned, misc cca/seed 1234 | [Record](data/integration/cca/job-44436792/.INFO), [log](logs/20261008-200926-integration-cca-44436792.log) |
| 44435481 / Oct 8 | Previous Harmony attempt | Not recorded here | CANCELLED by user | Run directories deleted | Historical cancellation record; no exact log path recorded |
| 44435612 / Oct 8 | Harmony, maximum 10; diagnostics on | 44416449 | COMPLETED 0:0 | No explicit convergence; validator 44436241 passed; **data/results later deleted by request**, logs retained | [Log](logs/20261008-184751-integration-harmony-44435612.log) |
| 44436227 / Oct 8 | Harmony, maximum 30; diagnostics off | 44416449 | CANCELLED by user | Stopped at iteration 1; no RDS/UMAP; run directories deleted | [Log](logs/20261008-192454-integration-harmony-44436227.log) |
| 44436475 / Oct 8 | Harmony, maximum 10 | 44416449 | CANCELLED | No convergence in 10; stopped during UMAP before saving | [Log](logs/20261008-195255-integration-harmony-44436475.log) |
| 44437024 / Oct 8 | Harmony, maximum 10; fresh restart | 44416449 | COMPLETED 0:0 | No convergence in 10; RDS + 6 inspected PNGs; no independent object checks | [Record](data/integration/harmony/job-44437024/.INFO) |
| 44437561 / Oct 8 | Harmony, maximum 30 | 44416449 | COMPLETED 0:0 | Converged at 11; RDS + 6 inspected PNGs; no independent object checks | [Record](data/integration/harmony/job-44437561/.INFO) |
| — / Oct 8 | Planned Harmony, maximum 50 | Intended 44416449 | **NOT SUBMITTED** | Skipped by user after 30-ceiling run converged | No run |

## QC sensitivity analyses

Both analyses reuse **44415373 labeled metrics, Scrublet calls and sample UMAPs**; neither creates a replacement clean RDS. CellBender 44436358 was explicitly excluded as the target.

| Job / date | Comparison | Execution / checks | Outputs |
|---|---|---|---|
| 44437522 / Oct 8 | Mitochondrial 15% cap, ribosomal +2 MAD, or both | COMPLETED 0:0; baseline reproduced; 92 sample/scenario rows checked; 18 PNGs inspected | [Record and tables](results/qc/sensitivity-job-44437522/), [log](logs/qc-sensitivity-44437522.log) |
| 44437686 / Oct 8 | Independent all-input thresholds; old vs new limits | COMPLETED 0:0; baseline/previous sensitivity reproduced; 92 rows checked; 12 PNGs inspected | [Record and tables](results/qc/sensitivity-job-44437686/), [log](logs/qc-independent-44437686.log) |

| Scenario | Analysis | Clean cells | Input retained |
|---|---|---:|---:|
| Historical baseline | 44415373 | 247,672 | 82.72% |
| Mitochondrial cap 15% only | 44437522 | 245,599 | 82.03% |
| Ribosomal +2 MAD only | 44437522 | 239,426 | 79.97% |
| Both, historical threshold reference | 44437522 | 237,410 | 79.30% |
| Independent all-input estimation, old limits | 44437686 | 247,282 | 82.59% |
| Independent all-input estimation, new combined limits | 44437686; adopted in 44438067 | **230,817** | **77.10%** |

All 23 samples remain. All-input estimation alone loses 413 and gains 23 cells; retained sets are **not necessarily nested**. Retention changes do not establish biological improvement.

## SoupX pilot and comparisons

Pilot input: REVO30-P4 raw RSEC archives, 13,422 original calls, 28,246 genes. Background: 0 < molecules <100, rank proxy >1,000, ≥2,000 non-called background barcodes, no called cells in range. Preliminary clustering: LogNormalize, 2,000 variable genes, PCs 1–20, k=20, Louvain 0.6, seed 1234. SoupX 1.6.2: automatic rho, subtraction, seeded integer rounding.

| Job / date | Purpose / parents | Execution / outcome | Availability / record |
|---|---|---|---|
| 44404957 / Oct 7 | First raw-input pilot | FAILED 1:0; gzip connection handling | Data folder later deleted; historical diagnostics retained |
| 44404963 / Oct 7 | Corrected pilot reader; raw input | COMPLETED 0:0; rho 0.031; 3.10% molecules removed; qualified marker review later passed | Data later deleted; [diagnostics](results/soupx/job-44404963/REVO30-P4/), [log](logs/soupx-44404963.log) |
| 44405863 / Oct 7 | Raw + SoupX 44404963 + CellBender 44400520 | COMPLETED 0:0; matched-cell comparison and marker tables checked; 3 figures inspected | [Review](results/soupx/job-44405863/REVO30-P4/) |
| 44408206 / Oct 7 | Raw + SoupX 44404963 + CellBender 44403564 | COMPLETED 0:0; updated comparator; counts/IDs/call flags checked; 3 figures inspected | [Review](results/soupx/job-44408206/REVO30-P4/) |
| 44409555 / Oct 7 | Initial 23-sample raw-input array | First four FAILED at plotting; remaining tasks cancelled | Eight data folders later deleted; historical diagnostics retained |

| REVO30-P4 method | Molecules retained | Removed | Positive-reference marker retention | Negative-reference marker removal |
|---|---:|---:|---:|---:|
| Raw | 84,811,621 | — | — | — |
| SoupX 44404963 / 44409807 | 82,183,607 | 3.09865% | 99.80% | 35.40% |
| CellBender 44400520 | 81,406,722 | 4.01466% | 99.77% | 40.01% |
| CellBender 44403564 | 81,056,997 | 4.42702% | 99.72% | 40.71% |

- CellBender 44400520 had convergence warnings; 44403564 (learning rate 5e-05, 150 epochs) passed reviewed convergence diagnostics. Neither was biologically accepted.
- Comparison uses original called cells only: 9,426 / 5,419 additional CellBender calls excluded respectively. In 44403564, original cell 55584626 is uncalled with zero corrected molecules, explicitly retained; jointly called sensitivity gives 4.42274% removal.
- Corrected-total correlations ≈0.9998 coexist with negative removed-total correlations (−0.5740 / −0.5644): corrections are not equivalent. More removal is not necessarily better.
- Marker references are post-hoc: 26 genes, four overlap rho fitting, six unresolved clusters; excluding fitting markers gives similar retention. SoupX half-maximum spans are not confidence intervals.
- All-sample 44409807 reproduces the pilot counts exactly for REVO30-P4. All 23 fits and representative marker figures were inspected; numerical status remains `validated_pending_qc` for biological acceptance.

## Historical raw-input runs

Legacy flat-path attribution uses logs/accounting/timestamps, not immutable job metadata. Unknown parents remain unknown; overwritten files cannot be reconstructed from a current path.

| Job / date | Stage / parent | Execution | Output / qualification | Log |
|---|---|---|---|---|
| 44176755 / Sep 29 | Normalization; exact parent not recorded | COMPLETED | [SCT retained](data/norm_feat/job-44176755/sct.rds); LogNormalize overwritten/unavailable | [Log](logs/20260929-0143-norm-feat-pipeline.log) |
| 44194877 / Sep 30 | QC; raw filtered MEX | COMPLETED | Labeled, DoubletFinder-clean and Scrublet-clean checkpoints; incomplete historical QC settings | [Log](logs/20260930-0006-qc-pipeline.log) |
| 44212906 / Sep 30 | LogNormalize; exact parent not recorded | COMPLETED | [Checkpoint](data/norm_feat/job-44212906/lognorm.rds) | [Log](logs/20260930-1325-norm-feat-pipeline.log) |
| 44276132 / Oct 1 | QC; 23 raw filtered MEX archives | COMPLETED 0:0 | 299,393 input → 231,993 QC-passing; DoubletFinder-clean 228,077; Scrublet-clean 226,738; cell sets checked; historical BD-rate clamp limitation | [Log](logs/20261001-172502-qc-pipeline-44276132.log) |
| 44280920 / Oct 1 | LogNormalize; 44276132 Scrublet-clean | COMPLETED 0:0 | **Wrong regression specification**: absent `percent_mito`; only nFeature_RNA regressed; preserved, not downstream input | [Log](logs/20261001-195329-norm-feat-pipeline-44280920.log) |
| 44281196 / Oct 1 | Corrected LogNormalize; 44276132 Scrublet-clean | COMPLETED 0:0 | Both percent.mt and nFeature_RNA regressed; historical downstream input | [Log](logs/20261001-203622-norm-feat-pipeline-44281196.log) |

### Historical integration

All used 23 samples and seed 1234. Historical Harmony maximum=10; scVI used original RNA counts, 3,000 selected genes, two layers, negative-binomial likelihood, automatic epoch limit (35 in 44281712). Legacy diagnostics sampled 20,000 cells, later runs 50,000; iLISI perplexity=30, Euclidean sample-label silhouettes.

| Job / date | Method | Parent / cells | PCA / integrated axes | Execution / availability | Log |
|---|---|---|---|---|---|
| 44184472 / Sep 29 | Legacy Harmony | Mutable flat path; 262,936 | 30 / 30 | CANCELLED; no uniquely attributable output | [Log](logs/20260929-1358-integration-pipeline.log) |
| 44190509 / Sep 29 | Legacy Harmony, LogNormalize + SCT | Mutable flat paths; 262,936 each | 30 / 30 | COMPLETED; shared `data/integration/lognorm.rds` and `sct.rds`; attribution qualified | [Log](logs/20260929-1717-integration-pipeline.log) |
| 44193502 / Sep 29 | Harmony | Exact producer unknown; 262,936 | 30 / 30 | COMPLETED | [Log](logs/20260929-213551-integration-harmony-44193502.log) |
| 44193503 / Sep 29 | scVI | Exact producer unknown; 262,936 | 30 / 30 | COMPLETED | [Log](logs/20260929-213551-integration-scvi-44193503.log) |
| 44215817 / Sep 30 | scVI | Exact producer unknown; 262,303 | 30 / 30 | COMPLETED | [Log](logs/20260930-154712-integration-scvi-44215817.log) |
| 44215822 / Sep 30 | Harmony | Exact producer unknown; 262,303 | 30 / 30 | COMPLETED | [Log](logs/20260930-154738-integration-harmony-44215822.log) |
| 44227281 / Sep 30 | scVI | Exact producer unknown; 262,303 | 30 / 20 | COMPLETED | [Log](logs/20260930-223724-integration-scvi-44227281.log) |
| 44281709 / Oct 1 | Harmony | 44281196; 226,738 | 20 / 20 | COMPLETED 0:0; diagnostics retained | [Log](logs/20261001-220536-integration-harmony-44281709.log) |
| 44281712 / Oct 1 | scVI | 44281196; 226,738 | 20 / 20 | COMPLETED 0:0; diagnostics retained | [Log](logs/20261001-220558-integration-scvi-44281712.log) |
| 44282203 / Oct 1–2 | CCA | 44281196; 226,738 | 20 / 20 | COMPLETED 0:0; 23h03m57s; diagnostics and UMAP fits logged; no independent full-object audit at completion | [Log](logs/20261001-234739-integration-cca-44282203.log) |
| 44445961 / Oct 9 | CCA UMAP exploration, 15-task array | 44282203: CCA checkpoint and saved default UMAP; 226,738 | — / 20 | COMPLETED 15/15, exit 0:0; 1h10m10s overall; all comparison PNGs present; no object persistence or plot inspection | [Record](results/integration/cca/job-44282203/.INFO), [results](results/integration/cca/job-44282203/) |

### Historical clustering

U=unintegrated PCA, H=Harmony, S=scVI. A–C: 262,303 cells; D/E: 226,738. All jobs below completed; 44310241's integration parent was not recorded in this history.

| Profile | Axes U/H/S | Neighbors | Resolutions |
|---|---|---:|---|
| A | 30/30/30 | 20 | 0.4, 0.6, 0.8, 1.0 |
| B | 30/30/30 | 60 | 0.1, 0.2, 0.3, 0.4 |
| C | 30/30/20 | 100 | 0.05, 0.1, 0.15, 0.2, 0.25, 0.3 |
| D | 20/20/20 | 100 | 0.05, 0.1, 0.15, 0.2, 0.25, 0.3 |
| E | 20/20/20 | 100 | 0.4, 0.6, 0.8, 1.0 |

A–E settings: Annoy Euclidean, 50 trees, SNN pruning 1/15; Louvain, 10 starts/iterations; seed 1234; 50,000 diagnostic cells. Cluster counts follow resolution order, not an optimal-resolution recommendation.

| Job | Method / profile | Parent / input reduction | Cluster counts | Record |
|---|---|---|---|---|
| 44221464 | U / A | 44215822: preserved PCA | 34/40/43/51 | [Log](logs/20260930-192711-clustering-unintegrated-44221464.log) |
| 44221465 | H / A | 44215822: Harmony | 33/37/44/52 | [Log](logs/20260930-192741-clustering-harmony-44221465.log) |
| 44221467 | S / A | 44215817: scVI | 38/49/55/60 | [Log](logs/20260930-192742-clustering-scvi-44221467.log) |
| 44222698 | U / B | 44215822: preserved PCA | 19/24/29/32 | [Log](logs/20260930-211635-clustering-unintegrated-44222698.log) |
| 44222697 | H / B | 44215822: Harmony | 19/24/28/30 | [Log](logs/20260930-211633-clustering-harmony-44222697.log) |
| 44222701 | S / B | 44215817: scVI | 26/33/34/37 | [Log](logs/20260930-211654-clustering-scvi-44222701.log) |
| 44228331 | U / C | 44215822: preserved PCA | 15/19/21/24/25/25 | [Log](logs/20260930-235042-clustering-unintegrated-44228331.log) |
| 44228335 | H / C | 44215822: Harmony | 15/18/20/22/22/24 | [Log](logs/20260930-235047-clustering-harmony-44228335.log) |
| 44228330 | S / C | 44227281: scVI | 17/23/27/30/31/31 | [Log](logs/20260930-235041-clustering-scvi-44228330.log) |
| 44282105 | U / D | 44281196: PCA directly | 11/14/16/17/18/19 | [Bundle](results/clustering/unintegrated/job-44282105/analysis.rds) |
| 44282104 | H / D | 44281709: Harmony | 12/15/16/16/17/20 | [Bundle](results/clustering/harmony/job-44282104/analysis.rds) |
| 44282108 | S / D | 44281712: scVI | 15/20/23/28/31/32 | [Bundle](results/clustering/scvi/job-44282108/analysis.rds) |
| 44282547 | U / E | 44281196: PCA directly | 21/24/27/30 | [Bundle](results/clustering/unintegrated/job-44282547/analysis.rds) |
| 44282546 | H / E | 44281709: Harmony | 20/24/28/30 | [Bundle](results/clustering/harmony/job-44282546/analysis.rds) |
| 44282549 | S / E | 44281712: scVI | 33/38/42/47 | [Bundle](results/clustering/scvi/job-44282549/analysis.rds) |
| 44310241 | CCA; 0.1/0.2/0.3/0.4/0.6/0.8/1.0 | Integration parent not recorded | 16/21/25/25/29/33/36 | [Bundle](results/clustering/cca/job-44310241/analysis.rds) |

D/E bundles passed provenance, cell-set, settings and diagnostic checks; most full RDS files were not reopened. scVI E was independently reopened by 44282993. Profile E **plots omit resolution 1.0 only** (47 scVI clusters exceed 45-color palette); objects/CSV tables retain all four resolutions.

### Marker discovery on CCA 44310241

Input: 226,738 cells, 20,916 genes, 23 normalized layers. RNA Wilcoxon/Presto, positive markers, min.pct=0.25, logfc.threshold=0.25, no cell cap; Bonferroni-adjusted p<0.05 for selection, ranked by decreasing log2FC.

| Job / date | Purpose / input | Execution / limit |
|---|---|---|
| 44394946 / Oct 6 | Preflight on 44310241 | COMPLETED 0:0; full checkpoint read, 50-gene all-cell API test only |
| 44395117 / Oct 6 | Software/figure validation on 44310241 | COMPLETED 0:0; normalization, joining, test statistics, ranking and figures checked on restricted features |
| 44395220 / Oct 6 | Seven-resolution production markers on 44310241 | CANCELLED by user; [partial output](results/find_markers/cca/job-44395220/), [log](logs/20261006-154439-find-markers-44395220.log) |
| 44395226 / Oct 6 | Plots from 44395220 + checkpoint 44310241 | CANCELLED by user |
| 44396189 / Oct 6 | Production plot check | CANCELLED by user |

No completed production marker analysis. Subsequent CSV-only simplification is **not covered by earlier validation**; new production/plot runs await explicit authorization.

## Supporting jobs: plotting, checks and failed attempts

Targets below are validation/comparison inputs, not new analysis checkpoints. Grouped IDs retain individual outcomes.

### Historical plotting and retired augmentation

| Job | Inputs / task | Recorded outcome | Record |
|---|---|---|---|
| 44215826; 44220807 | scVI 44215817 + Harmony 44215822 comparison | Both COMPLETED; second reused output directory | [First log](logs/20260930-163630-integration-plots-44215826.log), [second](logs/20260930-183438-integration-plots-44220807.log) |
| 44227466 | scVI 44227281 + Harmony 44215822 comparison | COMPLETED | [Log](logs/20260930-232055-integration-plots-44227466.log) |
| 44282037 | scVI 44281712 + Harmony 44281709 comparison | COMPLETED 0:0; seven figures | [Manifest/results](results/integration/comparison/scvi-job-44281712_harmony-job-44281709/) |
| 44310707 / 44310708 / 44310709 | Retired augmentation of scVI 44281712 / Harmony 44281709 / CCA 44282203 respectively | Last recorded RUNNING Oct 3 20:49; augmented copies later removed; terminal states unknown | `logs/20261003-204712-add-tsne-<ID>.log` |
| 44310723 | Comparison depending on all three augmentation jobs | Last recorded PENDING Oct 3; comparison directory later removed; terminal state unknown | [Log](logs/integration-plots-44310723.out) |
| 44221470 | Clustering A: U 44221464 / H 44221465 / S 44221467 | FAILED 1:0; >45-color guard; partial comparisons only | [Log](logs/20260930-195423-clustering-plots-44221470.log) |
| 44222706 | Clustering B: U 44222698 / H 44222697 / S 44222701 | COMPLETED | [Log](logs/20260930-215708-clustering-plots-44222706.log) |
| 44228340 | Clustering C: U 44228331 / H 44228335 / S 44228330 | COMPLETED | [Log](logs/20261001-010945-clustering-plots-44228340.log) |
| 44282442 | Clustering D: U 44282105 / H 44282104 / S 44282108 | COMPLETED 0:0; full comparisons | [Log](logs/20261002-004357-clustering-plots-44282442.log) |
| 44283792 | Clustering E: U 44282547 / H 44282546 / S 44282549 | COMPLETED 0:0; plots 0.4/0.6/0.8, tables all four resolutions | [Log](logs/20261002-065719-clustering-plots-44283792.log) |
| 44282993 | Full saved scVI E checkpoint 44282549 | COMPLETED 0:0; assignments, 20-D latent space and UMAP provenance matched | [Log](logs/scvi-checkpoint-verify-44282993.out) |

### SoupX and QC checks

| Job | Inputs / target | Recorded outcome | Record |
|---|---|---|---|
| 44405241 / 44405362 / 44405494 | Pilot review: SoupX 44404963 + CellBender 44400520 | Each FAILED 1:0: plot-order arguments / H5 latent-group reader / venv symlink respectively; replaced by 44405863 | `results/soupx/job-<ID>/REVO30-P4/` |
| 44408133 | Pilot review with CellBender 44403564 | FAILED 1:0; partial argument matching in dot-size settings; replaced by 44408206 | `results/soupx/job-44408133/REVO30-P4/` |
| 44409538; 44409543 | SoupX plotting preflight | First FAILED (Assay5 namespace); second COMPLETED 0:0, marker totals checked | `logs/soupx-review-<ID>.log` |
| 44409564 | Dependent check for initial array 44409555 | CANCELLED with array | Historical cancellation record |
| 44409746 | Final SoupX plotting-helper preflight | COMPLETED 0:0; restored nonempty Features label; figures inspected | [Log](logs/soupx-review-44409746.log) |
| 44409813 | All 23 saved count/marker outputs from 44409807 | COMPLETED 0:0; independent sparse-count and marker arithmetic checks passed | [Log](logs/soupx-verify-44409813.log) |
| 44410208 | Plot-only array on 44409807 | All 23 COMPLETED 0:0; widened plots, numerical summaries unchanged | `logs/soupx-plots-44410208_<task>.log` |
| 44411222 | SoupX-to-Seurat loader, all 23 samples | COMPLETED 0:0; exact counts/IDs/metadata, MT names, merge and raw-loader regression passed | [Log](logs/qc-soupx-io-test-44411222.log) |
| 44413862; 44413925 | Upper-ribosomal QC test | First FAILED on named-vector expectation; diagnostic confirmed matching values | Historical test record |
| 44413950; 44414023 | 4-MAD / final 3-MAD ribosomal tests on SoupX 44409807 | 3-MAD verification COMPLETED 0:0; 222 / 1,401 criterion flags respectively; not unique QC exclusions | [Final test log](logs/qc-ribo-mad3-test-44414023.log) |
| 44414192 | Two-sample QC pilot array on SoupX 44409807 | Two samples completed; parallel workflow withdrawn | Disposable outputs removed |
| 44414252; 44414254; 44414264 | Pilot collector / dependent check / read check | First two CANCELLED; third FAILED reading interrupted clean-RDS write | Disposable invalid output removed |
| 44414269 | Sequential QC on REVO26-C/REVO30-P4, SoupX 44409807 | COMPLETED 0:0; counts/flags/Scrublet/gene 9–10 boundary verified; clean 10,294 / 11,945 | Disposable test RDS removed |
| 44415118 | Saved QC 44415111 | FAILED 1:0; compared historical settings with mutable params; remaining validation incomplete, no demonstrated data defect | [Log](logs/qc-soupx-output-check-44415118.log) |
| 44415439 | Saved QC 44415373 | COMPLETED 0:0; exact input counts, saved settings, independent flags/retention/gene filtering passed | [Log](logs/qc-soupx-verify-44415439.log) |
| 44437549 | Sensitivity 44437522 tables/plots | COMPLETED 0:0; independent checks and UMAP legend refinement | [Log](logs/qc-sensitivity-plotcheck-44437549.log) |
| 44438061; 44438064 | Production QC preflight, REVO29-N4 | First FAILED on named-vector test; second COMPLETED 0:0; thresholds, equality, independence and sensitivity agreement passed | `logs/qc-preflight-<ID>.log` |
| 44438069 | Production QC 44438067 | FAILED: assumed historical Scrublet scores identical | [Log](logs/qc-combined-verify-44438069.log) |
| 44438245 | Score diagnostic: 44438067 vs historical scores | BD-rate table additions explain REVO30-P10 score change; **no doublet calls changed**; other 22 samples identical | [Log](logs/qc-score-check-44438245.log) |
| 44438254 | Production QC 44438067 | FAILED: assumed source/merged gene order identical; corrected check to align by IDs | [Log](logs/qc-combined-verify-44438254.log) |
| 44438256 | Production QC 44438067 | COMPLETED 0:0; saved parameters, thresholds, counts/IDs/metadata, flags and gene boundary passed; exact sensitivity agreement | [Log](logs/qc-combined-verify-44438256.log) |

### Integration checks

| Job | Inputs / target | Recorded outcome | Record |
|---|---|---|---|
| 44417083 | Plot-only CCA 44282203 | COMPLETED 0:0; six PNGs; inspected sample and faceted CCA plots; no new integration | [Log](logs/integration-plot-test-44417083.log) |
| 44417228 | CCA plotting rejection cases | COMPLETED 0:0; existing output directory and missing UMAP correctly rejected | [Log](logs/integration-plot-neg-44417228.log) |
| 44436241 | Saved Harmony 44435612 | COMPLETED 0:0; IDs/metadata/RNA/PCA, finite aligned embeddings and diagnostics checked; completed before cancellation request | [Log](logs/harmony-verify-44436241.log) |

## Reading and maintaining this index

- **Update a row, do not append a narrative.** Record job/date, stage, explicit parent/input role, status, one-line outcome and an existing evidence link. Keep lineage synchronized when a branch changes.
- Detailed parameters, seeds, timing, resources and checks belong in per-run `.INFO`, saved execution metadata and original logs. Current source settings are not historical evidence.
- Separate execution from output availability and computational checks from biological acceptance. Unknown parents/statuses stay unknown; file presence alone does not prove completion.
- No extra resource-consuming validation/check jobs for successful runs unless investigating an actual failure or explicitly requested. Coordinate jobs with their owning threads.
- Historical SVG copies were removed Oct 7; PNGs remain where outputs were retained. Old shared-path figures may mix runs. Links here are evidence locations, not a fresh audit of all artifacts.
- Clusters, mixing metrics, convergence and visually plausible UMAPs do not establish cell identities or biological preservation. Use [README](README.md) for current launch commands.
