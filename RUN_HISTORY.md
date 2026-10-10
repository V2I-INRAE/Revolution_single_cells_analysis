# Pipeline run history

Single-entry pipeline: `sbatch scripts/sbatch_pipeline.sh` runs `main.R`
(QC -> normalization -> integration -> clustering in one session), configured
entirely by `R/config.R`. Each run produces `data/job-<id>/.INFO` +
`clustered.rds` and `results/job-<id>/.INFO` + `qc_plots/`, `norm_plots/`,
`integration_plots/`, `clusters_plots/`.

History reset Oct 10, 2026 after the single-entry refactor. The four-stage
pipeline history (through Oct 10, 2026) is preserved in this file's git
history; legacy run folders and logs were removed by request.

## Runs

| Job / date | Config / route | Execution | Outcome and checks | Record |
|---|---|---|---|---|
| 44467784 / Oct 10 | Smoke test: 5 REVO26 samples (C, N4, P4, N10, P10); soupx corrected counts; harmony, dims 1–20, max_iter 10; fixed QC thresholds; 3,000 HVGs, 50 PCs; Leiden 0.1–0.4 | COMPLETED, ~22 min | 65,965 input -> 51,723 QC-passing -> 50,675 Scrublet-clean (76.82%); 19,684 genes. Harmony converged after 6 iterations. Leiden clusters 11/12/14/15 at 0.1/0.2/0.3/0.4. Validation (separate jobs): final object reopened independently — IDs, per-sample layers, 3,000x50,675 scale.data, finite aligned PCA/Harmony/UMAP reductions, all resolution columns, CSV agreement; old-vs-new QC flags identical on all input cells of REVO26-C and REVO26-P4; rhapsody and cellbender readers spot-checked; invalid-method and missing-input configs fail before run directories are created. Figures inspected (QC before/after, per-sample QC UMAPs, integration UMAPs, resolution UMAPs, clustree). Run folders and logs deleted by request after validation | none (removed) |
| 44468216 / Oct 10 | Full 23-sample production attempt (soupx, harmony) | CANCELLED by user ~2 min in, during QC start | No outputs produced; run folders removed | none |

## Reading and maintaining this index

- One row per run: job/date, config, execution, one-line outcome, evidence link.
- `.INFO` files in `data/job-<id>/` and `results/job-<id>/` hold the full
  resolved config, inputs, per-stage outcomes and status; link them as Record.
- Separate execution from checks and computational checks from biological
  acceptance; completed means the pipeline finished, not that clusters are
  validated cell identities.
- Keep pipeline evidence (logs, run folders) unless removal is requested;
  a row without a Record link should say where evidence went.
