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
| 44468374 / Oct 10 | Adaptive-QC validation: 5 REVO26 samples; soupx; harmony, dims 1–20, max_iter 10; per-sample median + 3 MAD capped (features 5000, mito 20%), ribo floor 2%, complexity >0.8 | COMPLETED, ~23 min | 65,965 input -> 52,186 QC-passing -> 51,026 clean (77.35%). Exclusion verified by independent reopen: no flagged or doublet cell in the final object; per-sample thresholds recorded; all 28 outputs present; QC panels with new captions inspected. Test artifacts retained | [INFO](../data/job-44468374/.INFO) |
| 44468410 / Oct 10 | Full production: 23 samples; soupx; harmony, dims 1–20, max_iter 10; adaptive QC as above | COMPLETED, ~1 h 48 min | 299,393 input -> 233,604 QC-passing (78.0%) -> 228,563 Scrublet-clean (76.34%); 20,984 genes. 3,000 consensus HVGs; 50 PCs (PC1–30 = 91.9%). Harmony converged after 8 iterations. Leiden 13/16/16/20 at 0.1–0.4. Independent reopen: exclusion verified (no flagged/doublet cells), per-sample thresholds recorded, all 46 outputs present, counts match .INFO; integration and cluster figures inspected. Final object: 8.0 GB clustered.rds | [INFO](../data/job-44468410/.INFO) |
| 44467784 / Oct 10 | Smoke test (fixed-threshold QC, superseded): 5 REVO26 samples; soupx; harmony, dims 1–20, max_iter 10 | COMPLETED, ~22 min | 65,965 input -> 51,723 QC-passing -> 50,675 Scrublet-clean (76.82%). Validation (separate jobs): final object reopened independently — IDs, per-sample layers, finite aligned reductions, all resolution columns, CSV agreement; old-vs-new QC flags identical on all input cells of REVO26-C and REVO26-P4; rhapsody and cellbender readers spot-checked; invalid-method and missing-input configs fail before run directories are created. Figures inspected. Run folders and logs deleted by request after validation | none (removed) |
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
