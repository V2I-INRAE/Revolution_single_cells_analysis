# Subset cell filtering report

| Subset | Before | T flags | B flags | Plasma flags | Removed once | Removed percentage | Retained |
|---|---:|---:|---:|---:|---:|---:|---:|
| T_cell_nk | 109494 | 0 | 78 | 12 | 90 | 0.082% | 109404 |
| B_cells | 29165 | 226 | 0 | 0 | 226 | 0.775% | 28939 |
| Myeloid | 54105 | 200 | 217 | 64 | 465 | 0.859% | 53640 |
| Granulocytes | 23532 | 54 | 23 | 7 | 83 | 0.353% | 23449 |
| Endothelial | 15389 | 124 | 73 | 25 | 219 | 1.423% | 15170 |
| Epithelial | 11315 | 212 | 90 | 28 | 327 | 2.890% | 10988 |
| Stromal | 12401 | 82 | 51 | 23 | 151 | 1.218% | 12250 |

The run removed 1561 unique cells from 255401 input cells and retained 253840 cells.

## Approved selection rules

Detection uses existing SoupX-corrected RNA counts >0; unavailable features are not biological zeros.
T pattern: CD3E and CD3D, together with CD5 or CD6.
B pattern: at least 3 of CD79A, CD79B, MS4A1, CD19, PAX5 and BANK1.
Plasma pattern: JCHAIN, MZB1 and TNFRSF17 together.
T patterns are allowed in T_cell_nk. B and plasma patterns are allowed in B_cells.
Remove the union of unexpected patterns. A cell matching several patterns is removed once.
The owner approved marker-based exclusions; flags are not confirmed cell identities or measured contamination rates.
No SoupX, normalization, feature filtering, HVG selection, PCA, integration or clustering is rerun.
Existing HVGs, scale layers, reductions, graphs and cluster labels are inherited, not refitted after removal.
Mixed has no saved subset RDS and is not filtered. Granulocytes_warning_check duplicates Granulocytes and is not processed.

## Inputs and outputs

The filtering producer is `job-44577467`. Original input RDS files were not overwritten.
- `data/job-44493014/subsets/T_cell_nk/subclustered.rds` becomes `data/job-44493014/subsets/T_cell_nk/job-44577467/subclustered.rds`.
- `data/job-44493014/subsets/B_cells/subclustered.rds` becomes `data/job-44493014/subsets/B_cells/job-44577467/subclustered.rds`.
- `data/job-44493014/subsets/Myeloid/subclustered.rds` becomes `data/job-44493014/subsets/Myeloid/job-44577467/subclustered.rds`.
- `data/job-44493014/subsets/Granulocytes/subclustered.rds` becomes `data/job-44493014/subsets/Granulocytes/job-44577467/subclustered.rds`.
- `data/job-44493014/subsets/Endothelial/subclustered.rds` becomes `data/job-44493014/subsets/Endothelial/job-44577467/subclustered.rds`.
- `data/job-44493014/subsets/Epithelial/subclustered.rds` becomes `data/job-44493014/subsets/Epithelial/job-44577467/subclustered.rds`.
- `data/job-44493014/subsets/Stromal/subclustered.rds` becomes `data/job-44493014/subsets/Stromal/job-44577467/subclustered.rds`.

Each subset has a cell_flags.csv file covering every input cell, with marker counts and exclusion reasons.
Each subset also has a removed_cells.csv file containing only excluded cell IDs and their evidence.
The by_sample.csv and by_parent_cluster.csv files report exclusion counts and retained counts.
The marker_features.csv file records the gene names and actual features used, including pig Ensembl ID resolutions.

## Validation and limitations

Every saved RDS was reopened. Retained cell IDs, metadata, identities, assays/count layers, reductions, graphs and neighbors matched the saved in-memory object.
All original input file sizes and modification times were unchanged. Required markers were available in every sample layer.
Flags are marker-based exclusion criteria, not experimentally confirmed identities or a measured contamination rate.
The inherited analysis and existing PCA figures still describe fits made before removal. Recompute HVGs, scaling and PCA from the new RDS files before choosing downstream PCs.
The original execution log is `logs/subset-filtering-44577467.log`.

