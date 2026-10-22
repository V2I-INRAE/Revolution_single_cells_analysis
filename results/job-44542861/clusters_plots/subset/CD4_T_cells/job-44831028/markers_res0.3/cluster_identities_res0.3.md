# CD4_T_cells Subclustering Identity Report (`job-44831028`, Resolution `0.3`)

- **Parent Run / Subset**: `job-44542861` $\rightarrow$ `CD4_T_cells` (`job-44831028`; cleaned baseline with `1,669` helper `CD4⁺` cells merged from `T_cells_unresolved` C4)
- **Clustering Parameters**: `harmony_snn_res.0.3` (`20,774` cells across `5` clusters; `2,000` HVGs, `dims 1:15`, `k_param = 100`)
- **QC & Parameter Assessment**:
  - **Keep `20,774` cells (Clusters 1–5). No clear mixed-lineage cluster identified from these marker exports.** `CD8A`/`CD8B` and B (`CD19`, `MS4A1`), plasma (`SDC1`, `JCHAIN`), monocyte (`S100A9`, `FCN1`, `SIRPA`, `CD14`, `CD163`), endothelial (`PECAM1`, `CLDN5`) and epithelial (`EPCAM`, `KRT19`) markers have no rows in the significant-marker export. C4 retains a partial CD8-associated contextual signal (`CCL5` #74, `CRTAM` #196, `RUNX3` #209) without `CD8A/B` or a coherent cytolytic program, so no residual CD8 cluster is established. Filtered exports cannot prove every cell is a singlet.
  - **Reassignment check**: prior run `job-44614236` (res `0.4`) resolved `25,888` cells in `6` clusters, of which C3 (`3,748`) and C4 (`3,035`) were co-clustered `CD8⁺` T cells reassigned to `CD8_T_cells`, plus a `1,669`-cell helper-`CD4⁺` merge from `T_cells_unresolved`. This cleaned run totals `19,105 + 1,669 = 20,774` cells in `5` clusters with no clear residual CD8-defined cluster detected; the reassignment worked operationally. Absence from the export does not prove no CD8-expressing cells remain.
  - **C5 (`1,167` cells) is unresolved.** Its top markers (`TXNRD1`, `ANXA1`, `CD9`, `SLC6A6`) carry small fold changes and low detection contrasts, with no convincing T-cell program; it is kept but not assigned an identity in this report.
  - **`k_param` Note**: `k_param = 100` at `res 0.3` resolves effector-memory, naive/central-memory, activated/Treg-Tfh-enriched, and quiescent `CD4⁺` programs. If resolving Treg/Tfh/cycling programs in C3 is desired, the superseding ledger recommends a `k=50` sensitivity comparison first, with targeted lower k only for cell-level-supported questions (untested here).
- **Master Index & References**: See `../../../NOTES_BARCODE_REASSIGNMENT_AND_CLEANUP.md` (Herrera-Uribe et al. 2021 `PMID: 34249103`; Gerner et al. 2022 `PMID: 34724393`; Hoog et al. 2022 `PMID: 35667468`).

---

## Summary Table (`harmony_snn_res.0.3`)

| Cluster | Cells (%) | Proposed Known Subset | QC / Doublet Flag | Key Named Top-25 Markers | Additional Canonical Markers (`markers_significant.csv`) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **1** | 10,376 (49.95%) | **Effector-Memory ($T_{EM}$) Helper `CD4⁺` T Cells (Th1-like, tissue residency provisional)** | Keep (Clean) | `FCGR3A`, `CCR5` (`34%`), `S100A4`, `CXCR3` (#4, `28%`), `AQP3`, `VIM`, `ITGB1` (#9, `77%`), `BHLHE40`, `PRDM1` (#12, `45%`), `GADD45B`, `AHR` (#25, `30%`), `CD40LG` (`47%`) | `TNF` (#27, `26%`), `GPR183` (#26, `78%`), `ID2` (#32, `75%`), `MAF` (#59, `61%`), `GATA3` (#121, `48%`) |
| **2** | 4,321 (20.80%) | **Naive / Central-Memory ($T_N / T_{CM}$) `CD4⁺` T Cells (`CCR7^{hi} SELL⁺ LEF1^{hi} TCF7^{hi}`)** | Keep (Clean) | `SELL` (`CD62L`, #1, `50%`), `CERS6`, `LEF1` (#10, `82%`), `GCNT4`, `CCR7` (#17, `85%`), `SATB1` (#23, `90%`), `IFIT1` (#25, `42%`), `IFNGR2` (#21, `47%`) | `TCF7` (#40, `89%`), `IL6ST` (#42, `72%`), `GBP5` (#70, `56%`), `KLF2` (#84, `65%`), `STAT1` (#120, `87%`), `BACH2` (#155, `50%`) |
| **3** | 3,273 (15.76%) | **Activated / Cycling Co-Stimulatory & Treg/Tfh-Enriched `CD4⁺` T Cells (`IL2RA⁺ TOP2A⁺ TNFRSF18/GITR⁺ IRF4⁺`)** | Keep (Clean) | `TIFA` (#1, `51%`), `IL2RA` (`CD25`, #2, `25%`), `TOP2A` (#3, `28%`), `LGMN`, `DAPP1`, `CD80` (#7, `39%`), `IRF4` (#10, `25%`), `TNFRSF18` (`GITR`, #16, `31%`), `S1PR1`, `IL6ST` (#13, `69%`) | `SOCS3` (#59, `75%`), `IL6R` (#67, `53%`), `ICOS` (#64, `50%`), `BCL6` (#104, `32%`), `MAF` (#101, `60%`), `CD69` (#63, `84%`), `CD28` (#89, `80%`); `FOXP3` absent from significant |
| **4** | 1,637 (7.88%) | **Distinct State, Quiescent / Resting Provisional (`BACH2⁺ BCL11B⁺ RUNX1⁺`)** | Keep (Clean) | `THSD7B` (#1), `SMYD3`, `MAML2` (#5, `61%`), `ATXN1`, `WWOX`, `PREX1`, `CHST11` (#18, `74%`), `MSI2` (#24, `66%`), `SKAP1` (#25, `83%`), `AFF3` (#22, `80%`) | `BCL11B` (#31, `42%`), `KLF12` (#33, `56%`), `RUNX1` (#58, `70%`), `BACH2` (#41, `45%`), `TOX` (#77, `65%`) |
| **5** | 1,167 (5.62%) | **Unresolved (kept, no identity assigned)** | Keep (Unresolved) | `TXNRD1`, `ANXA1`, `CD9`, `SLC6A6`, `GPCPD1` (all log2FC ≤ `2.12`, i.e. ~2.9–4.3-fold; detection contrasts ≤ `21.2` points; limited subtype specificity) | No coherent positively enriched T-cell subtype program (`CD3E`/`CD4` detected without positive enrichment) |

---

## Concise Cluster Rationale

* **C1 (`10,376`)**: Dominant effector-memory helper pool (`CD40LG⁺ PRDM1⁺ CCR5⁺ CXCR3⁺ MAF⁺ BHLHE40⁺ AHR⁺ ITGB1^{hi} TNF⁺ GATA3⁺`; tissue residency and finer Th subsetting provisional). `CD8A`/`CD8B` have no significant rows, so porcine `CD4⁺CD8α⁺` double-positive memory cannot be assessed from these CSVs.
* **C2 (`4,321`)**: Recirculating naive/central-memory `CD4⁺` T cells (`CCR7^{hi} SELL⁺ LEF1^{hi} TCF7^{hi} SATB1^{hi} KLF2⁺ BACH2⁺`), with interferon-response genes (`IFIT1`, `IFNGR2`, `GBP5`, `STAT1`).
* **C3 (`3,273`)**: Activated and proliferating (`IL2RA/CD25⁺ TOP2A⁺ IRF4⁺`) co-stimulatory `CD4⁺` T cells with Treg-associated (`TNFRSF18/GITR⁺ CD80⁺`) and Tfh-associated (`ICOS⁺ BCL6⁺ MAF⁺ CD28⁺ CD69⁺`) markers. `FOXP3` is absent from significant markers, so Treg identity is program-level, not master-regulator confirmed, and separate Treg/Tfh populations are not demonstrated.
* **C4 (`1,637`)**: Distinct transcriptional state, quiescent/resting provisional (`BACH2⁺ BCL11B⁺ RUNX1⁺ KLF12⁺ TOX⁺ SKAP1⁺`; `BACH2` detection is higher in C2, so BACH2 alone does not define this state).
* **C5 (`1,167`) — Unresolved**: Top markers show ~2.9–4.3-fold enrichment but limited detection coverage and weak subtype specificity, with no coherent positively enriched T-cell subtype program; kept in the object, excluded from identity claims pending further analysis.

> Note: pre-red-palette outputs (viridis dotplot + prior CSVs/heatmaps) are archived in `old/` in this directory; current plots use the red scale.
