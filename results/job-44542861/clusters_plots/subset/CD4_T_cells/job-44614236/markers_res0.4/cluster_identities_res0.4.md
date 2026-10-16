# CD4_T_cells Subclustering Identity Report (`job-44614236`, Resolution `0.4`)

- **Parent Run / Subset**: `job-44542861` $\rightarrow$ `CD4_T_cells` (`job-44614236`; Parent Cluster 2)
- **Clustering Parameters**: `harmony_snn_res.0.4` (`25,888` cells across `6` clusters; `2,000` HVGs, `dims 1:20`, `k_param = 100`)
- **QC & Reassignment Summary**: **`0` doublet clusters**.
  - **Keep `19,105` true `CD4⁺` T cells (Clusters 1, 2, 5, 6)** (+ merge `1,669` helper `CD4⁺` T cells from `T_cells_unresolved` C4 $\rightarrow$ `20,774` total `CD4⁺` T cells).
  - **Reassign `6,783` co-clustered `CD8⁺` T cells (Cluster 3 = `3,748` Naive/CM `CD8A⁺ CD8B⁺`; Cluster 4 = `3,035` `GZMK⁺ CCL5⁺ CD8A⁺` EM) $\rightarrow$ `CD8_T_cells`**.
  - **Parameter Recommendation**: Re-subcluster purified `CD4_T_cells` (`20,774` cells) with **`k_param = 30–40`** (`dims = 1:20`) to split `FOXP3⁺` Tregs from `ICOS⁺ BCL6⁺` Tfh and `TOP2A⁺` cycling cells in C5.
- **Master Index & References**: See `../../NOTES_BARCODE_REASSIGNMENT_AND_CLEANUP.md` (Herrera-Uribe et al. 2021 `PMID: 34249103`; Gerner et al. 2022 `PMID: 34724393`; Hoog et al. 2022 `PMID: 35667468`).

---

## Summary Table (`harmony_snn_res.0.4`)

| Cluster | Cells (%) | Proposed Known Subset | QC / Doublet Flag | Key Named Top-25 Markers | Additional Canonical Markers (`markers_significant.csv`) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **1** | 9,970 (38.51%) | **Effector-Memory ($T_{EM}$) / Tissue-Resident Th1/Th17-like `CD4⁺` T Cells (incl. porcine `CD4⁺CD8α⁺` DP memory)** | Keep (Clean) | `PRDM1` (`46%`), `CCR5` (`35%`), `CXCR3` (`27%`), `CD40LG` (`47%`), `MAF` (`63%`), `BHLHE40` (`36%`), `AHR` (`31%`), `ITGB1` (`78%`), `S100A4`, `VIM`, `GADD45B`, `ANXA2` | `TNF` (#26, `26%`), `GPR183` (#28, `79%`), `IL2RB` (#32, `51%`), `ID2` (#33, `74%`), `IL18R1` (#34, `47%`), `TNFRSF1B` (#48, `58%`), `GATA3` (#49, `49%`), `CD4` (#69, `38%`); `CD8A` in `23%` (porcine DP memory) |
| **2** | 4,466 (17.25%) | **Naive / Central-Memory ($T_N / T_{CM}$) `CD4⁺` T Cells (`CCR7^{hi} SELL⁺ LEF1^{hi} TCF7^{hi}`)** | Keep (Clean) | `CCR7` (`85%`), `SELL` (`CD62L`, `49%`), `LEF1` (`81%`), `SATB1` (`90%`), `IL6ST` (`72%`), `IFNGR2` (`47%`), `IFIT1` (`42%`), `GCNT4`, `SLC1A4`, `ICAM2` | `TCF7` (#34, `89%`), `GBP5` (#43, `56%`), `STAT1` (#108, `88%`), `ISG15` (#137, `62%`), `KLF2` (#169, `64%`), `BACH2` (#227, `49%`); `CD8A⁻` (`3.3%`) |
| **3** | 3,748 (14.48%) | **Naive / Central-Memory `CD8A⁺ CD8B⁺` T Cells (`SELL⁺ CCR7^{hi} LEF1^{hi} TCF7^{hi} CRTAM^{hi}`)** | **Reassign (Not Doublet) $\rightarrow$ `CD8_T_cells`** | `CD8B` (`36%` vs `5%`), `GZMH` (`57%`), `GZMB` (`37%`), `PRF1` (`31%`), `CRTAM` (`78%`), `SELL` (`64%`), `CD27` (`34%`), `BCL2` (`91%`), `NR4A1` (`64%`), `JAML`, `GPR18`, `LIX1`, `CERS6` | `CD8A` (#28, `60%`), `LEF1` (#32, `84%`), `CCR7` (#45, `85%`), `KLF2` (#79, `74%`), `TCF7` (#355, `84%`); **`CD4⁻` (`1.1%`) & `CD40LG⁻` (`1.0%`)** |
| **4** | 3,035 (11.72%) | **`GZMK⁺ CCL5⁺ CD8A⁺ CD8B⁺` Effector-Memory `CD8⁺` T Cells** | **Reassign (Not Doublet) $\rightarrow$ `CD8_T_cells`** | `CCL5` (`46%`), `SLAMF7` (`32%`), `CCL4` (`28%`), `GZMK` (`29%`), `JAML` (`42%`), `PRF1` (`27%`), `CD8A` (`55%`), `NKG7` (`29%`), `CTSW` (`54%`), `CRTAM` (`62%`), `IKZF2` (`43%`), `RUNX3` (`43%`), `HLA-DRA`, `SLA-DRB1`, `SLA-DQB1`, `IL2RB` | `CD74` (#31, `63%`), `EGR1` (#41, `33%`), `CCR5` (#50, `30%`); **`CD4⁻` (`3.6%`) & `CD40LG^{lo}` (`9.0%`)** |
| **5** | 2,806 (10.84%) | **Activated / Cycling Co-Stimulatory & Treg/Tfh-Enriched `CD4⁺` T Cells (`ICOS⁺ TNFRSF18/GITR⁺ IRF4⁺ TOP2A⁺`)** | Keep (Clean) | `TIFA` (`54%`), `TOP2A` (`28%`), `IL6ST` (`73%`), `CD80` (`42%`), `TNFRSF18` (`GITR`, `31%`), `IRF4` (`26%`), `S1PR1` (`71%`), `DAPP1`, `LGMN`, `HERC6`, `ISG20`, `MAP3K5`, `TBC1D4` | `IL6R` (#28, `55%`), `SOCS3` (#29, `77%`), `TACC3` (#37, `38%`), `MAF` (#39, `60%`), `ICOS` (#49, `51%`), `IL18R1` (#50, `51%`), `CD69` (#63, `85%`), `BCL6` (#66, `34%`), `SH2D1A` (#67, `29%`), `CD28` (#87, `82%`), `CD40LG` (#90, `41%`) |
| **6** | 1,863 (7.20%) | **Quiescent / Resting `BACH2⁺ BCL11B⁺ RUNX1⁺` `CD4⁺` T Cells** | Keep (Clean) | `THSD7B`, `MAML2` (`58%`), `SMYD3`, `ATXN1`, `WWOX`, `PREX1`, `SLC9A9`, `PDE8A`, `LPP`, `CHST11` (`70%`), `BCAS3`, `AHCYL2`, `SIK3`, `BTBD9`, `XYLT1` (`71%`), `CDK6`, `CAMK1D`, `PRKCB` | `MSI2` (#26, `63%`), `SKAP1` (#28, `78%`), `BCL11B` (#29, `39%`), `KLF12` (#36, `53%`), `AFF3` (#37, `75%`), `RUNX1` (#43, `67%`), `BACH2` (#82, `44%`), `TOX` (#88, `62%`), `CD4` (#431, `31%`) |

---

## Concise Cluster Rationale

* **True `CD4⁺` T-Cell Subsets (`Clusters 1, 2, 5, 6` = `19,105` cells)**:
  * **C1 (`9,970`)**: Dominant tissue-resident/effector-memory helper T-cell pool (`CD40LG⁺ PRDM1⁺ CCR5⁺ CXCR3⁺ MAF⁺ BHLHE40⁺ AHR⁺ ITGB1^{hi} TNF⁺`). Includes physiological porcine `CD4⁺ CD8αα⁺` (`CD8B⁻`) antigen-experienced memory helper T cells (`23%` `CD8A⁺`; Herrera-Uribe et al. 2021 `PMID: 34249103`).
  * **C2 (`4,466`)**: Recirculating naive and central-memory `CD4⁺ CD8α⁻` T cells (`CCR7^{hi} SELL⁺ LEF1^{hi} TCF7^{hi} SATB1^{hi} KLF2⁺`).
  * **C5 (`2,806`)**: Acutely activated and proliferating (`TOP2A⁺ TACC3⁺`) co-stimulatory `CD4⁺` T cells combining Treg (`TNFRSF18/GITR⁺ CD80⁺ IRF4⁺`) and Tfh (`ICOS⁺ BCL6⁺ MAF⁺ SH2D1A⁺ CD40LG⁺`; Hoog et al. 2022 `PMID: 35667468`) programs.
  * **C6 (`1,863`)**: Quiescent `CD4⁺` T cells enriched for chromatin/quiescence regulators (`BACH2⁺ BCL11B⁺ RUNX1⁺ KLF12⁺`).
* **Co-Clustered `CD8⁺` T-Cell Subsets (`Clusters 3 & 4` = `6,783` cells $\rightarrow$ Reassign to `CD8_T_cells`)**:
  * **C3 (`3,748`) & C4 (`3,035`)**: Co-clustered in Parent Cluster 2 due to shared `CCR7⁺ SELL⁺ LEF1⁺ TCF7⁺` naive/central-memory markers, but lack `CD4` (`1.1–3.6%`) and `CD40LG` (`1.0–9.0%`) and express `CD8A`/`CD8B`, `RUNX3`, `CRTAM`, `PRF1`, and `GZMH`/`GZMB` (C3) or `GZMK`/`CCL5` (C4).
