# Master Subset Annotation Index, Doublet/Reassignment Ledger & Handoff Guide (`job-44542861`)

This document is the central index for the 12 subclustered coarse compartments of `job-44542861` (`228,627` total cells across 23 porcine Control and EVLP lung samples). It links to all 12 subset identity reports and literature reviews, lists every doublet cluster and lineage reassignment, evaluates `k_param` / dimensionality settings per subset, provides execution instructions for the next agent, and indexes all primary literature references used across the annotations.

---

## 1. Directory of Subset Reports & Literature Dossiers

### A. Subset Identity Reports (`12 / 12` Complete)

| # | Coarse Subset | Subcluster Job | Selected Res | Total Cells | Clean Cells Kept | Doublets Removed | Reassigned Out / In | Subset Identity Report Path |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | **`B_cells`** | `job-44614229` | `0.4` (9 cl.) | `26,790` | `26,790` | `0` | `0` / `0` | `B_cells/job-44614229/markers_res0.4/cluster_identities_res0.4.md` |
| 2 | **`CD4_T_cells`** | `job-44614236` | `0.4` (6 cl.) | `25,888` | `20,774` | `0` | `-6,783` (C3, C4) / `+1,669` | `CD4_T_cells/job-44614236/markers_res0.4/cluster_identities_res0.4.md` |
| 3 | **`CD8_T_cells`** | `job-44614237` | `0.4` (6 cl.) | `20,552` | `38,518`* | `0` | `0` / `+17,966`* | `CD8_T_cells/job-44614237/markers_res0.4/cluster_identities_res0.4.md` |
| 4 | **`DC`** | `job-44614230` | `0.3` (4 cl.) | `7,436` | `7,436` | `0` | `0` / `0` | `DC/job-44614230/markers_res0.3/cluster_identities_res0.3.md` |
| 5 | **`Endothelial`** | `job-44614231` | `0.3` (6 cl.) | `14,154` | `13,836` | `318` (C6) | `0` / `0` | `Endothelial/job-44614231/markers_res0.3/cluster_identities_res0.3.md` |
| 6 | **`Epithelial`** | `job-44614232` | `0.4` (11 cl.) | `10,029` | `9,426` | `603` (C7) | `0` / `0` | `Epithelial/job-44614232/markers_res0.4/cluster_identities_res0.4.md` |
| 7 | **`Granulocytes`** | `job-44614233` | `0.3` (6 cl.) | `23,037` | `23,037` | `0` | `0` (`1,002` Mast opt.) | `Granulocytes/job-44614233/markers_res0.3/cluster_identities_res0.3.md` |
| 8 | **`Macrophage_Monocytes`** | `job-44614234` | `0.4` (10 cl.) | `42,214` | `41,524` | `690` (C9, C10) | `0` / `0` | `Macrophage_Monocytes/job-44614234/markers_res0.4/cluster_identities_res0.4.md` |
| 9 | **`NK_cells`** | `job-44614238` | `0.4` (7 cl.) | `24,342` | `19,389` | `0` | `-4,953` (C3, C4) / `0` | `NK_cells/job-44614238/markers_res0.4/cluster_identities_res0.4.md` |
| 10 | **`Stromal`** | `job-44614235` | `0.4` (7 cl.) | `11,398` | `11,079` | `319` (C6) | `0` / `0` | `Stromal/job-44614235/markers_res0.4/cluster_identities_res0.4.md` |
| 11 | **`unconventional_T_cells`** | `job-44614239` | `0.4` (8 cl.) | `20,852` | `22,374`* | `0` | `0` (`399` ILC2 opt.) / `+1,522` | `unconventional_T_cells/job-44614239/markers_res0.4/cluster_identities_res0.4.md` |
| 12 | **`T_cells_unresolved`** | `job-44614240` | `0.4` (5 cl.) | `11,935` | **`0` (Dissolved)** | `2,514` (C2) | `-9,421` (C1, C3, C4, C5) | `T_cells_unresolved/job-44614240/markers_res0.4/cluster_identities_res0.4.md` |
| **Total** | **12 Subsets $\rightarrow$ 11 Clean** | — | — | **`228,627`** | **`224,183` (`98.06%`)** | **`4,444` (`1.94%`)** | **`21,157` reassigned** | — |

*\*Note: If `NK_cells` C3 & C4 (`4,953` `CD3⁺ CD8αα⁺` NK-like T / NKT cells) are assigned to `unconventional_T_cells` instead of `CD8_T_cells`, `CD8_T_cells` = `33,565` cells and `unconventional_T_cells` = `27,327` cells.*

### B. Associated Literature & Parent Annotation Dossiers
* **Porcine Lung Macrophage & Monocyte Literature Review**: `results/literature_review/porcine_lung_macrophage_monocyte_subsets/literature_report.md` (with `paper_summary_table.csv` and `search_log.md`)
* **Whole-Lung Parent (`res 0.7`) Literature Scratchpad**: `results/literature_review/SCRATCHPAD.md`
* **Whole-Lung Parent (`res 0.7`) Cluster Identities**: `results/job-44542861/clusters_plots/markers/0.7/cluster_identities_res0.7.md`

---

## 2. All Identified Doublet Clusters (`REMOVE (Doublet)` — `4,444` Cells Total)

*(Note: `B_cells`, `CD4_T_cells`, `CD8_T_cells`, `DC`, `Granulocytes`, `NK_cells`, and `unconventional_T_cells` had **0** doublet clusters.)*

| Subset | Subcluster Job | Resolution & Cluster ID | Cells (% of Subset) | Doublet Composition & Key Markers | Exact Extraction Rule (`subclustered.rds`) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **`Endothelial`** | `job-44614231` | `res 0.3` — **Cluster `6`** | `318` (`2.25%`) | **Endothelial + Leukocyte Doublets** (`PTPRC/CD45` `86.2%`, `LCP1`, `CD53`, `HCST`, `TYROBP`, `LAPTM5`, `CD3E`, `CCL5`, `S100A12` + pan-EC RNA) | `harmony_snn_res.0.3 == "6"` |
| **`Epithelial`** | `job-44614232` | `res 0.4` — **Cluster `7`** *(or `res 0.3` — **Cluster `6`**; identical barcodes)* | `603` (`6.01%`) | **Epithelial + T/NK Leukocyte Doublets** (`PTPRC/CD45` `75.1%`, `CD2` `42.6%`, `CD3D` `38.1%`, `CD3E` `35.0%`, `GZMA` `25.0%`, `CCL5` `44.3%` + diluted `EPCAM`/`SFTPA1`/`SCGB1A1`) | `harmony_snn_res.0.4 == "7"` *(or `res.0.3 == "6"`)* |
| **`Macrophage_Monocytes`** | `job-44614234` | `res 0.4` — **Cluster `10`** | `245` (`0.58%`) | **Multi-Lineage Doublet Soup** (`COL1A1/2/3⁺` Fibroblast + `SFTPC/SCGB1A1⁺` Epithelial + `JCHAIN⁺` Plasma + `GNLY/CD247⁺` T/NK) | `harmony_snn_res.0.4 == "10"` |
| **`Macrophage_Monocytes`** *(Borderline EC-adherent complex)* | `job-44614234` | `res 0.4` — **Cluster `9`** | `445` (`1.05%`) | **Patrolling Non-Classical Monocyte (`FCGR3A^{hi} CD14^{lo} CSF1R⁺ MAFB⁺ POU2F2⁺ KLF2⁺ HES1⁺ SIGLEC5⁺ ADGRG1^{hi}`) + Capillary EC (`PLVAP⁺ KIT⁺ PECAM1⁺ HPGD⁺ CDH5⁺`) Doublet / Complex** (`MARCO⁻ TIMD4⁻ VSIG4⁻ MRC1⁻`) | `harmony_snn_res.0.4 == "9"` |
| **`Stromal`** | `job-44614235` | `res 0.4` — **Cluster `6`** | `319` (`2.80%`) | **Stromal + Leukocyte (T/NK/Myeloid) Doublets** (`PTPRC/CD45` `55.2%`, `SRGN` `77.7%`, `LAPTM5` `59.9%`, `CXCR4` `54.9%`, `RGS1` `53.9%`, `FCGR3A` `38.2%`, `CCL5` `37.0%`, `CD2` `33.5%`, `CD3D` `27.6%` + diluted `COL1A1`/`COL3A1`/`DCN`) | `harmony_snn_res.0.4 == "6"` |
| **`T_cells_unresolved`** | `job-44614240` | `res 0.4` — **Cluster `2`** | `2,514` (`21.06%`) | **T-Cell + Alveolar Macrophage Doublets** (`APOE` `28.0%`, `FCN2` `25.8%`, `BTC` `28.1%`, `CD163` `35.0%`, `LGALS3` `32.9%`, `C3` `26.7%`, `CD9` `33.2%`, `CYBB` `27.0%`, `LGMN` `26.7%`, `ANXA1` `42.6%`, `NCF2` `25.0%` + diluted `CD3D`/`CD3E`/`CD3G`/`CD2`) | `harmony_snn_res.0.4 == "2"` |

---

## 3. All Confirmed Lineage Reassignments (`Reassign (Not Doublet)` — `21,157` Mandatory + `1,401` Optional)

| Source Subset | Subcluster Job | Resolution & Cluster ID | Cells (% of Subset) | True Biological Identity & Key Markers | Target Coarse Subset |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **`CD4_T_cells`** | `job-44614236` | `res 0.4` — **Cluster `3`** | `3,748` (`14.48%`) | **Naive / Central Memory `CD8⁺` T Cells** (`CD8B` `95.0%`, `CD8A` `71.0%`, `CD4` `0.8%`, `SELL` `69.6%`, `CCR7` `53.0%`, `LEF1`, `TCF7`) | **`CD8_T_cells`** |
| **`CD4_T_cells`** | `job-44614236` | `res 0.4` — **Cluster `4`** | `3,035` (`11.72%`) | **`GZMK⁺ CCL5⁺` Effector-Memory `CD8⁺` T Cells** (`CD8A` `74.4%`, `CD8B` `62.2%`, `CD4` `4.1%`, `CCL5` `90.9%`, `GZMK` `61.3%`, `KLRK1` `57.2%`, `CXCR3`) | **`CD8_T_cells`** |
| **`NK_cells`** | `job-44614238` | `res 0.4` — **Cluster `3`** | `2,892` (`11.88%`) | **Activated `CD3⁺ CD8αα⁺ FCGR3A(CD16)⁺ GNLY^{hi} CSF1⁺ MHC-II⁺` Cytotoxic NK-like T / NKT Cells** (`CD3D` `75.1%`, `CD3E` `70.1%`, `CD3G` `63.9%`, `CD247` `65.9%`, `CD8A` `63.3%`, `CD8B⁻`, `TYROBP⁻` `8.3%`) | **`CD8_T_cells`** *(or `unconventional_T_cells`)* |
| **`NK_cells`** | `job-44614238` | `res 0.4` — **Cluster `4`** | `2,061` (`8.47%`) | **Resident-Memory `CD3⁺ CD8αα⁺ TCF7⁺ CRTAM^{hi} ZNF683⁺ GZMM⁺` Cytotoxic T / NKT Cells** (`CD3D` `80.5%`, `CD3E` `75.3%`, `CD3G` `72.1%`, `CD247` `72.3%`, `CD8A` `60.8%`, `CD8B⁻`, `TYROBP⁻` `7.7%`) | **`CD8_T_cells`** *(or `unconventional_T_cells`)* |
| **`T_cells_unresolved`** | `job-44614240` | `res 0.4` — **Cluster `1`** | `4,017` (`33.66%`) | **Parenchymal Tissue-Resident / Effector-Memory `CD8⁺` T Cells ($T_{EM}/T_{RM}$)** (`GZMK` `31.9%`, `GZMA` `57.8%`, `CXCR6` `39.2%`, `CCR5` `46.6%`, `CRTAM` `45.1%`, `CD96` `32.7%`, `JAML` `26.1%`, `KLRK1` `50.9%`, `CD8A` `56.9%`, `CCL5` `86.9%`, `CD69` `62.8%`, `RUNX3` `41.0%`) | **`CD8_T_cells`** |
| **`T_cells_unresolved`** | `job-44614240` | `res 0.4` — **Cluster `3`** | `2,213` (`18.54%`) | **Intravascular Cytolytic Effector `CD8⁺` T Cells ($T_{EMRA}/T_{EFF}$)** (`ADGRG1` `70.2%`, `GZMB` `64.8%`, `GNLY` `61.3%`, `GZMH` `47.9%`, `PRF1` `50.0%`, `NKG7` `63.7%`, `FCGR3A` `87.0%`, `TBX21` `47.7%`, `KLF2` `78.7%`, `KLF3` `37.7%`, `S1PR1` `34.0%`, `IL10` `26.3%`, `CD8A` `54.0%`, `CD2` `82.0%`, `CD3D` `64.0%`) | **`CD8_T_cells`** |
| **`T_cells_unresolved`** | `job-44614240` | `res 0.4` — **Cluster `4`** | `1,669` (`13.98%`) | **Helper / Effector-Memory `CD4⁺` T Cells** (`CD40LG` `33.9%`, `CD28` `43.3%`, `MAF` `47.9%`, `IL7R` `50.4%`, `IL6R` `27.4%`, `IL18R1` `33.1%`, `TCF7` `33.0%`, `CD5` `65.8%`, `CD6` `38.0%`, `GPR183` `68.0%`, `SOCS3` `46.3%`, `CD8A⁻` `18%`, `GNLY⁻` `2%`, `GZMB⁻` `4%`) | **`CD4_T_cells`** |
| **`T_cells_unresolved`** | `job-44614240` | `res 0.4` — **Cluster `5`** | `1,522` (`12.75%`) | **`CD2⁻ WC1⁺` Porcine $\gamma\delta$ T Cells** (`SOX13` `44.5%`, `RHEX` `33.4%`, `BLK` `37.4%`, `GATA3` `56.5%`, `MAF` `46.4%`, `ZBTB16` `28.2%`, `ID2` `86.0%`, `IL2RA` `34.2%`, `IL7R` `56.7%`, `AQP3` `51.9%`, `CGA` `39.4%`, `AHR` `32.0%`, `CD2⁻` `13%`, `LCK⁻` `4%`, `CD6⁻` `6%`) | **`unconventional_T_cells`** |
| **`Granulocytes`** *(Optional)* | `job-44614233` | `res 0.3` — **Cluster `4`** | `1,002` (`4.35%`) | **Resident Pulmonary Mast Cells** (`KIT` `46.9%`, `MS4A2` `38.4%`, `FCER1A` `30.3%`, `GATA2`, `PTGS1`, `HPGD`, `RAB27B`, `MLPH`, `VWA5A`) | Keep in `Granulocytes` (as `Mast_cells`) or promote to coarse `Mast_cells` |
| **`unconventional_T_cells`** *(Optional)* | `job-44614239` | `res 0.4` — **Cluster `8`** | `399` (`1.91%`) | **Pulmonary Group 2 Innate Lymphoid Cells (`ILC2`)** (`CD3D/E/G⁻`, `ID2^{hi}` `90%`, `GATA3^{hi}` `89%`, `KIT^{hi}` `84.7%`, `MAF^{hi}` `83%`, `IL2RA^{hi}` `64.2%`, `HLF^{hi}` `60.4%`, `CD2⁺` `60.2%`, `RORA⁺` `60%`, `ZBTB16⁺` `57%`, `HOPX^{hi}` `54.1%`, `IL7R⁺` `52%`, `IL1RL1/ST2⁺` `42%`, `PTGDR2/CRTH2⁺` `26.1%`) | Keep in `unconventional_T_cells` (as `ILC2`) or promote to coarse `ILC2` |

---

## 4. Assessment of `k_param`, `dims`, and Re-Subclustering Needs Across All 12 Subsets

### Why `k_param = 100` Affected Small/Medium Subsets
In `R/config.R`, `config$clustering$k_param = 100` was tuned for the `228,627`-cell whole-lung object. When `R/clusters/subcluster.R` runs on an individual subset (`7,436` to `42,214` cells), inheriting `k_param = 100` means each cell's SNN neighborhood spans `0.7%–1.35%` of the entire subset. While `k_param = 100` worked well for large compartments (`B_cells`, `Granulocytes`, `Macrophage_Monocytes`), it merges closely related sub-lineages or rare (<250-cell) states in smaller subsets.

### Subset-by-Subset Parameter Recommendations

| Subset (`job`) | Current Params (`n_hvg`, `dims`, `k_param`, `res`) | Should You Re-Subcluster? | Recommended `k_param` & `dims` If Re-Subclustering (and Why) |
| :--- | :--- | :--- | :--- |
| **`B_cells`** (`job-44614229`) | `2000` HVGs, `1:20`, `k=100`, `res 0.4` | **No (Keep as-is)** | `k=100` already cleanly resolved all 5 B-cell and 4 ASC/Plasma states (including `255` `AICDA⁺` GC B cells) with a 100% stable clustree plateau (`0.4–0.5`). |
| **`Granulocytes`** (`job-44614233`) | `500` HVGs, `1:10`, `k=100`, `res 0.3` | **No (Keep as-is)** | `k=100` already cleanly separated all 3 Neutrophil states, Eosinophils (`3,814`), Mast Cells (`1,002`), and Basophils (`271`). |
| **`DC`** (`job-44614230`) | `1500` HVGs, `1:20`, `k=100`, `res 0.3` | **Optional (`k_param` change recommended if splitting `cDC2` vs `DC3/MoDC`)** | **`k_param = 25–30`, `dims = 1:15`**: In `7,436` cells, `k=100` forced canonical `FCER1A⁺ CD1.1⁺ CD207⁺` `cDC2` and inflammatory `MRC1⁺ CD163⁺ FCGR3A⁺` `DC3/MoDC` into a single cluster (C2, `1,460` cells). Lowering `k_param` to `25–30` will cleanly split `cDC2` from `DC3/MoDC`. |
| **`Endothelial`** (`job-44614231`) | `2000` HVGs, `1:15`, `k=100`, `res 0.3` | **Optional after removing C6 (`318` doublets)** | **`k_param = 30`, `dims = 1:15`**: `res 0.3` (`k=100`) already resolved LEC, Venous, `gCap`, `aCap` (Aerocytes), and Arterial ECs. However, if re-running after removing C6 doublets, lowering `k_param` to `30` will allow separation of **Resting Venous (`ACKR1⁺ HDAC9⁺ SELE⁻`)** vs. **EVLP-Activated Venous (`ACKR1⁺ SELE⁺ VCAM1⁺`)** and **Bronchial Systemic Venous (`COL15A1⁺ VWA1⁺ PLVAP⁺`)** ECs. |
| **`Stromal`** (`job-44614235`) | `2000` HVGs, `1:15`, `k=100`, `res 0.4` | **Optional after removing C6 (`319` doublets)** | **`k_param = 25–30`, `dims = 1:15`**: At `k=100`, C2 (`2,751` cells) co-clusters `WIF1⁺ GREM1⁺ PDGFRA⁺` Alveolar Ductal Myofibroblasts with `MYH11^{hi} ACTA2^{hi}` Airway Smooth Muscle (`ASM`), and C7 (`188` cells, smaller than `2 × k_param`) merges `ABCC9⁺ KCNJ8⁺` Pericytes with `ACTA2⁺ MYH11⁺` Vascular Smooth Muscle (`VSMC`). Lowering `k_param` to `25–30` will split Pericytes from VSMCs and Ductal Myofibroblasts from ASM. |
| **`Epithelial`** (`job-44614232`) | `2000` HVGs, `1:20`, `k=100`, `res 0.4` | **No (Just drop C7 `603` doublets)** | `res 0.4` (`k=100`) already resolved all 10 clean epithelial populations including rare Ionocytes (`292`), Basal (`258`), Tuft-1 (`320`), and Tuft-2 (`448`). Only lower `k_param` to `25–30` if re-subclustering after C7 removal to isolate ultra-rare `<100`-cell PNECs or Serous gland cells. |
| **`Macrophage_Monocytes`** (`job-44614234`) | `2000` HVGs, `1:20`, `k=100`, `res 0.4` | **No (Just drop C10 `245` & C9 `445` doublets)** | In `42,214` cells, `k=100` (`0.24%` of subset) is well-scaled and cleanly resolved 3 PAM states, Interstitial Macrophages, 3 Monocyte/MoMac states, and `MARCO^{hi} TIMD4⁺` Macrophages. If re-subclustering after dropping C9/C10, `k_param = 50` is optimal. |
| **`CD4_T_cells`** (`job-44614236`) | `2000` HVGs, `1:20`, `k=100`, `res 0.4` | **Yes (After removing C3+C4 `6,783` `CD8⁺` & adding `1,669` from `T_cells_unresolved` C4)** | **`k_param = 30–40`, `dims = 1:20`**: With `20,774` pure `CD4⁺` T cells (no longer distorted by `CD8⁺` cytotoxic PCs), lowering `k_param` to `30–40` will cleanly separate `FOXP3⁺ IL2RA^{hi} CTLA4^{hi}` Tregs from `ICOS^{hi} BCL6⁺ MAF^{hi}` Tfh and Cycling `TOP2A⁺` `CD4⁺` T cells (currently merged in C5). |
| **`CD8_T_cells`** (`job-44614237`) | `2000` HVGs, `1:20`, `k=100`, `res 0.4` | **Yes (After adding `+17,966` reassigned `CD8⁺` T cells)** | **`k_param = 40–50`, `dims = 1:20`**: Re-subclustering the unified `38,518`-cell `CD8⁺` T-cell compartment (combining Naive/CM, `GZMK⁺` $T_{EM}$, `ITGAE⁺` $T_{RM}$, `CX3CR1⁺ ADGRG1⁺` $T_{EMRA}$, `IL23R⁺ RORC⁺` Tc17/MAIT, and `CD8αα⁺` NKT) with `k_param = 40–50` will give a complete, unified `CD8⁺` T-cell trajectory. |
| **`NK_cells`** (`job-44614238`) | `2000` HVGs, `1:15`, `k=100`, `res 0.4` | **Optional (After removing C3+C4 `4,953` `CD3⁺` T/NKT cells)** | **`k_param = 30–40`, `dims = 1:15`**: Removing C3 & C4 leaves `19,389` pure `CD3E⁻` NK/ILC1 cells across 5 already-clean clusters (C1, C2, C5, C6, C7). Re-subclustering with `k_param = 30–40` is optional for a T-cell–free NK UMAP. |
| **`unconventional_T_cells`** (`job-44614239`) | `2000` HVGs, `1:20`, `k=100`, `res 0.4` | **Optional (If merging `+1,522` from `T_cells_unresolved` C5)** | **`k_param = 30–40`, `dims = 1:20`**: Existing `res 0.4` already cleanly resolves 6 `CD2⁻ WC1⁺` $\gamma\delta$ T states, `CD2⁺ EOMES^{hi}` $\gamma\delta$ T cells (C7, `455`), and Pulmonary `ILC2s` (C8, `399`). Re-subcluster with `k_param = 30–40` if integrating `T_cells_unresolved` C5 (`1,522` cells). |
| **`T_cells_unresolved`** (`job-44614240`) | `2000` HVGs, `1:20`, `k=100`, `res 0.4` | **Dissolve completely (`0` cells remain)** | Do not re-run: C2 (`2,514` doublets) is removed, and C1 (`4,017`), C3 (`2,213`), C4 (`1,669`), and C5 (`1,522`) are reassigned to `CD8_T_cells`, `CD4_T_cells`, and `unconventional_T_cells`. |

---

## 5. Handoff Guidance for the Next Agent: How to Proceed

Every subclustered Seurat object (`data/job-44542861/subsets/<Subset>/<job_id>/subclustered.rds`) preserves the exact same cell barcodes (`colnames(seurat_obj)`) as the parent whole-lung object (`data/job-44542861/annotated_clustered.rds`).

### Option A (Recommended): Parent Object Cleanup + Selective Re-Subclustering
1. **Extract Barcode Lists**:
   * Load the subclustered `.rds` objects listed in **Section 2 (Doublets)** and **Section 3 (Reassignments)** and extract the cell barcodes (`colnames(sub_obj)`) matching the specified resolution and cluster IDs.
2. **Update and Filter the Parent Seurat Object (`data/job-44542861/annotated_clustered.rds`)**:
   * Remove all `4,444` doublet barcodes from the parent object (`subset(parent, cells = setdiff(colnames(parent), doublet_barcodes))`).
   * Reassign `parent$subset` (and `Idents(parent)`) for the `21,157` T-cell barcodes according to **Section 3** (`T_cells_unresolved` is eliminated).
   * Save the updated parent object (e.g., `data/job-44542861/annotated_clustered_cleaned.rds`).
3. **Run Selective Re-Subclustering (`scripts/sbatch_subcluster.sh` + `scripts/sbatch_subset_markers.sh`)**:
   * Add a `--k-param` CLI override (or set `k_param = 35` in `R/config.R` for subclustering) so subclustering uses `k_param = 30–50` instead of `100`.
   * Re-run subclustering on the modified T/NK subsets (`CD4_T_cells`, `CD8_T_cells`, and optionally `NK_cells` and `unconventional_T_cells`), plus any non-T subset where finer `k_param = 25–30` splitting is desired (`DC`, `Stromal`, `Endothelial`).

### Option B (Lightweight / Fast): Direct Barcode-Level Metadata Lookup Table (No Re-Clustering)
If new UMAPs are not required and you want to proceed immediately to compositional (`run_group_distribution.R`) or differential expression analyses:
1. Loop through all 12 `subclustered.rds` objects and build a master lookup table keyed by `barcode` with columns: `orig_subset`, `subcluster_id`, `curated_coarse_subset`, `curated_fine_subset`, and `qc_status` (`"Keep"` vs. `"REMOVE_Doublet"`).
2. Assign each `(orig_subset, subcluster_id)` to its `Proposed Known Subset` from the 12 `cluster_identities_res<res>.md` reports and its `Target Coarse Subset` from **Section 3**.
3. Exclude `qc_status == "REMOVE_Doublet"` (`4,444` cells), yielding **`224,183` clean, annotated single cells** across 11 coarse compartments.

---

## 6. Master Literature Reference Index (By Lineage)

### A. Porcine & Mammalian Monocytes, Alveolar Macrophages (PAMs), IMs, PIMs & DCs
1. **Álvarez et al. (2023)** *Cells* (`PMID: 37626913`, `PMC10453229`, DOI: `10.3390/cells12162103`) — Comprehensive review of porcine monocyte/macrophage markers; establishes PAMs as `CD203a(ENPP1)^{hi} CLEC12B⁺ CD206(MRC1)^{hi} CD169(SIGLEC1)^{hi} CD163^{hi} CD172a(SIRPA)^{hi} FCGR1A^{hi} CD14^{lo} CD16^{lo}` (`Macrophage_Monocytes` C1–C3).
2. **Poderoso et al. (2020)** *Dev Comp Immunol* (`PMID: 32535044`, DOI: `10.1016/j.dci.2020.103767`) — Demonstrates `CLEC12B` is specifically expressed on porcine alveolar macrophages (`Macrophage_Monocytes` C1–C3) and absent on blood monocytes.
3. **Shin et al. (2025)** *Nat Commun* (`PMID: 39875369`, `PMC11775223`, DOI: `10.1038/s41467-024-54676-2`) — `306,052`-cell porcine BALF scRNA-seq atlas; defines steady-state `MRC1-NREP^{high}` (`PON1⁺ PPBP⁺ CD9⁺`) PAMs (`Macrophage_Monocytes` C1), lipid-associated PAMs (`C2`), and `PDCD1LG2^{high}` macrophages (`C7`).
4. **Chen et al. (2025)** *Front Immunol* (`PMID: 40242774`, `PMC12000058`, DOI: `10.3389/fimmu.2025.1574120`) — `114,881`-cell pig macrophage atlas across 49 organs; validates `PLSCR1` + `FCGR1A` as lung-specific PAM markers (`C1–C3`), `FABP5/LPL/GPNMB/SPP1` lipid-remodeling module (`C2`), and `C4BPA` in pig macrophages (`C7`).
5. **Zhang et al. (2026)** *Virulence* (`PMID: 41843736`, `PMC13014559`, DOI: `10.1080/21505594.2026.2646800`) — `59,062`-cell whole piglet lung scRNA-seq atlas; distinguishes `S100A8/9⁺ CD9⁺` AMs (`C1–C3`), `MARCO⁺ CD163⁺ C1QA/B/C⁺ HMOX1⁺` IMs (`C4`, `C7`), and `CD14⁺ FCN1⁺` monocytes (`C5`, `C6`).
6. **Ke et al. (2023)** *Microbiol Spectr* (`PMID: 37078873`, `PMC10269733`, DOI: `10.1128/spectrum.04753-22`) — Shows `MARCO` is constitutively expressed on primary porcine PAMs and restricts PRRSV (`Macrophage_Monocytes` C7).
7. **Zong et al. (2025)** *Virulence* (`PMID: 40944943`, `PMC12459362`, DOI: `10.1080/21505594.2025.2560894`) — Identifies `SLAMF7^{high}` porcine BAL macrophage subpopulation (`Macrophage_Monocytes` C7).
8. **Muir et al. (2024)** *PLoS Pathog* (`PMID: 39024231`, `PMC11257366`, DOI: `10.1371/journal.ppat.1011910`) — Spectral flow and scRNA-seq of Babraham pig BAL vs. PBMC myeloid and lymphoid cells (`Macrophage_Monocytes`, `CD4_T_cells`, `CD8_T_cells`).
9. **Maisonnasse et al. (2016a,b), Bordet et al. (2018), Crisci et al. (2020)** (`PMID: 26530136`, `PMID: 27992536` / `PMC5167224`, `PMID: 29977043` / `PMC6033929`, `PMID: 33365028` / `PMC7750501`) — Dissect the porcine lung DC/macrophage network (`cDC1`, `cDC2`, `moDC`, `moMacro`, `AM`, and `AM-like` `PIM`), showing a `10:1` AM:PIM ratio and near-identical bulk RNA-seq profiles between sorted AMs and PIMs (`DC`, `Macrophage_Monocytes`).
10. **Auray et al. (2016) & Baillou et al. (2025)** (`PMID: 27798170`; *Front Immunol* `PMID: 41235223` / `PMC12606228`, DOI: `10.3389/fimmu.2025.1639553`) — Define porcine `XCR1⁺ CADM1⁺` `cDC1`, `FCER1A⁺ CD1.1⁺` `cDC2`, `CD163⁺` `DC3`, and `CD4^{hi} TCF4⁺ IRF8⁺` `pDC` (`DC` C1–C4).
11. **Fairbairn et al. (2013)** *J Leukoc Biol* (`PMID: 23667115`, DOI: `10.1189/jlb.0113030`) — Establishes the `CSF1`-driven maturation series from `CD14^{hi} CD163^{lo}` (`Macrophage_Monocytes` C6) to `CD14⁺ CD163^{hi}` (`C5`, `C8`) porcine monocytes.
12. **Gouin et al. (2023)** *Front Immunol* (`PMID: 37465668`, `PMC10351384`, DOI: `10.3389/fimmu.2023.1142228`) — Human `0h, 4h, 10h` EVLP lung scRNA-seq; shows `CD14⁺ FCN1⁺ VCAN⁺` `cMos` and `MoMacs` drive `IL1B/IL1R1/IL1R2/AREG/EREG` expression during EVLP (`Macrophage_Monocytes` C5, C6).
13. **Mulder et al. (2021)** *Immunity* (`PMID: 34297922`, DOI: `10.1016/j.immuni.2021.07.007`) — MoMac-VERSE cross-tissue atlas defining conserved `TIMD4⁺ MARCO⁺ VSIG4⁺ C1Q⁺` resident macrophages, `C1Q^{hi} MAF⁺` interstitial macrophages, and `IL1B⁺ AREG⁺ EREG⁺` `MoMacs` (`Macrophage_Monocytes` C4, C5, C7).

### B. Porcine T Cells ($\alpha\beta$, $\gamma\delta$), NK Cells & Innate Lymphoid Cells (`ILCs`)
14. **Herrera-Uribe et al. (2021)** *Front Genet* (`PMID: 34249103`, `PMC8261551`, DOI: `10.3389/fgene.2021.689406`) — Bulk and single-cell reference atlas of porcine PBMCs; defines `CD4⁺ CD8α⁻` naive vs. `CD4⁺ CD8α⁺` activated/memory porcine `CD4⁺` T cells (`CD4_T_cells` C1 vs. C2), `CD2⁻ WC1⁺` vs. `CD2⁺` $\gamma\delta$ T cells (`unconventional_T_cells` C1–C6 vs. C7), and `CD21⁺`/`CD21⁻` B cells and ASCs (`B_cells`).
15. **Wiarda et al. (2022) & Wiarda et al. (2026)** *Life Sci Alliance* (`PMID: 35995567`, `PMC9396248`, DOI: `10.26508/lsa.202201442`) & *Front Immunol* (`PMID: 41859083`, `PMC12995613`, DOI: `10.3389/fimmu.2026.1704257`) — Single-cell atlases of porcine intestinal, primary, and secondary lymphoid organs; define `SELL^{hi} ID3⁺ ITGB7⁺` $\gamma\delta$ T cells (`unconventional_T_cells` C1), `CD2⁺ CD3E⁻ ID2^{hi} KIT^{hi} IL7R⁺` porcine tissue ILCs (`unconventional_T_cells` C8 `ILC2`), and `AICDA⁺ RGS13⁺` GC B cells (`B_cells` C9).
16. **Rodríguez-Gómez et al. (2016, 2019)** *Dev Comp Immunol* (`PMID: 26920461`) & *Front Immunol* (`PMID: 30915070`, `PMC6421308`, DOI: `10.3389/fimmu.2019.00396`) — Establish that `~80–90%` of porcine `CD2⁻` $\gamma\delta$ T cells in blood and lung are `GATA3^{hi} CD27⁺ Perforin⁻` (`unconventional_T_cells` C1–C6), whereas `EOMES⁺` $\gamma\delta$ T cells (`~14%` in lung) are `CD2⁺ CD27⁺` (`unconventional_T_cells` C7).
17. **Gu et al. (2022) & Zhang et al. (2024)** *Cell Rep* (`PMID: 35793622`, `PMC9704770`, DOI: `10.1016/j.celrep.2022.111050`) & *Front Immunol* (`PMID: 38384475`, `PMC10879363`, DOI: `10.3389/fimmu.2024.1339787`) — Single-cell atlases of porcine thymopoiesis and peripheral T cells; confirm `SOX13, BLK, GATA3, ETV5, CD163L1(WC1)` as hub regulators of `CD2⁻` $\gamma\delta$ T cells (`unconventional_T_cells` C1–C6) and `IKZF2` in `CD2⁺` $\gamma\delta$ T cells (`C7`).
18. **Takamatsu et al. (2006) & Le Page et al. (2022)** *Anim Health Res Rev* (`PMID: 17389056`) & *Dev Comp Immunol* (`PMID: 34863955`, DOI: `10.1016/j.dci.2021.104330`) — Reviews of porcine unconventional T cells, documenting `CD80/CD86` and MHC-II upregulation on activated APC-like porcine $\gamma\delta$ T cells (`unconventional_T_cells` C3, C6).
19. **Mair et al. (2014, 2016) & Gerner et al. (2009, 2022)** *Front Immunol* (`PMID: 27471504`, `PMC4943943`), *Dev Comp Immunol* (`PMID: 18601948`), *Annu Rev Anim Biosci* (`PMID: 34724393`) — Characterize porcine `CD3⁻ NKp46(NCR1)^{high}` tissue NK cells (`NK_cells` C1), `CD3⁻ NKp46^{low} CD8αα⁺ CD16^{high}` cytolytic NK cells (`NK_cells` C2), and `CD3⁺ NKp46⁺ CD8⁺` NK-like T cells enriched in pig lungs (`CD8_T_cells` C3, `NK_cells` C3–C4).

### C. Pulmonary Endothelial, Epithelial, Stromal & Granulocyte Lineages
20. **Gillich et al. (2020) & Schupp et al. (2021)** *Nature* (`PMID: 33057196`) & *Circulation* (`PMID: 34097468`) — Define the mammalian alveolar capillary dichotomy between `GPIHBP1⁺ APLNR⁺ HPGD⁺ KIT⁺` General Capillaries (`gCap` / CAP1; `Endothelial` C3) and `APLN⁺ EDNRB⁺ PRX⁺ CA4⁺ KITLG⁺ VWF⁻ TEK⁻` Alveolar Aerocytes (`aCap` / CAP2; `Endothelial` C4), plus `ACKR1⁺` Venous (`C2`), `GJA5⁺ BMX⁺ HEY2⁺` Arterial (`C5`), and `PROX1⁺ CCL21⁺` Lymphatic (`C1`) ECs.
21. **Travaglini et al. (2020) & Sikkema et al. (2023)** *Nature* (`PMID: 33208946`) & *Nat Med* (`PMID: 37291214`, Human Lung Cell Atlas) — Reference atlases for pulmonary epithelial (`AT1`, `AT2`, `Club`, `Multiciliated`, `Basal`, `FOXI1⁺ CFTR⁺` `Ionocytes`, `POU2F3⁺` `Tuft`), stromal (`COL13A1⁺ ITGA8⁺ NPNT⁺` Alveolar Fibroblasts, `PI16⁺ MFAP5⁺` Adventitial Fibroblasts, `WIF1⁺ GREM1⁺ MYH11⁺` Ductal Myofibroblasts/ASM, `NOTCH3⁺ PDGFRB⁺` Pericytes), and immune compartments.
22. **Montoro et al. (2018) & Perniss et al. (2023)** *Nature* (`PMID: 30069044`) & *Eur Respir J* (`PMID: 37531421`) — Define airway `POU2F3⁺` `Tuft-1` (`CHAT^{hi} DCLK2⁺` cholinergic/neuronal; `Epithelial` C9) vs. `Tuft-2` (`SUCNR1⁺ ALOX5AP^{hi} IL17RB^{hi} GFI1B⁺ PTPRC⁺` succinate-sensing immunological; `Epithelial` C8) brush cells and `FOXI1⁺` pulmonary ionocytes (`Epithelial` C10).
23. **Lee et al. (2017) & Narvaez Del Pilar et al. (2022)** *Cell* (`PMID: 28886383`) & *Development* (`PMID: 35302583`) — Define the three-axis lung mesenchymal organization and `LGR5⁺ FGF10^{hi} HGF^{hi} SFRP1^{hi}` epithelial-niche mesenchymal cells (`Stromal` C5).
