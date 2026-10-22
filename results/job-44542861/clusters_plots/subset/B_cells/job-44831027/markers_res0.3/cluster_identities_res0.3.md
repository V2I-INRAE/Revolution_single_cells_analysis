# B_cells Subclustering Identity Report (`job-44831027`, Resolution `0.3`)

- **Parent Run / Subset**: `job-44542861` $\rightarrow$ `B_cells` (`job-44831027`; cleaned baseline; Parent Clusters 7 [`16,010` B cells] + 9 [`10,780` Plasma cells])
- **Clustering Parameters**: `harmony_snn_res.0.3` (`26,790` cells across `8` clusters; `1,500` HVGs, `dims 1:15`, `k_param = 100`)
- **QC & Parameter Assessment**:
  - **Keep all `26,790` cells (Clusters 1–8). No doublet cluster.** T-cell (`CD3D/E/G`), NK (`NCR1`, `GNLY`, `GZMB`) and monocyte (`S100A9`, `FCN1`, `SIRPA`, `CD14`, `CD163`) markers have no significant rows in any cluster. `PTPRC/CD45`, `LCP1`, `CD53`, `CXCR4` appear in significant rows with mixed directions (e.g. `PTPRC` C2 `+1.25`, C8 `+0.51`; C1 `−2.20`) on a broad B-lineage background, and `RGS1` enrichment marks the plasma branch (C1 `97%`); these lack foreign-lineage specificity and no cluster shows a coherent foreign-lineage program, so no cluster-level marker evidence currently warrants foreign-lineage reassignment or doublet-cluster exclusion. Filtered marker summaries cannot establish that every cell is a singlet.
  - **Reassignment check**: prior run `job-44614229` (res `0.4`) resolved `26,790` cells in `9` clusters with no doublet clusters and no reassignments. This cleaned run totals the same `26,790` cells in `8` clusters at res `0.3`; new C1 is marker-compatible with consolidation of the two previously annotated mature-plasma clusters (old C3 `5,595` + C4 `4,113` = `9,708` vs new C1 `9,770`), though exact cross-run membership is unverified, and all B-cell states are retained with no new doublet emerging.
  - **`k_param` Note**: `k_param = 100` at `res 0.3` resolves atypical-memory, activated follicular, co-stimulatory memory, naive-like, ER-stressed ASC, plasmablast, and GC B cells plus one merged mature-plasma cluster. A finer plasma subdivision may be explored at higher resolution if needed (unverified for this run; `LGALS1`/`RGS1` alone did not convincingly discriminate the old A/B split per the reviewed panels).
- **Master Index & References**: See `../../../NOTES_BARCODE_REASSIGNMENT_AND_CLEANUP.md` (Herrera-Uribe et al. 2021 `PMID: 34249103`; Wiarda et al. 2022 `PMID: 35995567`; Wiarda et al. 2026 `PMID: 41859083`).

---

## Summary Table (`harmony_snn_res.0.3`)

| Cluster | Cells (%) | Proposed Known Subset | QC / Doublet Flag | Key Named Top-25 Markers | Additional Canonical Markers (`markers_significant.csv`) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **1** | 9,770 (36.47%) | **Mature Plasma Cells (`JCHAIN^{hi} PRDM1⁺ TNFRSF17/BCMA⁺ SDC1⁺`)** | Keep (Clean) | `JCHAIN` (#1, `97%`), `IGLL1 (review)` (`50%`), `TNFRSF17` (`BCMA`, #10, `74%`), `TLCD4`, `SDC1` (`CD138`, #15, `41%`), `LGALS1`, `CCR2`, `FKBP11`, `PRDM1` (#22, `87%`) | `RGS1` (#24, `97%`), `MZB1` (#30, `71%`), `XBP1` (#56, `82%`) |
| **2** | 6,395 (23.87%) | **Atypical Memory (`FCRL5⁺ TBX21⁺`) / Marginal-Zone–like B Cells** | Keep (Clean) | `ZBTB32` (#1), `CD1.1 (review)` (#3, `48%`), `FCRL5` (#8, `74%`), `FCGR2B` (#10, `32%`), `CD72`, `ITGAD` (`CD11d`), `HCK`, `SCIMP`, `YBX3`, `SERPINB9` | `TBX21` (`T-bet`, #26, `39%`), `SPIB` (#27, `47%`), `TNFRSF13C` (`BAFF-R`, #36, `50%`), `FCRL3` (#49, `46%`), `CNR2` (#24, `53%`), `MS4A1` (#17, `93%`), `PAX5` (#147, `73%`), `CD19` (#211, `54%`) |
| **3** | 5,624 (20.99%) | **Activated Follicular / Pre-GC (`CCR7⁺ CD83⁺`, T–B Border) B Cells** | Keep (Clean) | `IRS2` (#2), `SATB1` (#3, `31%`), `ID3` (#4, `57.5%`), `CCR7` (#8, `74%`), `DUSP2` (#10), `MAFF`, `CD83` (#23, `82%`), `BCL6` (#21, `29%`), `ATF3` | `CXCR5` (#30, `44%`), `NOTCH2` (#41, `50%`), `FOXO1` (#43, `46%`), `CD40` (#121, `52%`), `BACH2` (#66, `47%`), `CD79B` (#100, `92%`) |
| **4** | 2,540 (9.48%) | **Co-Stimulatory Memory (`CD80⁺ CD86⁺`) / B1-like (`BHLHE41⁺ SPN/CD43⁺`) B Cells** | Keep (Clean) | `TCF7` (#1, `35%`), `FGD2` (#4), `BHLHE41` (#5, `51%`), `CD80` (#7, `49%`), `FCRLA` (#9), `SPN` (`CD43`, #10, `43%`), `TNFRSF1B`, `LTB`, `DUSP4`, `AHR`, `SEMA7A` | `CD86` (#37, `39%`), `IL2RB` (#40, `42%`), `ITGB1` (#42, `73%`), `MS4A1` (#234, `90%`), `BANK1` (#230, `79%`) |
| **5** | 1,145 (4.27%) | **Resting / Quiescent Naive-like B Cells (`EBF1^{hi} BACH2⁺ PAX5⁺ AFF3⁺`)** | Keep (Clean) | `FHIT` (#1), `AFF3` (#3, `97%`), `KCNQ5` (#4, `73%`), `TMCC3`, `SMYD3`, `SOX5`, `TOX`, `ARHGAP24`, `BACH2` (#15, `56%`) | `EBF1` (#25, `88%`), `CCND3` (#30, `82%`), `PAX5` (#67, `77%`), `BANK1` (#272, `78%`), `FCRL5` (#176, `54%`) |
| **6** | 597 (2.23%) | **ER-Stress / UPR-High (`ERN1/IRE1α⁺ EIF2AK3/PERK⁺`) ASCs** | Keep (Clean) | `FNDC3B` (#1, `88%`), `SCFD2`, `CADPS2`, `ABTB3`, `PITPNM2`, `GAB2`, `ABTB2`, `OSBPL10`, `ERN1` (`IRE1α`, #14, `75%`), `ENTPD1` (`CD39`) | `CREB3L2` (#30, `27%`), `SEC24A` (#31, `62%`), `EIF2AK3` (`PERK`, #42, `83%`), `IRF4` (#126, `73%`), `PRDM1` (#283, `68%`); `JCHAIN` detected (`89%`) without positive enrichment |
| **7** | 464 (1.73%) | **Proliferating / Newly Formed `IGHM⁺ BIRC5⁺` Plasmablasts** | Keep (Clean) | `TOM1L1`, `LAYN`, `IGHM` (#6, `88%`), `RFC3` (#7, `82%`), `BIRC5` (#9, `32.5%`), `AQP3` (#12, `97%`), `IL6R` (#14, `83%`), `TXNDC5` (#10, `100%`), `PDIA5` (#15, `69%`), `CRELD2` | `DNAJC3` (#37, `92%`), `HSP90B1` (#39, `100%`), `SLAMF1` (#40, `54%`), `MZB1` (#56, `88%`), `XBP1` (#244, `96%`) |
| **8** | 255 (0.95%) | **Germinal Center (GC) B Cells (`AICDA⁺ RGS13^{hi} S1PR2⁺ GCSAM⁺ BCL6⁺`; localization provisional)** | Keep (Clean) | `RGS13` (#1, `55%`), `AICDA` (`AID`, #2, `64%`), `NUGGC` (#3, `86%`), `S1PR2` (#4, `78%`), `NUAK1`, `BCL7A`, `MEF2B`, `GCSAM` (#16, `77%`), `NEIL1` (#22, `61%`), `FANCI` (#23, `28%`) | `ELL3` (#26, `90%`), `P2RY8` (#28, `45%`), `BCL6` (#29, `77%`), `PAX5` (#281, `88%`), `MS4A1` (#483, `92%`) |

---

## Concise Cluster Rationale

* **Mature & GC B-Cell Branch (representative `MS4A1⁺ PAX5⁺`, `15,959` cells; `MS4A1`/`PAX5` rows absent in some member clusters)**:
  * **C2 (`6,395`)**: Atypical memory (`FCRL5⁺ TBX21/T-bet⁺ ZBTB32⁺ HCK⁺ FCGR2B⁺ CD72⁺`) and marginal-zone–like (`CD1.1⁺ CNR2⁺`) B cells.
  * **C3 (`5,624`)**: Activated follicular / pre-GC-like B cells (provisional; `IRS2⁺ SATB1⁺ ID3⁺ CCR7^{hi} CD83^{hi} CXCR5⁺ CD40⁺ BCL6⁺ DUSP2⁺ ATF3⁺`). Tissue migration and stimulus are not demonstrated by these CSVs.
  * **C4 (`2,540`)**: Co-stimulatory / B1-like B cells (provisional; `BHLHE41⁺ SPN/CD43⁺ AHR⁺ TCF7⁺ CD80⁺ CD86⁺ TNFRSF1B⁺ LTB⁺ FCRLA⁺`). B1 identity and antigen-presenting function are not demonstrated by these CSVs.
  * **C5 (`1,145`)**: Naive-like B cells (provisional; core factors `EBF1⁺ BACH2⁺ PAX5⁺ AFF3^{hi}`; `CD83` and `CD80` remain detected, so globally low activation is not established).
  * **C8 (`255`)**: Germinal Center B cells expressing somatic hypermutation/CSR machinery (`AICDA`, `NUGGC`, `NEIL1`, `FANCI`) and GC-associated genes (`RGS13`, `S1PR2`, `GCSAM`, `P2RY8`, `MEF2B`, `BCL6`; tissue localization not established by these CSVs).
* **Antibody-Secreting Cell (ASC / Plasma) Branch (representative `JCHAIN^{hi} PRDM1⁺ TNFRSF17/BCMA⁺ SDC1/CD138⁺ XBP1⁺ MZB1⁺`, `10,831` cells; `TNFRSF17`/`SDC1` rows absent in some member clusters)**:
  * **C1 (`9,770`)**: Single mature-plasma cluster (`IGLL1 (review)⁺ CCR2⁺ FKBP11⁺ LGALS1⁺ RGS1^{hi}`) marker-compatible with consolidation of the old res-`0.4` plasma A/B programs; finer subdivision unverified for this run.
  * **C6 (`597`)**: ER-stress / UPR-associated ASCs (`ERN1/IRE1α⁺ EIF2AK3/PERK⁺ CREB3L2⁺`, ER-expansion factor `FNDC3B`).
  * **C7 (`464`)**: Plasmablast-like / cycling-enriched (`IGHM^{hi} IL6R^{hi} AQP3^{hi} BIRC5⁺ RFC3⁺ TXNDC5⁺ PDIA5⁺`); higher observed relative abundance in pooled N10/P10 samples without establishing absolute expansion or replicate-level treatment effects.

> Note: pre-red-palette outputs (viridis dotplot + prior CSVs/heatmaps) are archived in `old/` in this directory; current plots use the red scale.
