# B_cells Subclustering Identity Report (`job-44614229`, Resolution `0.4`)

- **Parent Run / Subset**: `job-44542861` $\rightarrow$ `B_cells` (`job-44614229`; Parent Clusters 7 [`16,010` B cells] + 9 [`10,780` Plasma cells])
- **Clustering Parameters**: `harmony_snn_res.0.4` (`26,790` cells across `9` clusters; `2,000` HVGs, `dims 1:20`, `k_param = 100`)
- **QC & Parameter Assessment**: **`0` doublet clusters, `0` reassignments**. Clustree has a 100% stable plateau from `0.4` to `0.5` (`in_prop = 1.00`). **Keep as-is (no re-subclustering needed).**
- **Master Index & References**: See `../../NOTES_BARCODE_REASSIGNMENT_AND_CLEANUP.md` (Herrera-Uribe et al. 2021 `PMID: 34249103`; Wiarda et al. 2022 `PMID: 35995567`; Wiarda et al. 2026 `PMID: 41859083`).

---

## Summary Table (`harmony_snn_res.0.4`)

| Cluster | Cells (%) | Proposed Known Subset | QC / Doublet Flag | Key Named Top-25 Markers | Additional Canonical Markers (`markers_significant.csv`) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **1** | 6,384 (23.83%) | **Atypical Memory (`FCRL5⁺ TBX21⁺`) / Marginal-Zone–like B Cells** | Keep (Clean) | `ZBTB32`, `FCRL5` (`74%`), `ITGAD` (`CD11d`), `FCGR2B`, `CD72`, `HCK`, `SCIMP`, `PTPN22`, `CD1.1`, `CNR2`, `MS4A1` | `TBX21` (`T-bet`, #26, `39%`), `SPIB` (#27, `47%`), `TNFRSF13C` (`BAFF-R`, #36, `50%`), `FCRL3` (#49, `46%`), `PAX5` (`73%`), `CD19` (`54%`) |
| **2** | 5,631 (21.02%) | **Activated Follicular / Pre-GC (`CCR7⁺ CD83⁺`, T–B Border) B Cells** | Keep (Clean) | `CCR7` (`74%`), `CD83` (`82%`), `BCL6`, `ID3`, `SATB1`, `DUSP2`, `ATF3`, `MAFF`, `IRS2`, `MYBL2`, `ISG15` | `CXCR5` (#31, `44%`), `NOTCH2` (#41, `50%`), `FOXO1` (#43, `46%`), `BACH2` (`47%`), `CD79B` (`92%`), `CD40` (`52%`) |
| **3** | 5,595 (20.88%) | **Mature Plasma Cells — Subset A (`IGLL1^{hi} JCHAIN^{hi} PRDM1⁺ SDC1⁺`)** | Keep (Clean) | `IGLL1` (`80%`), `JCHAIN` (`99%`), `PRDM1` (`Blimp-1`, `90%`), `TNFRSF17` (`BCMA`, `76%`), `SDC1` (`CD138`, `44%`), `FKBP11`, `TLCD4`, `CCR2`, `BST1` | `RGS1` (#27, `99%`), `AQP3` (#28, `78%`), `MZB1` (#31, `74%`), `SEC11C` (#35, `88%`), `BHLHA15` (`MIST1`, #38, `57%`), `XBP1` (#49, `84%`), `IGHM` (`87%`) |
| **4** | 4,113 (15.35%) | **Mature Plasma Cells — Subset B (`LGALS1⁺ RGS1^{hi} JCHAIN^{hi} PRDM1⁺`)** | Keep (Clean) | `JCHAIN` (`95%`), `PRDM1` (`85%`), `TNFRSF17` (`72%`), `SDC1` (`38%`), `FKBP11`, `SEC11C`, `SDF2L1`, `LGALS1` (`69%`), `RGS1` (`96%`), `HSP70.2`, `IER3` | `MZB1` (#26, `67%`), `CD69` (#28, `93%`), `HSP90B1` (#32, `98%`), `XBP1` (`79%`), `IGHM` (`77%`) |
| **5** | 2,552 (9.53%) | **Co-Stimulatory Memory (`CD80⁺ CD86⁺`) / B1-like (`BHLHE41⁺ SPN/CD43⁺`) B Cells** | Keep (Clean) | `BHLHE41` (`51%`), `SPN` (`CD43`, `43%`), `CD80` (`49%`), `TNFRSF1B`, `LTB`, `TCF7` (`35%`), `FCRLA`, `FGD2`, `DUSP4`, `AHR`, `SEMA7A`, `SEMA4D` | `CD86` (#37, `39%`), `IL2RB` (#40, `42%`), `ITGB1` (#43, `73%`), `MS4A1` (`90%`), `BANK1` (`79%`) |
| **6** | 1,145 (4.27%) | **Resting / Quiescent Naive-like B Cells (`EBF1^{hi} BACH2⁺ PAX5⁺`)** | Keep (Clean) | `EBF1` (`88%`), `BACH2` (`56%`), `AFF3` (`97%`), `SOX5`, `FHIT`, `KCNQ5`, `ARHGAP24`, `PREX1`, `TOX` | `CCND3` (#30, `82%`), `PAX5` (`77%`), `BANK1` (`78%`), `FCRL5` (`54%`) |
| **7** | 625 (2.33%) | **ER-Stress / UPR-High (`ERN1/IRE1α⁺ EIF2AK3/PERK⁺`) ASCs** | Keep (Clean) | `ERN1` (`IRE1α`, `75%`), `FNDC3B` (`87%`), `SCFD2`, `CADPS2`, `HDLBP`, `OSBPL10`, `ENTPD1` (`CD39`), `GAB2` | `SEC24A` (#28, `62%`), `CREB3L2` (#33, `27%`), `EIF2AK3` (`PERK`, #44, `82%`), `JCHAIN` (`89%`), `PRDM1` (`68%`), `IRF4` (`73%`) |
| **8** | 490 (1.83%) | **Proliferating / Newly Formed `IGHM⁺ BIRC5⁺` Plasmablasts** | Keep (Clean) | `IGHM` (`88%`), `BIRC5` (`31%`), `RFC3` (`81%`), `DUT` (`78%`), `H1-3`, `IL6R` (`82%`), `AQP3` (`97%`), `TXNDC5` (`100%`), `PDIA4` (`99%`), `PDIA5`, `CRELD2`, `SDF2L1` | `DNAJC3` (#34, `91%`), `HSP90B1` (#39, `100%`), `SLAMF1` (#40, `52%`), `XBP1` (`95%`), `MZB1` (`88%`), `JCHAIN` (`99%`) |
| **9** | 255 (0.95%) | **Germinal Center (GC) B Cells (`AICDA⁺ RGS13^{hi} S1PR2⁺ GCSAM⁺ BCL6⁺`)** | Keep (Clean) | `RGS13`, `AICDA` (`AID`), `NUGGC`, `S1PR2`, `GCSAM` (`HGAL`), `MEF2B`, `BCL7A`, `MYBL1`, `NEIL1`, `FANCI`, `RRM2B` | `ELL3` (#26, `90%`), `P2RY8` (#28, `45%`), `BCL6` (#29, `77%`), `BACH2` (`84%`), `PAX5` (`88%`), `MS4A1` (`92%`) |

---

## Concise Cluster Rationale

* **Mature & GC B-Cell Branch (`MS4A1⁺ PAX5⁺ CD19⁺`, `15,967` cells)**:
  * **C1 (`6,384`)**: Atypical memory (`FCRL5⁺ TBX21/T-bet⁺ ZBTB32⁺ ITGAD⁺ HCK⁺ FCGR2B⁺ CD72⁺`) and marginal-zone–like (`CD1.1⁺ CNR2⁺`) B cells.
  * **C2 (`5,631`)**: Antigen/cytokine-activated follicular B cells migrating to the T–B border (`CCR7^{hi} CD83^{hi} CXCR5⁺ CD40⁺ BCL6⁺ ID3⁺ DUSP2⁺ ATF3⁺`).
  * **C5 (`2,552`)**: Innate-like B1 (`BHLHE41⁺ SPN/CD43⁺ AHR⁺`) and APC-competent co-stimulatory memory B cells (`CD80⁺ CD86⁺ TNFRSF1B⁺ LTB⁺ FCRLA⁺`).
  * **C6 (`1,145`)**: Quiescent naive-like B cells maintaining core B-lineage transcription factors (`EBF1^{hi} BACH2⁺ PAX5⁺ AFF3^{hi}`) with low activation markers.
  * **C9 (`255`)**: BALT Germinal Center B cells expressing somatic hypermutation/CSR machinery (`AICDA`, `NUGGC`, `NEIL1`, `FANCI`) and GC confinement receptors (`RGS13`, `S1PR2`, `GCSAM`, `P2RY8`, `MEF2B`, `BCL6`).
* **Antibody-Secreting Cell (ASC / Plasma) Branch (`JCHAIN^{hi} PRDM1⁺ TNFRSF17/BCMA⁺ SDC1/CD138⁺ XBP1⁺ MZB1⁺`, `10,823` cells)**:
  * **C3 (`5,595`) & C4 (`4,113`)**: Terminally differentiated mature plasma cells distinguished by `IGLL1^{hi} CCR2⁺ BST1⁺` (C3) vs. `LGALS1⁺ RGS1^{hi} CD69^{hi}` (C4) and distinct Ig isotype/V-gene usage.
  * **C7 (`625`)**: Secretory-stressed ASCs with peak unfolded protein response (`ERN1/IRE1α⁺ EIF2AK3/PERK⁺ CREB3L2⁺`) and ER-expansion factor `FNDC3B`.
  * **C8 (`490`)**: Cycling `IGHM^{hi} IL6R^{hi} AQP3^{hi}` plasmablasts (`BIRC5⁺ RFC3⁺ DUT⁺ H1-3⁺`) that expand during 10h EVLP.
