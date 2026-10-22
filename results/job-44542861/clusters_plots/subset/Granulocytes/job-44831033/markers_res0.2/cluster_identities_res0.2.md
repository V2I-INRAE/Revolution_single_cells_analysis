# Granulocytes Subclustering Identity Report (`job-44831033`, Resolution `0.2`)

- **Parent Run / Subset**: `job-44542861` $\rightarrow$ `Granulocytes` (`job-44831033`; Parent Clusters 6 [`18,411` Neutrophils] + 15 [`4,626` Eosinophils/Mast/Basophils])
- **Clustering Parameters**: `harmony_snn_res.0.2` (`23,037` cells across `5` clusters; `500` HVGs, `dims 1:10`, `k_param = 100`)
- **QC & Parameter Assessment**: **`0` doublet-dominated clusters identified from this marker review.** The specified foreign-lineage screening markers (`CD3D/E`, `NCR1`, `CD19/MS4A1`, `CD68`, `PECAM1`, `EPCAM`, `JCHAIN`, `SDC1`) have no significant rows in any cluster. **Keep all `5` clusters as-is.**
  - **Reassignment check**: prior run at res `0.3` resolved the same `23,037` cells in `6` clusters. This res-`0.2` clustering is marker-compatible with merger of the two mature-neutrophil states and retention of the rest: new C1 (`17,253`) matches old inflammatory C1 (`11,228`) + old ISG-high C2 (`6,004`) = `17,232` vs `17,253`; new C2 (`3,809`) matches old eosinophil C3 (`3,814`); new C3 (`1,000`) matches old mast C4 (`1,002`); new C4 (`708`) matches old band-neutrophil C5 (`718`); new C5 (`267`) matches old basophil C6 (`271`). Cross-run membership is barcode-unverified; counts are compared, not traced.
- **Master Index & References**: See `../../../NOTES_BARCODE_REASSIGNMENT_AND_CLEANUP.md` (Travaglini et al. 2020 `PMID: 33208946`; Gouin et al. 2023 `PMID: 37465668`).

---

## Summary Table (`harmony_snn_res.0.2`)

| Cluster | Cells (%) | Proposed Known Subset | QC / Doublet Flag | Selected Identity Markers (`markers_significant.csv`) |
| :--- | :--- | :--- | :--- | :--- |
| **1** | 17,253 (74.89%) | **Mature Neutrophils, ISG-High / Inflammatory Mix (`IFIT2⁺ OAS1⁺ ISG15⁺ SELL⁺ CSF3R⁺ CXCL8⁺ ICAM1⁺`)** | Keep (Clean) | `IFIT2` (#1, `30.6%`), `OAS1` (#2, `27.9%`), `ISG15` (#3, `35.6%`), `OAS2` (#5, `55.8%`), `IL18` (#10, `30.6%`), `IFIT1` (#12, `43.9%`), `CSF1` (#23, `32.2%`), `CSF3R` (#48, `87.8%`), `SELL` (`CD62L`, #50, `66.7%`), `CXCL8` (#56, `88.0%`), `ICAM1` (#79, `61.9%`); `IFIT3`: no significant row |
| **2** | 3,809 (16.53%) | **Eosinophils (`EPX⁺ CCR3⁺ ZBTB16^{hi} ALOX5⁺ ADGRG3⁺ CXCR4⁺`)** | Keep (Clean) | `EPX` (#1, `31.0%` vs `0.4%`, `log2FC = 6.85`), `ZBTB16` (#4, `79.2%`), `CCR3` (#13, `31.9%`), `ALOX5` (#26, `69.1%`), `ADGRG3` (`GPR97`, #38, `30.3%`), `CXCR4` (#67, `93.8%`); `CSF3R`-depleted (#1065, `5.0%` vs `85.2%`) |
| **3** | 1,000 (4.34%) | **Resident Pulmonary Mast Cells (`KIT⁺ MS4A2⁺ FCER1A⁺ GATA2⁺ PTGS1⁺ HPGD⁺ RAB27B⁺ MLPH⁺`)** | Keep (Clean; optional promote to `Mast_cells`) | `MLPH` (#1, `log2FC = 10.21`), `CCN3` (#2, `8.48`), `RAB27B` (#3, `40.2%`), `KIT` (`CD117`, #8, `47.0%`), `MS4A2` (#15, `38.5%`), `GATA2` (#19, `25.3%`), `FCER1A` (#21, `30.4%`), `HRH4` (#32, `29.0%`), `RUNX2` (#174, `27.2%`) |
| **4** | 708 (3.07%) | **Immature / Band Neutrophils (`LTF^{hi} LCN2^{hi} MMP8^{hi} S100A8/9/12^{hi} FCN1⁺ CD24^{hi}`)** | Keep (Clean) | `LTF` (#1, `30.4%` vs `0.9%`, `log2FC = 8.17`), `S100A12` (#2, `73.0%`), `LCN2` (#3, `37.4%`), `MMP8` (#4, `28.5%`), `S100A8` (#5, `79.9%`), `MMP9` (#9, `38.4%`), `FCN1` (#10, `44.1%`), `S100A9` (#11, `69.6%`), `CD24` (#31, `83.8%`) |
| **5** | 267 (1.16%) | **Basophils (`FCER1A^{hi} MS4A2^{hi} GATA2^{hi} MMP9^{hi} HRH4⁺ ITGA2/CD49b⁺ RUNX2⁺`)** | Keep (Clean) | `TMPRSS11E` (#1, `61.4%`), `MMP9` (#3, `73.4%`), `FCER1A` (#4, `66.7%`), `GATA2` (#6, `50.9%`), `MS4A2` (#8, `71.5%`), `HRH4` (#13, `56.2%`), `ITGA2` (`CD49b`, `44.9%`), `RUNX2` (#33, `60.3%`); `CSF3R`-depleted (#504, `24.3%` vs `72.5%`); `KIT`: no significant row (not scored as negative) |

---

## Concise Cluster Rationale

* **Merged Mature Neutrophils (C1 = `17,253`)**: At res `0.2` the inflammatory and ISG-high programs co-cluster: interferon-stimulated genes (`IFIT2` #1, `OAS1` #2, `ISG15` #3, `OAS2` #5, `IFIT1` #12; `IFIT3` has no significant row) alongside inflammatory mediators (`CXCL8` #56 `88.0%`, `ICAM1` #79, `IL1A` #33, `CXCL2` #35, `NFKB2` #46, `OLR1` #44, `CSF1` #23, `EDN1` #61) on a mature-neutrophil backbone (`CSF3R` #48 `87.8%`, `SELL/CD62L` #50 `66.7%`, `FCGR3A` #63). Whether this merger groups two states or a continuum is not resolved at this resolution.
* **Eosinophils (C2 = `3,809`)**: `EPX` (`log2FC = 6.85`), `CCR3` (#13), `ZBTB16` (#4 `79.2%`), `ALOX5` (#26), `ADGRG3/GPR97` (#38), `CXCR4` (#67 `93.8%`), with `CSF3R` depletion (#1065).
* **Mast Cells (C3 = `1,000`)**: `KIT` (#8), `MS4A2`, `FCER1A`, `GATA2`, `PTGS1`, `HPGD`, `RAB27B`, `MLPH`, `HRH4`; `RUNX2` present at #174 (lower rank than in basophils).
* **Band Neutrophils (C4 = `708`)**: Secondary/specific (`LTF`, `LCN2`, `MMP8`) and tertiary (`MMP9`) granule-gene transcription plus `S100A8/9/12^{hi}`, `FCN1`, and `CD24^{hi}` (#31; compatible with neutrophil differentiation, not by itself evidence of DC contamination).
* **Basophils (C5 = `267`)**: Maximal `FCER1A`/`MS4A2`/`GATA2`, `MMP9` (#3 `73.4%`), `HRH4` (#13), `ITGA2/CD49b`, `RUNX2` (#33), `GCSAML`, `RHEX`, `CA8`; `KIT` without a significant row and `CSF3R`-depleted (#504).

> Note: pre-red-palette outputs (viridis dotplot + prior CSVs/heatmaps) are archived in `old/` in this directory; current plots use the red scale.
