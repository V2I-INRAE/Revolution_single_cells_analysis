# Granulocytes Subclustering Identity Report (`job-44831033`, Resolution `0.3`)

- **Parent Run / Subset**: `job-44542861` $\rightarrow$ `Granulocytes` (`job-44831033`; Parent Clusters 6 [`18,411` Neutrophils] + 15 [`4,626` Eosinophils/Mast/Basophils])
- **Clustering Parameters**: `harmony_snn_res.0.3` (`23,037` cells across `6` clusters; `500` HVGs, `dims 1:10`, `k_param = 100`)
- **QC & Parameter Assessment**: **`0` doublet-dominated clusters identified from this marker review**. The cleaned input reproduces the first pass exactly (`23,037` cells, identical counts): 3 Neutrophil states (`17,950` cells), Eosinophils (`3,814`), Resident Pulmonary Mast Cells (`1,002`), and Basophils (`271`); the specified foreign-lineage screening markers (`CD3D/E`, `NCR1`, `CD19/MS4A1`, `CD68`, `PECAM1`, `EPCAM`, `JCHAIN`, `SDC1`) have no significant rows in the screened clusters. **Keep as-is (no re-subclustering needed; C4 Mast Cells can optionally be promoted to a top-level `Mast_cells` label; consider whether C5 immature/band neutrophils merit a distinct state).**
- **Master Index & References**: See `../../../NOTES_BARCODE_REASSIGNMENT_AND_CLEANUP.md` (Travaglini et al. 2020 `PMID: 33208946`; Gouin et al. 2023 `PMID: 37465668`).

---

## Summary Table (`harmony_snn_res.0.3`)

| Cluster | Cells (%) | Proposed Known Subset | QC / Doublet Flag | Key Named Top-25 Markers | Selected Identity Markers (`markers_significant.csv`) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **1** | 11,228 (48.74%) | **Inflammatory / Activated Neutrophils (`ICAM1⁺ OLR1/LOX-1⁺ IL1A⁺ CXCL2⁺ CXCL8^{hi} NFKB2⁺ CSF1⁺ EDN1⁺`)** | Keep (Clean) | `CSF1`, `EDN1`, `IL1A`, `ICAM1` (`75.2%`), `CXCL2` (`67.2%`), `NFKB2`, `OLR1` (`LOX-1`, #17, `31.8%` vs `13.3%`) | `CXCL8` (#48, `92.0%`), `RETN` (#35, `71.6%`) |
| **2** | 6,004 (26.06%) | **Mature-Like ISG-High Neutrophils (`SELL^{hi} FCGR3A^{hi} CSF3R^{hi} IFIT1/2/3⁺ ISG15⁺`)** | Keep (Clean) | `IFIT2`, `IFIT3`, `VNN1`, `IFIT1`, `FCGR3A`, `HERC5`, `ISG15`, `BCL6`, `SIGLEC5`, `SELL` (`CD62L`, #18, `79.0%`), `OAS1`, `GPCPD1`, `IL6R` | `CSF3R` (#74, `93.5%`) |
| **3** | 3,814 (16.56%) | **Eosinophils (`EPX⁺ CCR3⁺ ZBTB16^{hi} ALOX5⁺ ADGRG3⁺`)** | Keep (Clean) | `EPX` (`30.9%` vs `0.4%`, `log2FC = 6.85`), `ZBTB16` (`79.0%`), `CCR3` (#13, `31.9%`, `log2FC = 4.05`), `STAP1`, `PAFAH2` | `ALOX5` (#26, `69.0%`), `ADGRG3` (`GPR97`, #38, `30.3%`), `CXCR4` (#67, `93.7%`) |
| **4** | 1,002 (4.35%) | **Resident Pulmonary Mast Cells (`KIT⁺ MS4A2⁺ FCER1A⁺ GATA2⁺ PTGS1⁺ HPGD⁺ RAB27B⁺ MLPH⁺ VWA5A⁺`)** | Keep (Clean; Optional promote to `Mast_cells`) | `MLPH` (`log2FC = 10.20`), `CCN3` (`8.47`), `RAB27B` (`40.1%`), `PDE3A`, `PTGS1` (`31.4%`), `HPGD` (`36.0%`), `KIT` (`CD117`, #8, `46.9%`, `log2FC = 6.09`), `MS4A2` (#15, `38.4%`), `GATA2` (`25.2%`), `FCER1A` (`30.3%`) | `HRH4` (#32, `28.9%`), `NKG7` (#25, `31.4%`; NK/T-associated transcript, alone not evidence of contamination) |
| **5** | 718 (3.12%) | **Immature / Band & Secondary-Granule–Rich Neutrophils (`LTF^{hi} LCN2^{hi} MMP8^{hi} MMP9⁺ S100A8/9/12^{hi} FCN1⁺ CD24^{hi}`)** | Keep (Clean) | `LTF` (`30.1%` vs `0.9%`, `log2FC = 8.15`), `S100A12` (`73.1%`), `LCN2` (`37.0%`), `MMP8` (`28.3%`), `S100A8` (`79.7%`), `MMP9` (`38.2%`), `FCN1` (`43.9%`), `S100A9` (`69.8%`) | `CSF3R` (#244, `89.7%`), `CD24` (#32, `83.7%`) |
| **6** | 271 (1.18%) | **Basophils (`FCER1A^{hi} MS4A2^{hi} GATA2^{hi} ITGA2/CD49b⁺ GCSAML⁺ RHEX⁺ CA8⁺ HRH4^{hi} MMP9^{hi} RUNX2⁺`)** | Keep (Clean) | `TMPRSS11E` (`60.5%`), `GPR82`, `MMP9` (`73.1%`), `FCER1A` (`65.7%`), `GCSAML` (`34.7%`), `GATA2` (#6, `50.2%`), `CA8` (`46.1%`), `MS4A2` (`71.2%`), `RHEX` (`43.5%`), `HRH4` (`55.4%`), `ITGA2` (`CD49b`, `44.3%`) | `RUNX2` (#32, `59.4%`); `KIT`: no significant row (not scored as negative); `CSF3R`-depleted (#496, `log2FC = −1.52`, `25.5%` vs `72.5%`) |

---

## Concise Cluster Rationale

* **Neutrophil States (`C1` = `11,228`, `C2` = `6,004`, `C5` = `718`)**:
  * **C1**: Inflammatory/activated neutrophil program upregulating `ICAM1` (`75%`), `OLR1 (LOX-1)` (#17), `IL1A`, `CXCL2`, `CXCL8` (`92%`), `NFKB2`, `CSF1`, `EDN1`, and `NLRP3` (working label: tissue-associated pro-inflammatory state; anatomical/condition assignment needs evidence beyond these exploratory markers).
  * **C2**: Mature-like neutrophils (`SELL/CD62L^{hi}` #18 `79%`, `CSF3R^{hi}` #74 `94%`, `FCGR3A^{hi}`, `SIGLEC5⁺`) with an interferon-stimulated signature (`IFIT1/2/3`, `ISG15`, `OAS1/2`, `MX1`).
  * **C5**: Immature/band neutrophils actively transcribing secondary/specific (`LTF`, `LCN2`, `MMP8`) and tertiary (`MMP9`) granule genes plus `S100A8/9/12^{hi}`, `FCN1`, and `CD24^{hi}` (`CD24` is compatible with neutrophil differentiation and is not, by itself, evidence of dendritic-cell contamination).
* **Eosinophil / Mast / Basophil Lineages (`C3` = `3,814`, `C4` = `1,002`, `C6` = `271`)**:
  * **C3 (Eosinophils)**: Defined by `EPX` (`log2FC = 6.85`), `CCR3` (#13, `log2FC = 4.05`), `ZBTB16 (PLZF)` (`79%`), `ALOX5`, and `ADGRG3 (GPR97)` (#38; `CSF3R`-depleted at #1063, `5.1%` vs `85.2%`).
  * **C4 (Mast Cells)**: Defined by `KIT (CD117)` (#8, `47%`), `MS4A2`, `FCER1A`, `GATA2`, `PTGS1/2`, `HPGD`, `MLPH`, `RAB27B`, and `VWA5A` (note `NKG7` #25 transcript present, alone not evidence of contamination).
  * **C6 (Basophils)**: `KIT` not represented in the significant-marker table and `CSF3R`-depleted (#496) cells expressing maximal `FCER1A` (`66%`), `MS4A2` (`71%`), `GATA2` (#6, `50%`), basophil integrin `ITGA2 (CD49b/DX5)` (`44%`), `GCSAML`, `RHEX`, `CA8`, `HRH4` (`55%`), `MMP9` (`73%`), and `RUNX2` (`59%`).

> Note: pre-red-palette outputs (viridis dotplot + prior CSVs/heatmaps) are archived in `old/` in this directory; current plots use the red scale.
