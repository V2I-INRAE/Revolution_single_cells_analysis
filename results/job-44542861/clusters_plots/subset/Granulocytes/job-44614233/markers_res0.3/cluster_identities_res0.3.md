# Granulocytes Subclustering Identity Report (`job-44614233`, Resolution `0.3`)

- **Parent Run / Subset**: `job-44542861` $\rightarrow$ `Granulocytes` (`job-44614233`; Parent Clusters 6 [`18,411` Neutrophils] + 15 [`4,626` Eosinophils/Mast/Basophils])
- **Clustering Parameters**: `harmony_snn_res.0.3` (`23,037` cells across `6` clusters; `500` HVGs, `dims 1:10`, `k_param = 100`)
- **QC & Parameter Assessment**: **`0` doublet clusters**. Cleanly resolves 3 Neutrophil states (`17,950` cells), Eosinophils (`3,814`), Resident Pulmonary Mast Cells (`1,002`), and Basophils (`271`). **Keep as-is (no re-subclustering needed; C4 Mast Cells can optionally be promoted to a top-level `Mast_cells` label).**
- **Master Index & References**: See `../../NOTES_BARCODE_REASSIGNMENT_AND_CLEANUP.md` (Travaglini et al. 2020 `PMID: 33208946`; Gouin et al. 2023 `PMID: 37465668`).

---

## Summary Table (`harmony_snn_res.0.3`)

| Cluster | Cells (%) | Proposed Known Subset | QC / Doublet Flag | Key Named Top-25 Markers | Additional Canonical Markers (`markers_significant.csv`) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **1** | 11,228 (48.74%) | **EVLP-Activated Pro-Inflammatory Tissue Neutrophils (`ICAM1⁺ OLR1/LOX-1⁺ IL1A⁺ CXCL2⁺ CXCL8^{hi} NFKB2⁺ CSF1⁺ EDN1⁺`)** | Keep (Clean) | `CSF1`, `EDN1`, `IL1A`, `ICAM1` (`75.2%`), `CCL5`, `C3`, `CCL4`, `CXCL2` (`67.2%`), `NFKB2`, `SLAMF7`, `NOD2`, `OLR1` (`LOX-1`), `TRAF1`, `BATF3`, `WDFY1`, `NLRP3`, `RALGDS`, `SNX10` | `RETN` (#35, `71.6%`), `CD14` (#36, `60.0%`), `TLR4` (#37, `66.2%`), `GPR84` (#43, `72.9%`), `CXCL8` (#48, `92.0%`), `NFKBIA` (#150, `98.7%`), `IL1R2` (#178, `45.2%`) |
| **2** | 6,004 (26.06%) | **Mature Circulating / Marginated & ISG-High Neutrophils (`SELL^{hi} FCGR3A^{hi} CSF3R^{hi} IFIT1/2/3⁺ ISG15⁺`)** | Keep (Clean) | `IFIT2`, `IFIT3`, `VNN1`, `IFIT1`, `FCGR3A`, `DHRS9`, `STEAP4`, `SDS` (`94.1%`), `CHI3L1`, `HERC5`, `CFD`, `ISG15`, `BCL6`, `SIGLEC5`, `SELL` (`CD62L`, `79.0%`), `OAS1`, `GPCPD1`, `HERC6`, `IL6R` | `OAS2` (#27, `61.4%`), `MX1` (#38, `30.1%`), `IRF7` (#55, `31.2%`), `CSF3R` (#74, `93.5%`), `STAT1` (#120, `58.3%`), `CD24` (#194, `66.3%`) |
| **3** | 3,814 (16.56%) | **Eosinophils (`EPX⁺ CCR3⁺ ZBTB16^{hi} ALOX5⁺ ADGRG3⁺`)** | Keep (Clean) | `EPX` (`30.9%` vs `0.4%`, `log2FC = 6.85`), `ZBTB16` (`79.0%`), `FAR1`, `PPA2`, `IGF1R`, `CCR3` (`31.9%`, `log2FC = 4.05`), `OR51E1`, `CLNK`, `STAP1`, `ALDH3A2`, `GRK3`, `PAFAH2`, `PYGL` | `ALOX5` (#26, `69.0%`), `PGD` (#30, `71.7%`), `MCTP1` (#31, `78.2%`), `ADGRG3` (`GPR97`, #37, `30.3%`), `RGS1` (#40, `89.3%`), `CXCR4` (#67, `93.7%`), `CD63` (#737, `63.1%`); `CSF3R⁻` (`5.1%`) |
| **4** | 1,002 (4.35%) | **Resident Pulmonary Mast Cells (`KIT⁺ MS4A2⁺ FCER1A⁺ GATA2⁺ PTGS1⁺ HPGD⁺ RAB27B⁺ MLPH⁺ VWA5A⁺`)** | Keep (Clean; Optional promote to `Mast_cells`) | `MLPH` (`log2FC = 10.20`), `CCN3` (`8.47`), `RAB27B` (`40.1%`), `PDE3A`, `PTGS1` (`31.4%`), `HPGD` (`36.0%`), `KIT` (`CD117`, `46.9%`, `log2FC = 6.09`), `ABCC3`, `GAB1`, `LAT`, `VWA5A`, `MS4A2` (`38.4%`), `PTGS2`, `GATA2` (`25.2%`), `FCER1A` (`30.3%`) | `SOX13` (#27), `HRH4` (#32, `28.9%`), `HDAC9` (#35), `RUNX2` (#174, `27.1%`), `ALOX5` (#241, `32.9%`) |
| **5** | 718 (3.12%) | **Immature / Band & Secondary-Granule–Rich Neutrophils (`LTF^{hi} LCN2^{hi} MMP8^{hi} MMP9⁺ S100A8/9/12^{hi} FCN1⁺ CD24^{hi}`)** | Keep (Clean) | `LTF` (`30.1%` vs `0.9%`, `log2FC = 8.15`), `S100A12` (`73.1%`), `LCN2` (`37.0%`), `MMP8` (`28.3%`), `S100A8` (`79.7%`), `CD52`, `S100A11`, `MMP9` (`38.2%`), `FCN1` (`43.9%`), `S100A9` (`69.8%`), `ANXA1`, `SERPINB1` | `IL1R2` (#28, `65.5%`), `CD24` (#32, `83.7%`), `ALOX5AP` (#38, `78.1%`), `RETN` (#47, `85.1%`), `CD14` (#148, `70.9%`), `SELL` (#191, `79.0%`), `CSF3R` (#244, `89.7%`) |
| **6** | 271 (1.18%) | **Basophils (`FCER1A^{hi} MS4A2^{hi} GATA2^{hi} ITGA2/CD49b⁺ GCSAML⁺ RHEX⁺ CA8⁺ HRH4^{hi} MMP9^{hi} RUNX2⁺ KIT⁻`)** | Keep (Clean) | `TMPRSS11E` (`60.5%`), `GPR82`, `MMP9` (`73.1%`), `FCER1A` (`65.7%`), `GCSAML` (`34.7%`), `GATA2` (`50.2%`), `CA8` (`46.1%`), `MS4A2` (`71.2%`), `PLBD1`, `RHEX` (`43.5%`), `HRH4` (`55.4%`), `SYTL3`, `ITGA2` (`CD49b`, `44.3%`) | `RUNX2` (#32, `59.4%`), `IKZF2` (#34, `33.2%`), `ALOX5` (#72, `56.1%`), `IKZF1` (#212, `38.4%`), `CXCR4` (#316, `60.9%`); `KIT⁻` & `CSF3R⁻` |

---

## Concise Cluster Rationale

* **Neutrophil States (`C1` = `11,228`, `C2` = `6,004`, `C5` = `718`)**:
  * **C1**: EVLP-expanded tissue-infiltrating neutrophils upregulating `ICAM1` (`75%`), `OLR1 (LOX-1)`, `IL1A`, `CXCL2`, `CXCL8` (`92%`), `NFKB2`, `CSF1`, `EDN1`, and `NLRP3`.
  * **C2**: Mature intravascular/marginated neutrophils (`SELL/CD62L^{hi}` `79%`, `CSF3R^{hi}` `94%`, `FCGR3A^{hi}`, `SIGLEC5⁺`) with an interferon-stimulated signature (`IFIT1/2/3`, `ISG15`, `OAS1/2`, `MX1`).
  * **C5**: Immature/band neutrophils actively transcribing secondary/specific (`LTF`, `LCN2`, `MMP8`) and tertiary (`MMP9`) granule genes plus `S100A8/9/12^{hi}`, `FCN1`, and `CD24^{hi}`.
* **Non-Neutrophil $\text{GATA2}^+$ Lineages (`C3` = `3,814`, `C4` = `1,002`, `C6` = `271`)**:
  * **C3 (Eosinophils)**: Defined by `EPX` (`log2FC = 6.85`), `CCR3` (`log2FC = 4.05`), `ZBTB16 (PLZF)` (`79%`), `ALOX5`, and `ADGRG3 (GPR97)` (`CSF3R⁻`).
  * **C4 (Mast Cells)**: Defined by `KIT (CD117)` (`46.9%`), `MS4A2`, `FCER1A`, `GATA2`, `PTGS1/2`, `HPGD`, `MLPH`, `RAB27B`, and `VWA5A`.
  * **C6 (Basophils)**: `KIT⁻ CSF3R⁻` cells expressing maximal `FCER1A` (`66%`), `MS4A2` (`71%`), `GATA2` (`50%`), basophil integrin `ITGA2 (CD49b/DX5)` (`44%`), `GCSAML`, `RHEX`, `CA8`, `HRH4` (`55%`), `MMP9` (`73%`), and `RUNX2` (`59%`).
