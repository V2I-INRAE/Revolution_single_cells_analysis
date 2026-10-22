# Endothelial Subclustering Identity Report (`job-44831031`, Resolution `0.3`)

- **Parent Run / Subset**: `job-44542861` $\rightarrow$ `Endothelial` (`job-44831031`; cleaned baseline; Parent Clusters 12 [`7,986` Blood EC] + 13 [`6,168` Lymphatic EC])
- **Clustering Parameters**: `harmony_snn_res.0.3` (`13,836` cells across `5` clusters; `2,000` HVGs, `dims 1:20`, `k_param = 100`)
- **QC & Parameter Assessment**:
  - **Keep all `13,836` clean endothelial cells (Clusters 1–5)** spanning the five major pulmonary vascular/lymphatic EC classes.
  - **No doublet cluster.** Leukocyte markers `PTPRC/CD45`, `LCP1`, `CD53`, `CXCR4`, `RGS1` have no rows in `markers_significant.csv` in any cluster, and no cluster shows a coherent leukocyte–endothelial doublet signature in the exported markers; no further cluster-level removal is indicated. Residual individual doublets are not excluded by this analysis.
  - **Reassignment check**: prior run `job-44614231` resolved `14,154` cells in `6` clusters including C6 (`318` cells, leukocyte–EC doublets, removed per regrouping proposal). This cleaned run totals `13,836 = 14,154 − 318` cells in `5` clusters with no new doublet emerging; the removal worked.
  - **`k_param` Note**: `k_param = 100` at `res 0.3` resolves LEC, Venous, `gCap`, `aCap` (Aerocytes), and Arterial ECs. If subdividing Resting Venous (`ACKR1⁺ HDAC9⁺ SELE⁻`) from EVLP-Activated Venous (`ACKR1⁺ SELE⁺ VCAM1⁺`) and Systemic Bronchial Venous (`COL15A1⁺ VWA1⁺ PLVAP⁺`) ECs is desired, consider testing lower **`k_param` (`30`)** with the same `dims = 1:20` baseline for comparison (untested proposal, not a data-derived optimum).
- **Master Index & References**: See `../../../NOTES_BARCODE_REASSIGNMENT_AND_CLEANUP.md` (Gillich et al. 2020 `PMID: 33057196`; Schupp et al. 2021 `PMID: 34030460`).

---

## Summary Table (`harmony_snn_res.0.3`)

| Cluster | Cells (%) | Proposed Known Subset | QC / Doublet Flag | Key Named Top-25 Markers | Additional Canonical Markers (`markers_significant.csv`) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **1** | 6,154 (44.48%) | **Lymphatic Endothelial Cells (`LEC`, `CCL21^{hi} MMRN1^{hi} PROX1⁺ PDPN^{hi} FLT4⁺ LYVE1⁺`)** | Keep (Clean) | `SCG3`, `CCL21` (`84%`), `HOXD9`, `SEMA3D` (`79%`), `TNFAIP8L3`, `PTX3`, `CCDC3`, `TRPC6`, `PROX1` (`53%`), `MMRN1` (`94%`), `PDPN` (`88%`), `RELN` (`79%`), `NTN1`, `HS3ST1` | `FLT4` (`VEGFR3`, #40, `68%`), `TBX1` (#36, `29%`), `ACKR2` (#58, `64%`), `NR2F2` (#155, `80%`), `LYVE1` (#195, `76%`; detected at `47%` in Venous C2 without positive enrichment) |
| **2** | 2,634 (19.04%) | **Activated Venous / Post-Capillary Venular EC (`ACKR1^{hi} VCAM1^{hi} SELE⁺ SELP⁺ VWF^{hi}`)** | Keep (Clean) | `ACKR1` (`90%`), `IL6` (`56%`), `MATN2`, `OLR1`, `LHX6`, `CCL2`, `ADAMTS4`, `SELE` (`59%`), `F3`, `VCAM1` (`84%`), `SELP` (`65%`), `CXCL9` | `CXCL10` (#29, `80%`), `SERPINE1` (`PAI-1`, `92%`), `VWF` (#52, `94%`), `EDN1` (`89%`), `PLAT` (`93%`), `PLVAP` (`39%`), `TEK` (#152, `88%`) |
| **3** | 1,919 (13.87%) | **Alveolar General Capillary EC (`gCap` / CAP1, `GPIHBP1^{hi} APLNR^{hi} HPGD^{hi} MCAM⁺ KIT⁺`)** | Keep (Clean) | `GPIHBP1` (#1, `32%`), `ANK3`, `HPGD` (#3, `28%`), `APLNR` (#4, `26%`), `CLEC1B` (`56%`), `TNC`, `CD93` (`71%`), `MCAM` (`CD146`, `64%`), `CADM1` | `KIT` (#32, `40%`), `ITGA1` (#34, `81%`), `EPAS1` (#41, `94%`), `CLIC5` (#48, `84%`), `TMEM100` (#57, `52%`), `SPARCL1` (#74, `65%`) |
| **4** | 1,746 (12.62%) | **Alveolar Aerocytes / Gas-Exchange Capillary EC (`aCap` / CAP2, `APLN^{hi} EDNRB^{hi} PRX⁺ DPP4⁺ KITLG⁺ CA4⁺`)** | Keep (Clean) | `APLN` (#1, `29%`), `TRNP1` (`74%`), `ITGA4` (`50%`), `EDNRB` (`99%`), `MFGE8` (`56%`), `DPP4` (`58%`), `PDGFD`, `ITGA3`, `PRX` (`52%`), `GJC1` | `KITLG` (#26, `48%`), `CA4` (#30, `28%`), `RGCC` (#31, `90%`), `FOXF1` (#62, `57%`), `KDR` (`VEGFR2`, #73, `72%`), `AQP1` (#190, `84%`); **low `VWF` (`5%`) & `TEK` (`8%`) detection** |
| **5** | 1,383 (10.00%) | **Pulmonary Arterial EC (`GJA5/Cx40^{hi} BMX^{hi} HEY2⁺ SEMA3G^{hi} CXCL12⁺ SOX17⁺`)** | Keep (Clean) | `SLC52A3` (#1, `26%`), `MFAP5` (#2, `73%`), `HEY2` (#4, `32%`), `BMX` (`77%`), `FN1`, `GJA5` (`Cx40`, `86%`), `CXCL12`, `SSTR1`, `PDE3A`, `SEMA3G` (`94%`) | `HTR2B` (#30, `46%`), `PTGIS` (#36, `84%`), `FBLN5` (#45), `HEY1` (#48, `46%`), `IL33` (#49, `85%`), `SOX17` (#70, `69%`), `EFNB2` (#87, `88%`), `NOS3` (#115, `65%`) |

---

## Concise Cluster Rationale

* **C1 (`6,154`) — `LEC`**: Expresses `PROX1`, `CCL21` (`84%`), `MMRN1` (`94%`), `PDPN` (`88%`), `RELN`, `FLT4 (VEGFR3)`, `TBX1`, and `LYVE1` (`76%`; Venous C2 also expresses moderate `LYVE1` at `47%`).
* **C2 (`2,634`) — Venous EC**: Defined by pan-venular `ACKR1 (DARC)` (`90%`), `VWF` (`94%`), and leukocyte adhesion/inflammatory molecules (`VCAM1`, `SELE`, `SELP`, `IL6`, `CCL2`, `CXCL9/10`, `SERPINE1`). `SELE`/`SELP`/`VCAM1` are shared with Arterial C5 and were kept out of the discriminative core per the reviewed panel.
* **C3 (`1,919`) & C4 (`1,746`) — Alveolar Capillary Dichotomy (`gCap` vs. `aCap`)**: Marker programs consistent with the established paracrine pair (Gillich et al. 2020 `PMID: 33057196`): **C3 (`gCap`/CAP1)** expresses `GPIHBP1`, `HPGD`, `MCAM`, and receptors `APLNR` + `KIT`, while **C4 (`aCap`/Aerocyte CAP2)** expresses ligands `APLN` + `KITLG`, plus `EDNRB` (`99%`), `PRX`, `DPP4`, `CA4`, with low `VWF` (`5%`) and `TEK` (`8%`) detection.
* **C5 (`1,383`) — Arterial EC**: Expresses Notch targets `HEY2`/`HEY1`, `BMX` (`77%`), `GJA5 (Cx40)` (`86%`), `SEMA3G` (`94%`), `SOX17`, `EFNB2`, `CXCL12`, `PTGIS`, `NOS3`, and elastic-lamina ECM (`MFAP5`, `FBLN5`, `FN1`).
* **Doublet check (`REMOVE` none)**: No leukocyte–EC doublet cluster in this cleaned run; the prior C6 doublet population is gone and no coherent replacement leukocyte program emerged in the exported markers.

> Note: pre-red-palette outputs (viridis dotplot + prior CSVs/heatmaps) are archived in `old/` in this directory; current plots use the red scale.
