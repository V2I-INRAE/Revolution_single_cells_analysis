# Endothelial Subclustering Identity Report (`job-44614231`, Resolution `0.3`)

- **Parent Run / Subset**: `job-44542861` $\rightarrow$ `Endothelial` (`job-44614231`; Parent Clusters 12 [`7,986` Blood EC] + 13 [`6,168` Lymphatic EC])
- **Clustering Parameters**: `harmony_snn_res.0.3` (`14,154` cells across `6` clusters; `2,000` HVGs, `dims 1:15`, `k_param = 100`)
- **QC & Parameter Assessment**:
  - **Keep `13,836` clean endothelial cells (Clusters 1–5)** spanning all 5 canonical pulmonary vascular/lymphatic beds.
  - **Remove Cluster 6 (`318` cells, `2.25%`, `PTPRC/CD45⁺ LCP1⁺ CD53⁺` Leukocyte–EC doublets)**.
  - **`k_param` Note**: `k_param = 100` at `res 0.3` already resolves LEC, Venous, `gCap`, `aCap` (Aerocytes), and Arterial ECs. If re-subclustering after removing C6 doublets to split Resting Venous (`ACKR1⁺ HDAC9⁺ SELE⁻`) from EVLP-Activated Venous (`ACKR1⁺ SELE⁺ VCAM1⁺`) and Systemic Bronchial Venous (`COL15A1⁺ VWA1⁺ PLVAP⁺`) ECs, lower **`k_param` to `30`** (`dims = 1:15`).
- **Master Index & References**: See `../../NOTES_BARCODE_REASSIGNMENT_AND_CLEANUP.md` (Gillich et al. 2020 `PMID: 33057196`; Schupp et al. 2021 `PMID: 34097468`).

---

## Summary Table (`harmony_snn_res.0.3`)

| Cluster | Cells (%) | Proposed Known Subset | QC / Doublet Flag | Key Named Top-25 Markers | Additional Canonical Markers (`markers_significant.csv`) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **1** | 6,172 (43.61%) | **Lymphatic Endothelial Cells (`LEC`, `CCL21^{hi} MMRN1^{hi} PROX1⁺ PDPN^{hi} FLT4⁺ LYVE1⁺`)** | Keep (Clean) | `CCL21` (`84%`), `MMRN1` (`93%`), `PROX1` (`53%`), `PDPN` (`88%`), `RELN` (`79%`), `SEMA3D`, `NTN1`, `PTX3`, `HOXD8`, `HOXD9`, `HS3ST1` | `TBX1` (#34), `FLT4` (`VEGFR3`, #37, `68%`), `ACKR2` (#57), `NR2F2` (#146, `80%`), `LYVE1` (#191, `76%` vs `47%` in Venous C2) |
| **2** | 2,647 (18.70%) | **Activated Venous / Post-Capillary Venular EC (`ACKR1^{hi} VCAM1^{hi} SELE⁺ SELP⁺ VWF^{hi}`)** | Keep (Clean) | `ACKR1` (`90%`), `VCAM1` (`83%`), `SELE` (`59%`), `SELP` (`65%`), `IL6`, `CCL2`, `CXCL9`, `OLR1`, `F3`, `F2R`, `ADGRG6`, `OMD`, `ADAMTS4` | `CXCL10` (#29, `80%`), `SERPINE1` (`PAI-1`, #32, `92%`), `VWF` (#50, `94%`), `EDN1` (#57, `89%`), `PLVAP` (#65, `39%`), `PLAT` (#105, `93%`) |
| **3** | 1,863 (13.16%) | **Alveolar General Capillary EC (`gCap` / CAP1, `GPIHBP1^{hi} APLNR^{hi} HPGD^{hi} MCAM⁺ KIT⁺`)** | Keep (Clean) | `GPIHBP1`, `APLNR`, `HPGD`, `MCAM` (`CD146`, `65%`), `CD93` (`72%`), `ANK3`, `TNC`, `CLEC1B`, `CLEC1A`, `CADM1` | `KIT` (#27, `41%`), `ITGA1` (#35, `81%`), `EPAS1` (#39, `94%`), `CLIC5` (#45, `85%`), `TMEM100` (#51, `53%`), `SPARCL1` (#70, `66%`) |
| **4** | 1,759 (12.43%) | **Alveolar Aerocytes / Gas-Exchange Capillary EC (`aCap` / CAP2, `APLN^{hi} EDNRB^{hi} PRX⁺ DPP4⁺ KITLG⁺ CA4⁺`)** | Keep (Clean) | `APLN`, `EDNRB` (`98.5%`), `PRX` (`52%`), `DPP4` (`58%`), `KITLG` (`48%`), `ITGA4`, `ITGA3`, `MFGE8`, `LMO7`, `PDGFD`, `GJC1` | `CA4` (#29, `28%`), `RGCC` (#32, `91%`), `FOXF1` (#60, `56%`), `KDR` (`VEGFR2`, #71, `71%`), `AQP1` (#182, `84%`); **`VWF⁻` (`5%`) & `TEK⁻` (`8%`)** |
| **5** | 1,395 (9.86%) | **Pulmonary Arterial EC (`GJA5/Cx40^{hi} BMX^{hi} HEY2⁺ SEMA3G^{hi} CXCL12⁺ SOX17⁺`)** | Keep (Clean) | `SLC52A3`, `HEY2`, `BMX` (`77%`), `GJA5` (`Cx40`, `85%`), `SEMA3G` (`93%`), `CXCL12`, `KCNE3`, `SSTR1`, `PDE3A`, `FBLN1`, `MFAP5`, `FN1`, `PI16` | `HTR2B` (#30), `PTGIS` (#36, `84%`), `FBLN5` (#44), `HEY1` (#48), `IL33` (#50, `84%`), `SOX17` (#69, `68%`), `EFNB2` (#86, `88%`), `NOS3` (#107) |
| **6** | 318 (2.25%) | **Leukocyte–Endothelial Doublets (`PTPRC/CD45⁺ LCP1⁺ CD53⁺ CXCR4⁺ RGS1⁺` + Pan-EC)** | **REMOVE (Doublet)** | `RGS1`, `CXCR4`, `PTPRC` (`CD45`, `60%`), `FCGR3A`, `CD69`, `CORO1A`, `CXCL8`, `IL2RG`, `IKZF3`, `PLEK`, `LCP2`, `SLA`, `LCP1`, `CD53` | Diluted pan-EC transcripts (`PECAM1` `77%`, `CLDN5` `48%`, `CDH5` `43%`) co-detected with pan-leukocyte markers |

---

## Concise Cluster Rationale

* **C1 (`6,172`) — `LEC`**: Expresses `PROX1`, `CCL21`, `MMRN1`, `PDPN`, `RELN`, `FLT4 (VEGFR3)`, `TBX1`, and `LYVE1` (`76%`; ranked #191 within ECs because Venous C2 also expresses moderate `LYVE1` at `47%`).
* **C2 (`2,647`) — Venous EC**: Defined by pan-venular `ACKR1 (DARC)` (`90%` vs `1.7%`), `VWF` (`94%`), and EVLP-induced leukocyte adhesion/inflammatory molecules (`VCAM1`, `SELE`, `SELP`, `IL6`, `CCL2`, `CXCL9/10`, `F3`, `SERPINE1`).
* **C3 (`1,863`) & C4 (`1,759`) — Alveolar Capillary Dichotomy (`gCap` vs. `aCap`)**: Textbook paracrine pair (Gillich et al. 2020 `PMID: 33057196`): **C3 (`gCap`/CAP1)** expresses `GPIHBP1`, `HPGD`, `MCAM`, and receptors `APLNR` + `KIT`, while **C4 (`aCap`/Aerocyte CAP2)** expresses ligands `APLN` + `KITLG`, plus `EDNRB` (`98.5%`), `PRX`, `DPP4`, `CA4`, and lacks `VWF` (`5%`) and `TEK` (`8%`).
* **C5 (`1,395`) — Arterial EC**: Expresses Notch targets `HEY2`/`HEY1`, `BMX` (`77%`), `GJA5 (Cx40)` (`85%`), `SEMA3G` (`93%`), `SOX17`, `EFNB2`, `CXCL12`, `PTGIS`, `NOS3`, and elastic-lamina ECM (`MFAP5`, `FBLN5`).
* **C6 (`318`) — Doublets (`REMOVE`)**: Leukocyte (`PTPRC`, `LCP1`, `CD53`, `CXCR4`, `RGS1`) + Endothelial (`PECAM1`, `CLDN5`, `CDH5`) doublets.
