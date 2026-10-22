# Epithelial Subclustering Identity Report (`job-44831032`, Resolution `0.3`)

- **Parent Run / Subset**: `job-44542861` $\rightarrow$ `Epithelial` (`job-44831032`; cleaned baseline; Parent Clusters 17, 18, 22, 23, 24)
- **Clustering Parameters**: `harmony_snn_res.0.3` (`9,426` cells across `7` clusters; `2,000` HVGs, `dims 1:20`, `k_param = 100`)
- **QC & Parameter Assessment**:
  - **Keep all `9,426` cells (Clusters 1–7).** No new cluster reproduces the previously removed T/NK–epithelial marker program (`CD2`/`CD3D`/`CD3E`/`GZMA`/`GZMB`/`CCL5`/`NKG7`/`GNLY` have no significant rows in any cluster; no cluster pairs leukocyte markers with diluted `EPCAM`/`SCGB1A1`/`SFTPA1`). These exports do not prove every retained cell is a singlet or exclude every other doublet combination.
  - **`PTPRC` in C4 is compatible with tuft-2-like biology.** `PTPRC` has a significant row only in C4 (`+4.43`, `82.7%` vs `8.4%`) alongside the full tuft program, consistent with reported `Ptprc` enrichment in mouse tracheal tuft-2 cells; without the accompanying T/NK program (`CD2/CD3/GZMA/CCL5` without significant rows), it is not by itself grounds for removing C4.
  - **Reassignment check**: prior run `job-44614232` (res `0.4`) resolved `10,029` cells in `11` clusters including C7 (`603` cells, epithelial–T/NK doublets, removed per regrouping proposal). This cleaned run totals `10,029 − 603 = 9,426` cells in `7` clusters with no emergent doublet; the removal worked.
  - **Resolution note**: res `0.3` gives one cluster per canonical epithelial type (merged secretory, merged ciliated, merged tuft), matching the old res-`0.3` node structure. The old res-`0.4` splits (distal vs canonical secretory; apical-tip vs dynein-enriched ciliated; Tuft-1 vs Tuft-2) are available at higher resolution if needed.
- **Master Index & References**: See `../../../NOTES_BARCODE_REASSIGNMENT_AND_CLEANUP.md` (Travaglini et al. 2020 `PMID: 33208946`; Montoro et al. 2018 `PMID: 30069044`; Basil et al. 2022 `PMID: 35355013`; Perniss et al. 2023 `PMID: 37531421`).

---

## Summary Table (`harmony_snn_res.0.3`)

| Cluster | Cells (%) | Proposed Known Subset | QC / Doublet Flag | Selected Identity Markers | Additional Canonical Markers (`markers_significant.csv`) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **1** | 3,171 (33.64%) | **Airway Secretory / Club Cells (`SCGB1A1^{hi} SCGB3A1^{hi} BPIFB1⁺ TFF3⁺`)** | Keep (Clean) | `SCGB1A1` (#1, `88%`), `SCGB3A1` (#2, `59%`), `CRISP3` (#3, `26%`), `BPIFB1` (#4, `56%`), `ATP13A5` (#5, `44%`), `TFF3` (#8, `96%`), `ALPL` (#12, `49%`), `PIGR` (#23, `89%`), `SCGB3A2` (#25, `88%`), `WFDC2` (#34, `91%`) | `SFTPB` (#141, `75%`), `SFTPA1` (#287, `91%`), `S100A2` (#47, `80%`), `KDR` (#28, `39%`), `CLIC6` (#33, `62%`); `PGC` without rows in significant or prefiltered exports |
| **2** | 2,982 (31.64%) | **Airway Multiciliated Cells (`FOXJ1⁺ CFAP⁺ DNAH⁺ RFX⁺`)** | Keep (Clean) | `DNAH2` (#8, `50%`), `CFAP54` (#11, `79%`), `ANKFN1` (#18, `63%`), `RSPH14`, `SNTN` (#72, `72%`), `CAPS` (#103, `74%`), `FOXJ1` (#207, `74%`), `TPPP3` (#208, `81%`) | `CFAP126` (#168, `72%`), `CFAP90` (#144, `56%`), `PIERCE1` (#167, `73%`), `DNAH7` (#115, `46%`), `RFX2` (#275, `65%`), `RFX3` (#378, `45%`) |
| **3** | 1,287 (13.65%) | **Alveolar Type 2 (AT2) Pneumocytes (`SFTPC^{hi} LAMP3^{hi} ABCA3⁺ NAPSA^{hi}`)** | Keep (Clean) | `IL13RA2` (#1, `40%`), `KCNE5` (#2, `35%`), `SFTPC` (#3, `100%`), `LAMP3` (#4, `96%`), `GJB1` (#5, `46%`), `SLC34A2` (#23, `90%`), `ETV5` (#34, `49%`), `NAPSA` (#35, `93%`), `ABCA3` (#31, `63%`) | `SFTPD` (#71, `96%`), `SFTPA1` (#86, `100%`), `SFTPB` (#118, `93%`), `CFTR` (#166, `56%`) |
| **4** | 776 (8.23%) | **Tuft (Brush) Cells, merged Tuft-1/Tuft-2 (`POU2F3^{hi} POU2AF2⁺ PTPRC/CD45⁺`)** | Keep (Clean) | `HMX3` (#1, `52%`), `POU2AF2` (#3, `49%`), `GNG13` (#4, `72%`), `CHAT` (#5, `34%`), `SUCNR1` (#6, `29%`), `SH2D7` (#7, `96%`), `IL17RB` (#8, `76%`), `ALOX5` (#15, `94%`), `POU2F3` (#17, `95%`), `RGS21` (#19, `87%`), `GFI1B` (#21, `28%`), `ALOX5AP` (#25, `88%`), `SPIB`, `PTPRC` (`83%`, tuft-2 hallmark) | `AVIL` (#72, `84%`), `NAV2` (#71, `71%`), `DCLK2` (#82, `50%`); `HTR4` without a significant row; Tuft-1/Tuft-2 split was resolved in the previous run and can be reassessed at higher resolution |
| **5** | 663 (7.03%) | **Alveolar Type 1 (AT1) Pneumocytes (`COL4A3/4^{hi} AGER^{hi} RTKN2^{hi} CAV1⁺`)** | Keep (Clean) | `COL4A3` (#1, `74%`), `COL4A4` (#2, `78%`), `UPK3B` (#4, `43%`), `ANKRD1` (#5, `44%`), `CAV1` (#8, `59%`), `RTKN2` (#13, `84%`) | `AGER` (#29, `98%`), `EMP2` (#68, `81%`), `AQP5` (#124, `78%`), `HOPX` (#116, `97%`), `PDPN` (#252, `90%`) |
| **6** | 292 (3.10%) | **Pulmonary Ionocytes (`FOXI1^{hi} ATP6V1G3^{hi} ATP6V1C2^{hi} CFTR⁺ PDE1C⁺`)** | Keep (Clean) | `PGF` (#2, `78%`), `FOXI1` (#3, `92%`), `ATP6V1G3` (#4, `90%`), `SLC9A2` (#5, `70%`), `ADGRF1` (#8, `92%`), `PDE1C` (#10, `94%`), `ATP7B` (#12, `64%`), `STAP1` (#13, `90%`), `SLC4A11`, `ATP6V1C2` (#41, `91%`) | `CFTR` (#104, `80%`); `ASCL3` without rows in significant or prefiltered exports; `CFAP54` depleted (#2142, `2%`) |
| **7** | 255 (2.71%) | **Airway Basal Progenitor Cells (`TP63⁺ KRT5⁺ KRT17⁺ COL17A1^{hi}`)** | Keep (Clean) | `COL17A1` (#2, `78%`), `KRT5` (#3, `62%`), `TP63` (#4, `62%`), `DLK2` (#5, `34%`), `PKP1` (#11, `39%`), `HAS3` (#13, `57%`), `KRT17` (#14, `52%`) | `ITGA2` (#40, `66%`), `KRT19` (#151, `86%`), `SFN` (#21, `65%`) |

---

## Concise Cluster Rationale

* **C1 (`3,171`) — Secretory**: Canonical club/secretory program (`SCGB1A1` `88%`, `SCGB3A1` `59%`, `BPIFB1` `56%`, `TFF3` `96%`, `WFDC2` `91%`, `PIGR` `89%`) with distal-like leads (`ATP13A5`, `ALPL`). Distal vs canonical split unresolved at this resolution (`CRISP3` was significant in both old secretory clusters).
* **C2 (`2,982`) — Multiciliated**: Ciliary program (`FOXJ1`, `CFAP54/126/90`, `TPPP3`, `SNTN`, `CAPS`, `PIERCE1`, `RSPH14`, `DNAH2`, `ANKFN1`, dynein/TF `RFX2/3`). Apical-tip vs dynein-enriched split unresolved at this resolution.
* **C3 (`1,287`) — AT2**: Definitive AT2 (`SFTPC` `100%`, `LAMP3` `96%`, `NAPSA` `93%`, `SLC34A2` `90%`, `ABCA3`, `SFTPD`/`SFTPA1`/`SFTPB`, `ETV5`, `GJB1`).
* **C4 (`776`) — Tuft (merged)**: Pan-tuft (`POU2F3` `95%`, `POU2AF2`, `SH2D7` `96%`, `ALOX5`/`ALOX5AP`) with both Tuft-2 (`IL17RB`, `SPIB`, `SUCNR1`, `PTPRC/CD45` `83%`) and Tuft-1 (`CHAT`, `NAV2`, `DCLK2`, `RGS21`; `HTR4` without a significant row) features. Split was resolved in the previous run and can be reassessed at higher resolution.
* **C5 (`663`) — AT1**: Basement-membrane AT1 (`COL4A3/4`, `CAV1`, `RTKN2`, `AGER` `98%`, `HOPX` `97%`, `PDPN` `90%`, `AQP5`, `EMP2`).
* **C6 (`292`) — Ionocytes**: V-ATPase/CFTR program (`FOXI1` `92%`, `ATP6V1G3` `90%`, `ATP6V1C2`, `CFTR` `80%`, `PDE1C` `94%`, `SLC9A2`, `ADGRF1`, `ATP7B`, `STAP1`, `PGF`).
* **C7 (`255`) — Basal**: Progenitor program (`COL17A1`, `TP63`, `KRT5`, `KRT17`, `DLK2`, `PKP1`, `HAS3`, `SFN`, `ITGA2`, `KRT19`).

> Note: pre-red-palette outputs (viridis dotplot + prior CSVs/heatmaps) are archived in `old/` in this directory; current plots use the red scale.
