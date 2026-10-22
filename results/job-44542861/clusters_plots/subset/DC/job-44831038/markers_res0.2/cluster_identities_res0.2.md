# DC Subclustering Identity Report (`job-44831038`, Resolution `0.2`)

- **Parent Run / Subset**: `job-44542861` $\rightarrow$ `DC` (`job-44831038`; cleaned sensitivity run, `k=50`; Parent Clusters 16 [`4,542` pDC] + 19 [`2,894` cDC])
- **Clustering Parameters**: `harmony_snn_res.0.2` (`7,436` cells across `4` clusters; `1,500` HVGs, `dims 1:10`, `k_param = 50`)
- **QC & Parameter Assessment**: **No coherent doublet-cluster signature identified in these exports; no reassignment indicated.**
  - **Resolution comparison**: at res `0.3` the same `7,436` cells resolved into `5` clusters with the pDC compartment split in two (old C1 `2,841` + C2 `1,756` = `4,597`). At res `0.2` the pDC program reunites into a single cluster (new C1 `4,487`, count-compatible with `4,597`), while the cDC2/DC3-like, cDC1, and mregDC-like programs persist as new C2 (`1,459`), C3 (`1,160`), and C4 (`330`). Cross-run membership is barcode-unverified; counts are compared, not traced.
  - T-cell (`CD3D/E`), B-cell (`MS4A1/CD19`), plasma (`JCHAIN`), monocyte (`S100A9`), and epithelial (`EPCAM`) markers have no significant rows in any cluster. `PECAM1` has positive rows in C2 (#565, `+0.68`) and C3 (#438, `+0.90`) with a depletion row in C1; no endothelial program accompanies it (no `VWF`/`CDH5` rows), so `PECAM1` alone does not warrant endothelial reassignment. Filtered marker summaries cannot establish that every cell is a singlet.
- **Master Index & References**: See `../../../NOTES_BARCODE_REASSIGNMENT_AND_CLEANUP.md` (Auray et al. 2016 `PMID: 27837108`; Maisonnasse et al. 2016 `PMID: 26530136`; Mair et al. 2022 `PMID: 35371050`; Baillou et al. 2025 `PMID: 41235223`).

---

## Summary Table (`harmony_snn_res.0.2`)

| Cluster | Cells (%) | Proposed Known Subset | QC / Doublet Flag | Selected Identity Markers (`markers_significant.csv`) |
| :--- | :--- | :--- | :--- | :--- |
| **1** | 4,487 (60.34%) | **Plasmacytoid DCs (`TCF4⁺ IRF7⁺ BCL11A⁺ TLR7/9/3⁺ CD4⁺ BLNK⁺ NCR1⁺ CMKLR1⁺ VPREB1⁺ MGP⁺`; reunited at res `0.2`)** | Keep (Clean) | `VPREB1` (#1, `76.6%`), `NCR1` (#7, `60.3%`), `CD4` (#11, `87.8%`), `MGP` (#14, `90.0%`), `BLNK` (#23, `97.4%`), `CMKLR1` (#36, `69.0%`), `IRF7` (#135, `48.4%`), `TCF4` (#145, `94.7%`), `BCL11A` (#230, `92.5%`), `TLR7` (#327, `66.5%`), `TLR9` (#444, `30.8%`) |
| **2** | 1,459 (19.62%) | **cDC2-like DCs with Inflammatory / DC3-like Features** | Keep (Clean) | `FCER1A` (#1, `31.7%`), `IL1B` (#3, `37.7%`), `CD1.1 (review)` (#4, `63.9%`), `CD207` (#5, `25.4%`), `MS4A2` (#9, `88.2%`), `TLR4` (#7, `36.8%`), `MRC1` (#14, `46.4%`), `CD163` (#16, `29.7%`), `KLF4` (#71, `76.3%`), `SIRPA` (#184, `57.2%`) |
| **3** | 1,160 (15.60%) | **Classical Type 1 DCs (`cDC1`)** | Keep (Clean) | `XIRP2` (#1, `53.0%`), `CLNK` (#2, `49.3%`), `ANPEP` (`CD13`, #3, `49.2%`), `CADM1` (#7, `77.3%`), `XCR1` (#8, `79.5%`), `SNX22` (#11), `BATF3`, `LY75`, `ZNF366`, `FLT3` (#312, `76.6%`) |
| **4** | 330 (4.44%) | **Mature / Migratory `mregDC`-like State** | Keep (Clean) | `TMEM150C` (#1, `72.4%`), `IDO1` (#3, `75.8%`), `FSCN1` (#6, `81.2%`), `CCR7` (#42, `95.8%`), `CD274` (#51, `56.7%`), `CD40` (#67, `84.2%`), `LAMP3` (#173, `93.0%`) |

---

## Concise Cluster Rationale

* **C1 (`4,487`) — reunited pDCs**: The full porcine pDC program (`VPREB1` #1, `NCR1` #7, `CD4` #11, `MGP` #14, `BLNK` #23 `97.4%`, `CMKLR1` #36, `IRF7` #135, `TCF4` #145 `94.7%`, `BCL11A` #230, `TLR7` #327, `TLR9` #444) in one cluster, marker-compatible with reunion of the two res-`0.3` pDC states (old C1 `2,841` + C2 `1,756` = `4,597` vs `4,487`). Whether the res-`0.3` split reflected biological states or technical variation is not established by these exports.
* **C2 (`1,459`) — cDC2 with DC3-like features**: `FCER1A⁺ CD1.1⁺ CD207⁺ MS4A2⁺ TLR4⁺ MRC1⁺ CD163⁺ IL1B⁺` with `KLF4` (#71) and `SIRPA` (#184); count-compatible with res-`0.3` C3 (`1,460`). Finer cDC2/DC3 identity unresolved.
* **C3 (`1,160`) — `cDC1`**: `XIRP2⁺ CLNK⁺ ANPEP/CD13⁺ CADM1^{hi} XCR1^{hi} SNX22⁺ BATF3⁺ LY75⁺ ZNF366⁺ FLT3⁺`; count-compatible with res-`0.3` C4 (`1,049`) plus membership shifts.
* **C4 (`330`) — `mregDC`-like**: identical cell count to res-`0.3` C5 (`TMEM150C⁺ IDO1⁺ FSCN1⁺ CCR7^{hi} CD274⁺ CD40⁺ LAMP3⁺`); whether the same cells are retained is barcode-unverified.

> Note: pre-red-palette outputs (viridis dotplot + prior CSVs/heatmaps) are archived in `old/` in this directory; current plots use the red scale.
