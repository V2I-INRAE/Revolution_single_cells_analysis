# CD8_T_cells Subclustering Identity Report (`job-44831029`, Resolution `0.2`)

- **Parent Run / Subset**: `job-44542861` $\rightarrow$ `CD8_T_cells` (`job-44831029`; cleaned baseline with incoming reassignments; `33,565` cells)
- **Clustering Parameters**: `harmony_snn_res.0.2` (`33,565` cells across `5` clusters; `2,000` HVGs, `dims 1:15`, `k_param = 100`)
- **QC & Parameter Assessment**:
  - **Keep all `33,565` cells (Clusters 1–5). No doublet cluster.** B (`CD19`, `MS4A1`), plasma (`SDC1`, `JCHAIN`), monocyte (`S100A9`, `FCN1`, `CD14`), endothelial (`PECAM1`) and epithelial (`EPCAM`) markers have no significant rows in any cluster, and `CD4` has no significant row. No cluster shows a coherent foreign-lineage signature; filtered exports cannot prove every cell is a singlet.
  - **Reassignment check**: res `0.4` resolved the same `33,565` cells in `6` clusters. This res-`0.2` clustering is marker-compatible with consolidation into `5`: new C1 (`17,510`) matches old `GZMK⁺` effector-memory C1 (`8,942`) plus old `CD96⁺` activated/resident-like C2 (`8,361`) = `17,303` vs `17,510`; new C2 (`7,153`) matches old cytotoxic-effector C3 (`6,587`) plus membership shifts; new C3 (`4,081`) matches old naive/CM-like C5 (`4,057`); new C4 (`3,340`) is marker-compatible with old innate-like C4 (`4,081`) at reduced membership; new C5 (`1,481`) matches old Tc17-like C6 (`1,537`). Cross-run membership is barcode-unverified; counts are compared, not traced.
- **Master Index & References**: See `../../../NOTES_BARCODE_REASSIGNMENT_AND_CLEANUP.md` (Herrera-Uribe et al. 2021 `PMID: 34249103`; Mair et al. 2016 `PMID: 27471504`; Xiao et al. 2019 `PMID: 31552029`).

---

## Summary Table (`harmony_snn_res.0.2`)

| Cluster | Cells (%) | Proposed Known Subset | QC / Doublet Flag | Selected Identity Markers (`markers_significant.csv`) |
| :--- | :--- | :--- | :--- | :--- |
| **1** | 17,510 (52.17%) | **`GZMK⁺ CXCR6⁺ CCR5⁺` Effector-Memory / Activated `CD8⁺` T Cells** | Keep (Clean) | `CXCR6` (#1, `37.9%`), `GZMK` (#2, `39.7%`), `TNFSF9` (#3, `33.1%`), `CD96` (#5, `35.2%`), `CCR5` (#6, `48.5%`), `GZMA` (#22, `57.0%`), `SLAMF7` (#21, `31.9%`), `GZMB`-depleted (#505), `CD8A` (#83, `63.2%` vs `53.1%`); `CD8B`: no significant row |
| **2** | 7,153 (21.31%) | **Vascular `CX3CR1⁺ ADGRG1/GPR56⁺ GZMB^{hi}` Cytotoxic Effector `CD8⁺` T Cells** | Keep (Clean) | `CX3CR1` (#1, `27.2%`), `ADGRG1` (`GPR56`, #3, `76.9%`), `KLF3` (#6, `40.7%`), `GZMB` (#8, `68.6%`), `TBX21` (`T-bet`, #16, `53.6%`), `KLF2` (#18, `82.9%`), `GZMH` (#19, `57.4%`), `GNLY` (#15, `62.4%`), `PRF1` (#110, `50.8%`), `NKG7` (#36, `68.6%`) |
| **3** | 4,081 (12.16%) | **Naive / Central-Memory `CD8⁺` T Cells (`SELL^{hi} CCR7^{hi} LEF1^{hi} TCF7⁺ BACH2⁺`)** | Keep (Clean) | `SELL` (`CD62L`, #4, `61.3%`), `CCR7` (#5, `84.1%`), `LEF1` (#6, `82.1%`), `BACH2` (#17, `56.9%`), `TCF7` (#27, `83.4%`), `BCL2` (#50, `89.3%`), `CD27` (#54, `34.6%`), `IL7R` (#42, `81.8%`); `GZMB`/`GZMH`: no significant rows |
| **4** | 3,340 (9.95%) | **Innate-like `NCR1⁺ KIT⁺ FCER1G⁺ KLRB1⁺ ZNF683⁺` NK-like T / ILC1-like Cells** | Keep (Clean) | `NCR1` (#1, `28.0%`), `KIT` (#2, `31.1%`), `FCER1G` (#7, `40.7%`), `KLRB1` (`CD161`, #5, `31.1%`), `ZNF683` (`Hobit`, #9, `29.3%`), `XCL1` (#11, `53.6%`), `HOPX` (#15, `53.3%`), `IKZF2` (`Helios`, #16, `64.7%`), `TCF7` (#22, `69.4%`), `IL2RB` (#37, `87.1%`); `CD8A`-depleted (#587), `CD3D`-depleted (#290), `CD3E`-depleted (#514) |
| **5** | 1,481 (4.41%) | **Type-17 Cytotoxic (`Tc17`) / MAIT-like `IL23R⁺ RORC⁺ KLRB1⁺` T Cells** | Keep (Clean) | `IL23R` (#1, `28.4%`), `RORC` (#2, `28.0%`), `KLRB1` (`CD161`, #3, `39.2%`), `TNF` (#5, `62.7%`), `EOMES` (#7, `32.0%`), `IL7R` (#15, `71.3%`), `GNLY` (#16, `64.6%`), `IL18R1` (#18, `51.9%`) |

---

## Concise Cluster Rationale

* **Merged Effector-Memory/Activated (C1 = `17,510`)**: At res `0.2` the res-`0.4` `GZMK⁺` effector-memory and `CD96⁺` activated/resident-like programs co-cluster (`CXCR6` #1, `GZMK` #2, `TNFSF9` #3, `CD96` #5, `CCR5` #6, `GZMA` #22, `SLAMF7` #21; `CD8A` #83 with `63.2%` detection). `GZMB` is depleted (#505), consistent with a non-terminal-effector state.
* **Cytotoxic Effectors (C2 = `7,153`)**: `CX3CR1⁺ ADGRG1/GPR56⁺` vascular effectors with `GZMB^{hi}`/`GZMH`/`GNLY`/`PRF1`/`NKG7` and `TBX21`/`KLF2`/`KLF3`.
* **Naive/CM (C3 = `4,081`)**: `SELL`/`CCR7`/`LEF1`/`TCF7`/`BACH2`/`BCL2`/`CD27`/`IL7R` program; `GZMB`/`GZMH` have no significant rows, and `GZMK`/`GZMA` appear only as depletion rows.
* **Innate-like (C4 = `3,340`)**: `NCR1`-led NK-like program (`KIT`, `FCER1G`, `KLRB1/CD161`, `ZNF683/Hobit`, `XCL1`, `HOPX`, `IKZF2/Helios`, `IL2RB`) recovered at res `0.2` where res-`0.4` C4 lacked an `NCR1` row; `CD8A`/`CD3D`/`CD3E` depleted rows consistent with a non-conventional-T identity (depletion rows describe relative expression, not absent expression).
* **Tc17-like (C5 = `1,481`)**: `IL23R`/`RORC`/`KLRB1` with `TNF`, `EOMES`, `IL7R`, `GNLY`, `IL18R1`; `ZBTB16` has no significant row.

> Note: pre-red-palette outputs (viridis dotplot + prior CSVs/heatmaps) are archived in `old/` in this directory; current plots use the red scale.
