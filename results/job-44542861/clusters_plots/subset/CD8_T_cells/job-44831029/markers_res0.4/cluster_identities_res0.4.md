# CD8_T_cells Subclustering Identity Report (`job-44831029`, Resolution `0.4`)

- **Parent Run / Subset**: `job-44542861` $\rightarrow$ `CD8_T_cells` (`job-44831029`; cleaned baseline with incoming reassignments: `6,783` co-clustered `CD8⁺` cells from `CD4_T_cells`, `4,017 + 2,213` cells from `T_cells_unresolved`)
- **Clustering Parameters**: `harmony_snn_res.0.4` (`33,565` cells across `6` clusters; `2,000` HVGs, `dims 1:15`, `k_param = 100`)
- **QC & Parameter Assessment**:
  - **Keep all `33,565` cells (Clusters 1–6). No doublet cluster.** B (`CD19`, `MS4A1`), plasma (`SDC1`, `JCHAIN`), monocyte (`S100A9`, `FCN1`, `CD14`, `CD163`), endothelial (`PECAM1`) and epithelial (`EPCAM`, `KRT19`) markers have no rows in the significant-marker export, and `CD4` has no significant row in any cluster (`PTPRC` positive in C3, negative in C5). No cluster shows a coherent foreign-lineage signature in these exports; filtered exports cannot prove every cell is a singlet.
  - **Reassignment check**: prior run `job-44614237` resolved `20,552` cells in `6` clusters and proposed incoming merges totalling `+17,966`. Executed CD8 transfers added `13,013` cells (`3,748 + 3,035` from `CD4_T_cells` plus `4,017 + 2,213` from `T_cells_unresolved`), yielding `33,565`. The earlier `+17,966` scenario additionally included `4,953` NK-derived cells, which were instead assigned to `Unconventional_T_and_ILCs`. No emergent doublet program appears in these exports; barcode sets were not independently rechecked here.
  - **Caveat — compartment shift vs old run**: the six identities do not map one-to-one onto the old run. Notably the old `ITGAE⁺ ITGA1⁺ ZNF683⁺` intraepithelial $T_{RM}$ program and the old quiescent `BCL11B⁺` cluster are not recovered as such: `ITGAE`/`ITGA1`/`ZNF683` have no significant rows in new C2, and new C5 is a `SELL⁺ CCR7⁺ LEF1⁺` naive/CM-like cluster absent from the old run.
  - **`k_param` Note**: `k_param = 100` at `res 0.4` resolves the six programs below. If finer cytotoxic or innate-like subdivisions are desired, consider a `k=40–50` sensitivity comparison (untested proposal, not a data-derived optimum).
- **Master Index & References**: See `../../../NOTES_BARCODE_REASSIGNMENT_AND_CLEANUP.md` (Herrera-Uribe et al. 2021 `PMID: 34249103`; Mair et al. 2016 `PMID: 27471504`; Xiao et al. 2019 `PMID: 31552029`).

---

## Summary Table (`harmony_snn_res.0.4`)

| Cluster | Cells (%) | Proposed Known Subset | QC / Doublet Flag | Key Named Top-25 Markers | Additional Canonical Markers (`markers_significant.csv`) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **1** | 8,942 (26.64%) | **`GZMK⁺ CCR5⁺` Effector-Memory `CD8⁺` T Cells (parenchymal $T_{EM}/T_{RM}$-like)** | Keep (Clean) | `GZMK` (#1, `58%`), `TNFSF9`, `GZMA` (#3, `78%`), `TIGIT` (#4, `29%`), `LAG3` (#6, `28%`), `CXCR6` (#7, `38%`), `CCR5` (#8, `58%`), `CSF1`, `CCL4` (#14, `65%`), `PRDM1` (#16, `50%`), `SLAMF7`, `CD8B` (#23, `29%`) | `CD8A` (#41, `69%`); `CXCR3`/`CD69` without significant rows (not assessed) |
| **2** | 8,361 (24.91%) | **`CD96⁺ XCL1⁺ CXCR6⁺` Activated / Resident-like `CD8⁺` T Cells (intraepithelial $T_{RM}$ program not recovered)** | Keep (Clean) | `CD96` (#4, `44%`), `XCL1` (#7, `29%`), `CXCR6` (#10, `38%`), `NR4A1` (#16, `37%`), `JAML` (#17, `34%`), `FOS`/`FOSB`/`EGR1` activation module | `CD69` (#36, `67%`); `KLRK1`, `ITGAE`, `ITGA1`, `ZNF683`, `CRTAM` without significant rows (prior $T_{RM}$ program not recovered as such) |
| **3** | 6,587 (19.62%) | **Vascular `CX3CR1⁺ ADGRG1(GPR56)⁺ GZMB^{hi}` Cytotoxic Effector ($T_{EMRA}/T_{EFF}$) `CD8⁺` T Cells** | Keep (Clean) | `CX3CR1` (#1, `28%`), `IL10` (#3, `30.5%`), `ADGRG1` (`GPR56`, #4, `77%`), `KLF3` (#6, `41%`), `GZMB` (#8, `70%`), `TBX21` (`T-bet`, #16, `54%`), `KLF2` (#20, `83%`), `GNLY` (#22, `61%`) | `GZMH` (#15, `60%`), `NKG7` (#41, `69%`), `S1PR1` (#44, `38%`), `PRF1` (#136, `50%`); `CD8A` without a significant row |
| **4** | 4,081 (12.16%) | **Innate-like `KIT⁺ FCER1G⁺ KLRB1⁺ ZNF683⁺` NK-like T / ILC1-like Cells (`NCR1` not recovered)** | Keep (Clean) | `KIT` (#1, `26%`), `FCER1G` (#3, `39%`), `KLRB1` (`CD161`, #4, `27%`), `ZNF683` (`Hobit`, #7, `26.5%`), `XCL1` (#8, `49%`), `IKZF2` (`Helios`, #13, `62%`), `TXK`, `HOPX`, `TCF7` (#22, `61%`), `KLRK1` (#25, `79%`) | `IL2RB` (#28, `87%`), `GATA3` (#42, `45%`); `NCR1` absent from significant; `CD8A` depleted (#503, `33%`), `CD5`/`CD6` depleted |
| **5** | 4,057 (12.09%) | **Naive / Central-Memory ($T_N / T_{CM}$) `CD8⁺` T Cells (`SELL^{hi} CCR7^{hi} LEF1^{hi}`)** | Keep (Clean) | `LIX1` (#2, `40%`), `CERS6` (#3, `77%`), `SELL` (`CD62L`, #4, `62%`), `CCR7` (#5, `84%`), `LEF1` (#6, `82%`), `BACH2` | `TCF7` (#27, `84%`), `KLF2` (#437, `74%`), `BCL2` (#49, `89%`), `CD27` (#55, `35%`), `CRTAM` (#370, `77%`); `GZMB`/`GZMH` without significant rows, `PRF1` depleted (#945, `30%`); `CD8A` detected (`58%`) but significantly lower in expression |
| **6** | 1,537 (4.58%) | **Type-17 Cytotoxic (`Tc17`) / MAIT-like `IL23R⁺ RORC⁺ KLRB1⁺` T Cells (`ZBTB16` not recovered)** | Keep (Clean) | `IL23R` (#1, `28%`), `RORC` (#2, `27%`), `KLRB1` (`CD161`, #3, `38%`), `TNF` (#5, `62%`), `EOMES` (#7, `31%`), `NFIL3`, `IL7R` (#15, `71%`), `GNLY` (#16, `64%`), `IL18R1` (#18, `51%`), `SH2D1A` | `MAF` (#30, `34%`), `FCGR3A` (#21, `98%`); `ZBTB16` (`PLZF`) without a significant row; `CD4`/`CD8A` without significant rows |

---

## Concise Cluster Rationale

* **C1 (`8,942`)**: Dominant `GZMK⁺` effector-memory pool (`GZMA⁺ CCR5⁺ CXCR6⁺ LAG3⁺ TIGIT⁺ CCL4⁺ PRDM1⁺ SLAMF7⁺`) with `CD8B`/`CD8A` support. `CXCR3` and `CD69` lack significant rows, so finer $T_{RM}$ confirmation is not available from these exports (absence from the export does not establish absent expression).
* **C2 (`8,361`)**: Activated/resident-like program (`CD96⁺ XCL1⁺ CXCR6⁺ JAML⁺ NR4A1⁺ CD69⁺`) on an AP-1 activation module (`FOS/FOSB/EGR1/ATF3/JUN`). The prior intraepithelial $T_{RM}$ integrin/Hobit program is not recovered here.
* **C3 (`6,587`)**: Terminal cytolytic effector program (`CX3CR1⁺ ADGRG1/GPR56^{hi} TBX21⁺ KLF2⁺ KLF3⁺ GZMB^{hi} GZMH⁺ GNLY⁺ NKG7⁺ PRF1⁺ IL10⁺`, `S1PR1⁺`).
* **C4 (`4,081`)**: Innate-like lymphocytes (`KIT⁺ FCER1G⁺ KLRB1/CD161⁺ ZNF683/Hobit⁺ XCL1⁺ IKZF2/Helios⁺ TXK⁺ HOPX⁺ KLRK1⁺ IL2RB⁺`) with depleted `CD8A`/`CD5`/`CD6`. `NCR1 (NKp46)` lacks a significant row, so NKp46 positivity is not established in this run (absence from the export does not establish absent expression).
* **C5 (`4,057`)**: Naive/central-memory program (`SELL^{hi} CCR7^{hi} LEF1^{hi} TCF7⁺ BACH2⁺ BCL2⁺ CD27⁺ KLF2⁺ CRTAM⁺`) with `GZMB`/`GZMH` lacking significant rows and `PRF1` depleted — no counterpart in the old run's six identities.
* **C6 (`1,537`)**: Type-17 cytotoxic program (`IL23R⁺ RORC⁺ KLRB1/CD161⁺ EOMES⁺ NFIL3⁺ MAF⁺ IL7R^{hi} IL18R1⁺ TNF⁺ GNLY⁺ SH2D1A⁺`). `ZBTB16 (PLZF)` lacks a significant row, so the MAIT-like label rests on `IL23R/RORC/KLRB1` rather than PLZF confirmation (absence from the export does not establish absent expression).

> Note: pre-red-palette outputs (viridis dotplot + prior CSVs/heatmaps) are archived in `old/` in this directory; current plots use the red scale.
