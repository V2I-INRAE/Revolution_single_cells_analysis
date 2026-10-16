# CD8_T_cells Subclustering Identity Report (`job-44614237`, Resolution `0.4`)

- **Parent Run / Subset**: `job-44542861` $\rightarrow$ `CD8_T_cells` (`job-44614237`; Parent Cluster 5)
- **Clustering Parameters**: `harmony_snn_res.0.4` (`20,552` cells across `6` clusters; `2,000` HVGs, `dims 1:20`, `k_param = 100`)
- **QC & Reassignment Summary**: **`0` doublet clusters, `0` outgoing reassignments (Keep all `20,552` cells)**.
  - **Incoming Reassignments (`+17,966` cells $\rightarrow$ `38,518` total `CD8_T_cells`)**: Merge `CD4_T_cells` C3 (`3,748` Naive/CM `CD8⁺`) & C4 (`3,035` `GZMK⁺` EM `CD8⁺`), `T_cells_unresolved` C1 (`4,017` `GZMK⁺ CXCR6⁺` $T_{EM}/T_{RM}$) & C3 (`2,213` `ADGRG1^{hi} GZMB^{hi}` $T_{EMRA}$), and optionally `NK_cells` C3 (`2,892`) & C4 (`2,061`) (`CD3⁺ CD8αα⁺` NK-like T/NKT).
  - **Parameter Recommendation**: Re-subcluster the unified `CD8_T_cells` compartment with **`k_param = 40–50`** (`dims = 1:20`).
- **Master Index & References**: See `../../NOTES_BARCODE_REASSIGNMENT_AND_CLEANUP.md` (Herrera-Uribe et al. 2021 `PMID: 34249103`; Mair et al. 2016 `PMID: 27471504`; Xiao et al. 2019 `PMID: 31552029`).

---

## Summary Table (`harmony_snn_res.0.4`)

| Cluster | Cells (%) | Proposed Known Subset | QC / Doublet Flag | Key Named Top-25 Markers | Additional Canonical Markers (`markers_significant.csv`) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **1** | 8,459 (41.16%) | **Parenchymal `GZMK⁺ CXCR6⁺ CCR5⁺` Effector-Memory / Tissue-Resident ($T_{EM}/T_{RM}$) `CD8αβ⁺` T Cells** | Keep (Clean) | `GZMK` (`52%`), `TNFSF9` (`41%`), `CXCR6` (`43%`), `LAG3` (`28%`), `CXCR3` (`26%`), `CCR5` (`56%`), `TIGIT` (`27%`), `GZMA` (`73%`), `CD8B` (`28%`), `JAML`, `SLAMF7`, `GPR183` (`61%`), `LGALS1`, `CSF1` | `PRDM1` (`Blimp-1`, #27, `48%`), `SLA-DQB1` (#28, `76%`), `SLA-DRB1` (#33, `69%`), `CD69` (#36, `60%`), `CCL4` (#37, `61%`), `CD8A` (#44, `68%`), `FOXO1` (#54, `36%`) |
| **2** | 5,284 (25.71%) | **Vascular `CX3CR1⁺ ADGRG1(GPR56)⁺ GZMB^{hi} KLF2⁺` Cytotoxic Effector ($T_{EMRA}/T_{EFF}$) `CD8⁺` T Cells** | Keep (Clean) | `CX3CR1` (`29%` vs `1%`), `ADGRG1` (`GPR56`, `78%`), `GZMB` (`68%`), `GZMH` (`60%`), `KLF2` (`84%`), `KLF3` (`42%`), `S1PR1` (`38%`), `IL10` (`29%`), `RASA3` (`51%`), `ITGB2` (`58%`), `ADD3`, `FLNA`, `VCL` | `TBX21` (`T-bet`, #27, `55%`), `GNLY` (#30, `61%`), `S100A4` (#47, `60%`), `NKG7` (#58, `70%`), `CD5` (#69, `63%`), `CD6` (#75, `50%`), `PRF1` (#199, `51%`); `CXCR6⁻` (`3%`) |
| **3** | 2,731 (13.29%) | **Innate-like `NCR1(NKp46)⁺ FCER1G⁺ KIT⁺ ZNF683(Hobit)⁺` NK-like T / ILC1-like Cells** | Keep (Clean) | `KIT` (`31%`), `P2RY14` (`30%`), `NCR1` (`NKp46`, `30%`), `KLRB1` (`34%`), `FCER1G` (`43%`), `TCF7` (`70%`), `ZNF683` (`Hobit`, `30%`), `XCL1` (`58%`), `IKZF2` (`Helios`, `65%`), `TXK` (`62%`), `SPRY2`, `FGL2`, `SATB1`, `HOPX` | `GATA3` (#57, `50%`), `CD96` (#60, `50%`), `KLRK1` (#66, `83%`), `IL2RB` (#67, `90%`); `CD5^{lo}` (`14%`), `CD6^{lo}` (`10%`), `CD8A^{lo}` (`29%`) |
| **4** | 2,004 (9.75%) | **Intraepithelial `ITGAE(CD103)⁺ ITGA1(CD49a)⁺ ZNF683⁺ CD96^{hi}` Tissue-Resident Memory ($T_{RM}$) `CD8⁺` T Cells** | Keep (Clean) | `ITGAE` (`CD103`, `34%`), `ITGA1` (`CD49a`, `34%`), `CD96` (`76%`), `XCL1` (`55%`), `ZNF683` (`Hobit`, `29%`), `IKZF2` (`64%`), `TCF7` (`54%`), `NR4A1` (`41%`), `FCER1G` (`26%`), `GCSAM`, `RHOB`, `ATF3`, `FOSB` | `DUSP4` (#26, `38%`), `KLRK1` (#38, `81%`), `JAML` (#40, `33%`), `CRTAM` (#48, `65%`), `CXCR6` (#50, `45%`), `CD69` (#99, `65%`), `CD3G` (#139, `65%`) |
| **5** | 1,390 (6.76%) | **MAIT-like / Type-17 Cytotoxic (`Tc17`) `IL23R⁺ RORC⁺ KLRB1⁺ IL7R^{hi} IL18R1⁺` T Cells** | Keep (Clean) | `IL23R` (`28%` vs `1%`), `RORC` (`28%` vs `2%`), `KLRB1` (`CD161`, `41%`), `IL7R` (`72%`), `IL18R1` (`53%`), `TNF` (`60%`), `GNLY` (`60%`), `LTB` (`37%`), `EOMES` (`32%`), `NFIL3` (`60%`), `SH2D1A` (`55%`), `EGR1`, `MYC`, `AQP3`, `GPR183` | `MAF` (#29, `35%`), `FCGR3A` (#30, `97%`), `CCL4` (#34, `74%`), `CD69` (#37, `63%`), `LAG3` (#39, `39%`), `S1PR1` (#53, `33%`), `ZBTB16` (`PLZF`, #63, `26%`) |
| **6** | 684 (3.33%) | **Quiescent / Resting `BCL11B⁺ RUNX1⁺ SKAP1⁺` `CD8⁺` T Cells** | Keep (Clean) | `THSD7B`, `BCL11B`, `IMMP2L`, `SMYD3`, `WWOX`, `MAML2`, `TSPAN5`, `ATXN1`, `PREX1`, `BCAS3`, `XYLT1` (`70%`), `CDK6`, `AFF3` (`56%`), `TRPS1`, `MSI2` (`57%`), `BTBD9`, `SSH2` (`82%`), `SKAP1` (`73%`), `CHST11` (`73%`) | `RUNX1` (#32, `70%`), `KLF12` (#38, `31%`), `TOX` (#57, `43%`), `CD226` (#304, `33%`), `RUNX3` (#338, `49%`) |

---

## Concise Cluster Rationale

* **C1 (`8,459`)**: Classical $\alpha\beta$ (`CD8A⁺ CD8B⁺`) parenchymal $T_{EM}/T_{RM}$ expressing `GZMK^{hi} GZMA^{hi} CXCR6⁺ CCR5⁺ CXCR3⁺ CD69⁺ PRDM1⁺` and porcine activation MHC class II (`SLA-DQB1`, `SLA-DRB1`).
* **C2 (`5,284`)**: Intravascular/marginated terminal cytolytic effector $T_{EMRA}/T_{EFF}$ (`CX3CR1⁺ ADGRG1/GPR56^{hi} KLF2^{hi} KLF3⁺ S1PR1⁺ TBX21⁺ GZMB^{hi} GZMH^{hi} GNLY^{hi} PRF1⁺ IL10⁺`, `CXCR6⁻ CD69^{lo}`).
* **C3 (`2,731`)**: Porcine pulmonary `CD3⁺ NKp46(NCR1)⁺ FCER1G⁺ KIT⁺ ZNF683(Hobit)⁺ XCL1⁺` innate-like T cells (Mair et al. 2016 `PMID: 27471504`).
* **C4 (`2,004`)**: Intraepithelial mucosal $T_{RM}$ defined by epithelial retention integrins `ITGAE (CD103)` and `ITGA1 (CD49a)` plus `ZNF683 (Hobit)`, `CD96^{hi}`, `CRTAM`, `JAML`, and `XCL1`.
* **C5 (`1,390`)**: Porcine MAIT-like / `Tc17` cells co-expressing `IL23R`, `RORC`, `KLRB1 (CD161)`, `ZBTB16 (PLZF)`, `MAF`, `IL7R^{hi}`, `IL18R1`, `TNF`, and `GNLY` (Xiao et al. 2019 `PMID: 31552029`).
* **C6 (`684`)**: Quiescent `BCL11B⁺ RUNX1⁺ KLF12⁺` `CD8⁺` T cells.
