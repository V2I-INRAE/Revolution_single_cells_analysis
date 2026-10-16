# NK_cells Subclustering Identity Report (`job-44614238`, Resolution `0.4`)

- **Parent Run / Subset**: `job-44542861` $\rightarrow$ `NK_cells` (`job-44614238`; Parent Cluster 3)
- **Clustering Parameters**: `harmony_snn_res.0.4` (`24,342` cells across `7` clusters; `2,000` HVGs, `dims 1:15`, `k_param = 100`)
- **QC & Reassignment Summary**: **`0` doublet clusters**.
  - **Keep `19,389` true `CD3E⁻` NK / ILC1 cells (Clusters 1, 2, 5, 6, 7)**.
  - **Reassign `4,953` co-clustered `CD3D/E/G⁺ CD247⁺ CD8A⁺ CD8B⁻` NK-like T / NKT cells (Cluster 3 = `2,892` activated `FCGR3A⁺ GNLY^{hi} CSF1⁺ MHC-II⁺`; Cluster 4 = `2,061` resident-memory `TCF7⁺ CRTAM^{hi} ZNF683⁺`) $\rightarrow$ `CD8_T_cells` (or `unconventional_T_cells`)**.
  - **Parameter Recommendation**: Removing C3 & C4 leaves 5 clean `CD3E⁻` NK clusters; optional re-subclustering with **`k_param = 30–40`** (`dims = 1:15`).
- **Master Index & References**: See `../../NOTES_BARCODE_REASSIGNMENT_AND_CLEANUP.md` (Gerner et al. 2009 `PMID: 18601948`; Mair et al. 2016 `PMID: 27471504`; Herrera-Uribe et al. 2021 `PMID: 34249103`; Wiarda et al. 2026 `PMID: 41859083`).

---

## Summary Table (`harmony_snn_res.0.4`)

| Cluster | Cells (%) | Proposed Known Subset | QC / Doublet Flag | Key Named Top-25 Markers | Additional Canonical Markers (`markers_significant.csv`) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **1** | 8,573 (35.22%) | **Tissue-Resident / Cytokine-Responsive `NCR1(NKp46)^{hi} KLRB1^{hi} GZMK⁺ CXCR6⁺ CD69⁺ EOMES^{hi} IL18R1⁺` NK Cells (`CD3E⁻`)** | Keep (Clean) | `GZMK` (`56.2%`), `IL18RAP`, `CXCR6` (`60.1%`), `EOMES` (`53.9%`), `CD69` (`72.2%`), `KLRB1` (`93.5%`), `ISG20`, `TNFSF9`, `CLEC1B`, `IL18R1` (`58.0%`), `NCR1` (`NKp46`, `89.9%`), `SLAMF7`, `SH2D1A`, `BCL2A1` | `CCL4` (#27, `85.8%`), `GZMA` (#64, `98.3%`), `TYROBP` (#96, `53.2%`), `FCER1G` (#209, `62.9%`), `XCL1` (#266, `71.2%`), `FCGR3A` (#267, `91.6%`), `NKG7` (#305, `83.4%`); **`CD3E⁻` (`5.5%`)** |
| **2** | 7,080 (29.09%) | **Mature Cytolytic Intravascular `CX3CR1⁺ S1PR5⁺ KLF2^{hi} GZMB^{hi} GNLY^{hi} ITGAM(CD11b)⁺ KLRG1⁺` NK Cells (`CD3E⁻`)** | Keep (Clean) | `CX3CR1` (`35.7%`), `PTH1R`, `GZMH` (`48.7%`), `PTAFR`, `S1PR5` (`43.2%`), `GZMB` (`71.2%`), `GNLY` (`94.3%`), `RASA3`, `VCL`, `KLRG1` (`35.6%`), `KLF2` (`83.1%`), `ITGAM` (`CD11b`, `63.4%`), `CD226` (`45.0%`) | `ZNF683` (#45, `34.7%`), `TBX21` (#67, `66.4%`), `CD8A` (#142, `46.6%`, porcine NK `CD8αα`), `FCGR3A` (#167, `90.5%`), `PRF1` (#233, `74.0%`), `IFNG` (#240, `25.6%`), `TYROBP` (#296, `42.8%`), `ID2` (`89.5%`); **`CD3E⁻` (`7.0%`)** |
| **3** | 2,892 (11.88%) | **Activated `CD3⁺ CD8αα⁺ FCGR3A(CD16)⁺ GNLY^{hi} CSF1⁺ MHC-II⁺` Cytotoxic NK-like T / NKT Cells** | **Reassign (Not Doublet) $\rightarrow$ `CD8_T_cells` / `unconventional_T_cells`** | `CD3E` (`70.1%` vs `14.3%`, `log2FC = 2.50`), `CSF1` (`35.8%`), `FRMD4B`, `SLA-DQB1` (`71.0%`), `HLA-DRA` (`44.9%`), `CD8A` (`63.3%`), `SPOCK2`, `CD3G` (`63.9%`), `ITGB7` (`46.0%`), `SLA-DRB1` (`67.5%`), `CD52`, `GRK3` | `TNF` (#26, `26.5%`), `CD74` (#39, `69.6%`), `CD96` (#40, `46.3%`), `CD3D` (#57, `75.1%`), `CD247` (#105, `65.9%`), `PRDM1` (#122, `51.4%`), `CRTAM` (#142, `58.1%`), `GNLY` (`77.8%`), `FCGR3A` (`77.4%`); **`TYROBP⁻` (`8.3%`), `CD8B⁻`** |
| **4** | 2,061 (8.47%) | **Resident-Memory `CD3⁺ CD8αα⁺ TCF7⁺ CRTAM^{hi} ZNF683⁺ GZMM⁺` Cytotoxic T / NKT Cells** | **Reassign (Not Doublet) $\rightarrow$ `CD8_T_cells` / `unconventional_T_cells`** | `CD3E` (`75.3%` vs `15.9%`), `ABTB2`, `TCF7` (`59.6%`), `P2RY14` (`45.8%`), `GZMM` (`29.2%`), `CD3G` (`72.1%`), `ZNF683` (`Hobit`, `35.6%`), `CRTAM` (`81.3%`), `SPOCK2`, `SEMA7A`, `BCL2` (`80.3%`) | `CD8A` (#29, `60.8%`), `XCL1` (#40, `79.7%`), `CD247` (#69, `72.3%`), `CD96` (#71, `50.8%`), `IKZF2` (#130, `77.3%`), `SLAMF6` (#140, `28.9%`), `CD3D` (#145, `80.5%`); **`TYROBP⁻` (`7.7%`), `CD8B⁻`** |
| **5** | 1,278 (5.25%) | **Basal / Low-Complexity Cytotoxic `KLF2⁺ GZMB⁺ PRDM1⁺ ANXA1⁺ BTC⁺` NK State** | Keep (Clean) | `BTC` (`28.1%`), `ANXA1` (`34.0%`), `TXNRD1` (`26.3%`), `AGPS`, `RGCC`, `ERC1`, `SLC6A6`, `BACH1`, `MED13L`, `ELMO1`, `ELL2`, `KLF9`, `VIM`, `DOCK8` | `KLF2` (#37, `54.6%`), `GZMB` (#42, `47.6%`), `PRDM1` (#44, `53.3%`), `NCR1` (`46.6%`) |
| **6** | 1,246 (5.12%) | **Quiescent `BACH2⁺ ZEB1⁺ ZBTB16(PLZF)⁺ IKZF2(Helios)^{hi} NCR1⁺` Innate NK / ILC Subset** | Keep (Clean) | `GPR39`, `THSD7B`, `SMYD3`, `KCNQ5`, `WWOX`, `AFF3` (`68.9%`), `PIK3CD`, `IMMP2L`, `BTBD9`, `BCAS3`, `BACH2` (`38.9%`), `CHST11` (`90.3%`), `MAML2`, `KLF12` (`54.6%`), `ZEB1` (`57.8%`), `PREX1` (`58.7%`) | `ZBTB16` (`PLZF`, #78, `47.8%`), `IKZF2` (`Helios`, #240, `82.5%`), `CD96` (#575, `41.4%`), `NCR1` (#611, `60.6%`), `IL18R1` (#622, `37.6%`); `CD3G^{lo}` (`18.0%`) |
| **7** | 1,212 (4.98%) | **Immature / Progenitor `TCF7^{hi} NCR1^{hi} GATA3⁺ CD96^{hi} IL18R1⁺ CD69⁺ XCL1⁺ GNLY⁻ GZMB⁻` NK / ILC1-like Cells** | Keep (Clean) | `TCF7` (`67.8%`, `log2FC = 2.26`), `SPOCK2` (`54.8%`), `SLC38A1`, `SERPINB1`, `RBPJ` (`50.3%`), `EEIG1`, `VIM`, `P2RY14` (`53.6%`), `ZNF683` (`32.3%`), `CD96` (`61.4%`), `GPR183` (`68.9%`) | `IL18R1` (#33, `54.5%`), `GATA3` (#45, `56.8%`), `NCR1` (#87, `83.0%`), `CD69` (#126, `64.5%`), `XCL1` (#132, `71.9%`), `ZBTB16` (#167, `33.7%`), `IL2RB` (`97.0%`), `TYROBP` (`47.9%`); `GZMB⁻` (`10.1%`) |

---

## Concise Cluster Rationale

* **True `CD3E⁻` NK / ILC1 Subsets (`Clusters 1, 2, 5, 6, 7` = `19,389` cells)**:
  * **C1 (`8,573`) vs. C2 (`7,080`)**: Textbook porcine/mammalian NK dichotomy (Mair et al. 2014/2016 `PMID: 27471504`; Wiarda et al. 2026 `PMID: 41859083`): **C1** is tissue-resident/cytokine-responsive `NKp46(NCR1)^{hi} KLRB1^{hi} GZMK⁺ GZMA^{hi} CXCR6⁺ CD69⁺ EOMES^{hi} IL18R1⁺ XCL1⁺` NK, whereas **C2** is mature intravascular cytolytic `CX3CR1⁺ S1PR5⁺ KLF2^{hi} TBX21^{hi} ITGAM(CD11b)⁺ KLRG1⁺ GZMB^{hi} GNLY^{hi} PRF1^{hi} CD8αα⁺` NK.
  * **C5 (`1,278`), C6 (`1,246`), C7 (`1,212`)**: C5 is a low-UMI basal `KLF2⁺ GZMB⁺` NK state; C6 is a quiescent `BACH2⁺ ZEB1⁺ ZBTB16⁺ IKZF2^{hi} NCR1⁺` innate subset; C7 is an immature/progenitor `TCF7^{hi} NCR1^{hi} GATA3⁺ CD96^{hi} XCL1⁺ GZMB⁻ GNLY⁻` NK/ILC1 subset.
* **Co-Clustered `CD3⁺ CD8αα⁺` NK-like T / NKT Subsets (`Clusters 3 & 4` = `4,953` cells $\rightarrow$ Reassign)**:
  * **C3 (`2,892`) & C4 (`2,061`)**: Express full TCR/CD3 subunits (`CD3D` `75–81%`, `CD3E` `70–75%`, `CD3G` `64–72%`, `CD247` `66–72%`) and `CD8A` (`61–63%`, `CD8B⁻`), while lacking `TYROBP (DAP12)` (`8%`). C3 is acutely activated (`FCGR3A^{hi} GNLY^{hi} CSF1⁺ TNF⁺ MHC-II^{hi}`); C4 is resident-memory (`TCF7⁺ CRTAM^{hi} ZNF683⁺ GZMM⁺ XCL1^{hi}`).
