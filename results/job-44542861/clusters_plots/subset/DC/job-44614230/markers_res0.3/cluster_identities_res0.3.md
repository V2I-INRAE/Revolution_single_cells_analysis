# DC Subclustering Identity Report (`job-44614230`, Resolution `0.3`)

- **Parent Run / Subset**: `job-44542861` $\rightarrow$ `DC` (`job-44614230`; Parent Clusters 16 [`4,542` pDC] + 19 [`2,889` cDC])
- **Clustering Parameters**: `harmony_snn_res.0.3` (`7,436` cells across `4` clusters; `1,500` HVGs, `dims 1:20`, `k_param = 100`)
- **QC & Parameter Assessment**: **`0` doublet clusters, `0` reassignments**.
  - **`k_param` Note**: Because `k_param = 100` spans `1.35%` of this `7,436`-cell subset, Cluster 2 (`1,460` cells) merges canonical `FCER1A⁺ CD1.1⁺ CD207⁺` `cDC2` with inflammatory `MRC1⁺ CD163⁺ FCGR3A⁺ TLR4⁺ IL1B⁺` `DC3 / MoDC`. If splitting `cDC2` from `DC3/MoDC` is desired, re-run `subcluster.R` with **`k_param = 25–30`** and **`dims = 1:15`**.
- **Master Index & References**: See `../../NOTES_BARCODE_REASSIGNMENT_AND_CLEANUP.md` (Auray et al. 2016 `PMID: 27798170`; Maisonnasse et al. 2016 `PMID: 26530136`; Mair et al. 2022 `PMID: 35371050`; Baillou et al. 2025 `PMID: 41235223`).

---

## Summary Table (`harmony_snn_res.0.3`)

| Cluster | Cells (%) | Proposed Known Subset | QC / Doublet Flag | Key Named Top-25 Markers | Additional Canonical Markers (`markers_significant.csv`) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **1** | 4,499 (60.50%) | **Plasmacytoid Dendritic Cells (`pDC`, `CD4^{hi} TCF4^{hi} IRF8^{hi} BCL11A^{hi} VPREB1⁺ BLNK^{hi} NCR1⁺ TLR7/9⁺`)** | Keep (Clean) | `VPREB1` (`76%`), `IGLL1` (`73%`), `SMIM5` (`76%`), `ATP13A2` (`95%`), `CD4` (`88%`), `CD8B` (`52%`), `BLNK` (`97%`), `NCR1` (`60%`), `TDRD1` (`69%`), `MGP` (`90%`), `COBLL1` (`80%`), `CD93` (`46%`), `BHLHA15` (`34%`) | `CMKLR1` (#38, `69%`), `IRF7` (#134, `48%`), `TCF4` (`E2-2`, #150, `95%`), `TLR3` (#176, `58%`), `CLEC12A` (#221, `68%`), `BCL11A` (#238, `92%`), `TLR7` (#334, `66%`), `IRF4` (#363, `70%`), `PLD4` (#409, `83%`), `TLR9` (#443, `31%`), `IRF8` (#476, `99%`), `RUNX2` (#493, `81%`) |
| **2** | 1,460 (19.63%) | **Classical Type 2 DCs (`cDC2`, `FCER1A⁺ CD1.1⁺ CD207⁺ SIRPA⁺ KLF4⁺`) + Inflammatory `DC3 / MoDC` (`MRC1⁺ CD163⁺ FCGR3A⁺ IL1B⁺`)** | Keep (Clean) | `FCER1A` (`32%` vs `1%`), `CD1.1` (`64%` vs `2%`), `CD207` (`Langerin`, `25%` vs `1%`), `S100A4` (`71%`), `MS4A2` (`88%`), `FCGR3A` (`79%`), `MRC1` (`CD206`, `46%`), `CD163` (`30%`), `IL1B` (`38%`), `IL1A` (`49%`), `TLR4` (`37%`), `CXCL2` (`55%`), `AREG` (`65%`), `AXL` (`36%`), `CD300C` (`48%`) | `CEBPB` (#30, `78%`), `FCGR2B` (#32, `63%`), `KIT` (#47, `38%`), `KLF4` (#70, `76%`), `CD80` (#90, `56%`), `SLA-DQB1` (`100%`), `SLA-DRB1` (`99%`), `SIRPA` (#185, `57%`), `HLA-DRA` (`98%`), `NR4A3` (#271, `70%`), `NOTCH2` (#464, `59%`) |
| **3** | 1,160 (15.60%) | **Classical Type 1 Dendritic Cells (`cDC1`, `XCR1^{hi} CADM1^{hi} ANPEP/CD13⁺ CLNK⁺ DNASE1L3⁺ BATF3⁺ WDFY4^{hi} LY75⁺`)** | Keep (Clean) | `XCR1` (`79%`), `CADM1` (`78%`), `ANPEP` (`CD13`, `50%`), `CLNK` (`49%`), `DNASE1L3` (`62%`), `SNX22` (`33%`), `SLAMF8` (`30%`), `RAB32` (`72%`), `ARHGAP18` (`79%`), `RAB7B` (`40%`), `SCIMP` (`27%`), `XIRP2`, `ME1`, `GLO1` | `ZNF366` (`DC-SCRIPT`, #40, `49%`), `LY75` (`DEC-205`, #49, `56%`), `BATF3` (#53, `62%`), `CIITA` (#59, `65%`), `WDFY4` (#73, `79%`), `SLA-DRB1` (`97%`), `ITGAE` (`CD103`, #109, `39%`), `ID2` (#160, `75%`), `FLT3` (#321, `76%`), `IRF8` (`88%`) |
| **4** | 317 (4.26%) | **Mature / Migratory Immunoregulatory DCs (`mregDC`, `CCR7^{hi} FSCN1^{hi} LAMP3^{hi} IDO1^{hi} IL4I1⁺ CD274/PD-L1⁺ PDCD1LG2/PD-L2⁺`)** | Keep (Clean) | `IDO1` (`78%` vs `2%`), `FSCN1` (`84%` vs `4%`), `S1PR1` (`62%` vs `3%`), `S1PR5` (`47%`), `TNFRSF11A` (`RANK`, `73%`), `CD200` (`47%`), `IL4I1` (`52%`), `PDCD1LG2` (`PD-L2`, `49%`), `CCL22` (`41%`), `CACNB3` (`35%`), `TMEM150C`, `SLCO5A1`, `ANKRD33B` | `STAT4` (#26, `43%`), `CCR7` (#40, `97%`), `IL21R` (#50, `73%`), `CD274` (`PD-L1`, #54, `58%`), `CD40` (#68, `86%`), `CD80` (#142, `71%`), `LAMP3` (`DC-LAMP`, #174, `94%`), `RELB` (#200, `77%`), `CD83` (#488, `96%`), `FLT3` (#571, `77%`) |

---

## Concise Cluster Rationale

* **C1 (`4,499`) — `pDC`**: Expresses core pDC master regulators (`TCF4/E2-2` `95%`, `IRF8` `99%`, `BCL11A` `92%`, `RUNX2` `81%`), endosomal TLRs (`TLR7` `66%`, `TLR9` `31%`), and porcine-specific pDC surface markers `CD4^{hi}` (`88%`), `CD8B⁺` (`52%`), and `NCR1 (NKp46)⁺` (`60%`; Mair et al. 2022 `PMID: 35371050`).
* **C2 (`1,460`) — `cDC2` + `DC3/MoDC`**: Captures canonical porcine lung `cDC2` (`FCER1A⁺ CD1.1⁺ CD207/Langerin⁺ MS4A2⁺ SIRPA⁺ KLF4⁺ NOTCH2⁺ MHC-II^{hi}`) together with inflammatory `CD163⁺ MRC1⁺ FCGR3A⁺ TLR4⁺ IL1B⁺ AREG⁺` `DC3 / MoDC` (Baillou et al. 2025 `PMID: 41235223`).
* **C3 (`1,160`) — `cDC1`**: Cross-presenting `XCR1^{hi} CADM1^{hi} ANPEP(CD13)⁺ CLNK⁺ SNX22⁺ DNASE1L3⁺ WDFY4^{hi} LY75(DEC-205)⁺ ITGAE(CD103)⁺ BATF3⁺ ZNF366⁺ FLT3⁺` `cDC1` (Auray et al. 2016 `PMID: 27798170`).
* **C4 (`317`) — `mregDC`**: Activated migratory `CCR7^{hi} FSCN1^{hi} LAMP3^{hi} CD83^{hi} CD40^{hi} CD80^{hi} S1PR1⁺` DCs upregulating tryptophan/phenylalanine catabolism (`IDO1` `78%`, `IL4I1` `52%`) and PD-1 ligands (`CD274/PD-L1` `58%`, `PDCD1LG2/PD-L2` `49%`).
