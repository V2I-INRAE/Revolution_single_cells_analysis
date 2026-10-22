# NK_cells Subclustering Identity Report (`job-44831035`, Resolution `0.2`)

- **Parent Run / Subset**: `job-44542861` $\rightarrow$ `NK_cells` (`job-44831035`; Parent Cluster 3; `19,389` cells)
- **Clustering Parameters**: `harmony_snn_res.0.2` (`19,389` cells across `3` clusters; `2,000` HVGs, `dims 1:15`, `k_param = 100`)
- **QC & Reassignment Summary**: **`0` doublet-dominated clusters identified from this marker review.** The reassigned NK-like T/NKT populations do not reappear as clusters (`CD3E`/`CD4`/`FOXP3` have no significant rows anywhere; `IFNG` has no significant row). Residual T-lineage signal persists in C2 and needs flagging: `CD3G` is a significant row (#30, `log2FC = +1.91`, `44.6%` vs `15.2%`) alongside weak `CD3D` (#223, `+0.48`); `CD3E` has no significant row in any cluster, but absence from the filtered table does not establish absent expression. C2 additionally carries `PLVAP` (#20, `+2.45`, `43.9%` vs `11.3%`), worth checking for a broader foreign-lineage program — not sufficient alone to call doublets. **Keep all `3` clusters pending expression-level review of CD3/TCR coexpression; no re-subclustering prescribed.**
  - **Resolution comparison**: res `0.3` resolved the same `19,389` cells in `5` clusters. This res-`0.2` clustering is marker-compatible with consolidation into `3`: new C1 (`10,006`) matches old tissue-resident C1 (`8,732`) plus membership shifts; new C2 (`8,252`) matches old cytolytic C2 (`7,077`) plus membership shifts; new C3 (`1,131`) matches old quiescent C5 (`1,155`). The res-`0.3` immature/progenitor C3 (`1,216`) and basal C4 (`1,209`) programs do not reappear as separate clusters at res `0.2` (no significant `IL2RA`/`IL7R`/`SPOCK2`/`BTC`/`ANXA1` rows in new C3). Cross-run membership is barcode-unverified; counts are compared, not traced.
- **Master Index & References**: See `../../../NOTES_BARCODE_REASSIGNMENT_AND_CLEANUP.md` (Gerner et al. 2009 `PMID: 18601948`; Mair et al. 2016 `PMID: 27471504`; Herrera-Uribe et al. 2021 `PMID: 34249103`; Wiarda et al. 2026 `PMID: 41859083`).

---

## Summary Table (`harmony_snn_res.0.2`)

| Cluster | Cells (%) | Proposed Known Subset | QC / Doublet Flag | Selected Identity Markers (`markers_significant.csv`) |
| :--- | :--- | :--- | :--- | :--- |
| **1** | 10,006 (51.61%) | **Tissue-Resident / Cytokine-Responsive `GZMK⁺ CXCR6⁺ CD69⁺ EOMES⁺ NCR1⁺ KLRB1⁺` NK Cells** | Keep (pending CD3 review; no `CD3E` significant row) | `GZMK` (#1, `51.6%`, `log2FC = 3.81`), `CXCR6` (#2, `55.8%`), `EOMES` (#3, `49.7%`), `CD69` (#4, `71.0%`), `IL18RAP` (#6, `29.3%`), `IL18R1` (#9, `57.3%`), `NCR1` (`NKp46`, #15, `89.0%`), `KLRB1` (#17, `91.4%`), `TCF7` (#39, `27.1%`), `GZMA` (#87, `94.7%`); `GZMB`-depleted (#834), `GZMH`-depleted (#840) |
| **2** | 8,252 (42.56%) | **Mature Cytolytic Intravascular `CX3CR1⁺ GZMB^{hi} GNLY^{hi} KLF2^{hi}` NK Cells** | Keep w/ `CD3G` + `PLVAP` flag (pending expression review) | `CX3CR1` (#1, `33.0%`), `GZMH` (#2, `45.8%`), `S1PR5` (#4, `39.8%`), `GZMB` (#8, `68.3%`), `KLF2` (#9, `79.6%`), `GNLY` (#10, `91.0%`), `KLRG1` (#11, `33.7%`), `KLF3` (#14, `26.9%`), `PLVAP` (#20, `43.9%`), `ZNF683` (#28, `32.5%`), `CD8A` (#32, `45.3%`, porcine NK `CD8αα`), `ITGAM` (`CD11b`, #40, `60.3%`), `TBX21` (#69, `64.5%`); `CD3G` #30 (`44.6%` vs `15.2%`, flagged); `PRF1` (#244, `73.4%`); `NCR1`-depleted (#900, `31.0%`) |
| **3** | 1,131 (5.83%) | **Quiescent `BACH2⁺ ZEB1⁺ ZBTB16(PLZF)⁺ IKZF2(Helios)^{hi}` Innate NK / ILC Subset** | Keep (pending CD3 review; no `CD3E` significant row) | `SMYD3` (#2, `42.4%`), `KCNQ5` (#4, `33.2%`), `AFF3` (#10, `70.0%`), `BACH2` (#13, `40.5%`), `ZEB1` (#21, `58.4%`), `ZBTB16` (`PLZF`, #115, `48.2%`), `IKZF2` (`Helios`, #195, `82.8%`); `NCR1`: no significant differential enrichment (prefiltered `58.7%` vs `62.8%`, adj `p = 1`); `CD3D`-depleted (#1074, `39.3%`) |

---

## Concise Cluster Rationale

* **Mature NK Dichotomy (C1 = `10,006`, C2 = `8,252`)**: Textbook porcine/mammalian NK split: **C1** is tissue-resident/cytokine-responsive (`GZMK` #1, `CXCR6` #2, `EOMES` #3, `CD69` #4, `IL18R1/RAP` #9/#6, `NCR1` #15 `89.0%`, `KLRB1` #17 `91.4%`, `TCF7` #39, `GZMA` #87 `94.7%`, `TYROBP`/`FCER1G` adaptor rows) with `GZMB`/`GZMH` depletion; **C2** is mature intravascular cytolytic (`CX3CR1⁺ S1PR5⁺ KLF2^{hi} KLF3⁺ KLRG1⁺ GZMB^{hi} GZMH⁺ GNLY^{hi} PRF1⁺ CD8αα⁺ ITGAM/CD11b⁺ TBX21⁺`) with `NCR1` depletion (#900) — but `CD3G` (#30) and `PLVAP` (#20) are flagged in QC and need expression-level follow-up before this cluster is called clean. `CD3E` has no significant row in any cluster.
* **Quiescent Subset (C3 = `1,131`)**: `BACH2` (#13) / `ZEB1` (#21) / `ZBTB16` (#115) / `IKZF2` (#195 `82.8%`) with `SMYD3` (#2), `AFF3` (#10); count-compatible with res-`0.3` C5 (`1,155`). The res-`0.3` immature/progenitor (`IL2RA`/`IL7R`/`SPOCK2`) and basal (`BTC`/`ANXA1`) programs have no corresponding significant rows here and are not recovered as separate clusters at this resolution.
* **Reassignment check**: no cluster shows a coherent T-cell program; biological purity of the kept clusters is not established from this marker review alone (residual `CD3G` in C2, `PLVAP` in C2); direct CD3/TCR coexpression checks are needed.

> Note: pre-red-palette outputs (viridis dotplot + prior CSVs/heatmaps) are archived in `old/` in this directory; current plots use the red scale.
