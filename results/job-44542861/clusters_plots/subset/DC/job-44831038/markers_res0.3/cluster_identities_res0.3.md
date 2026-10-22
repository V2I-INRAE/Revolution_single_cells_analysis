# DC Subclustering Identity Report (`job-44831038`, Resolution `0.3`)

- **Parent Run / Subset**: `job-44542861` $\rightarrow$ `DC` (`job-44831038`; cleaned sensitivity run, `k=50`; Parent Clusters 16 [`4,542` pDC] + 19 [`2,894` cDC])
- **Clustering Parameters**: `harmony_snn_res.0.3` (`7,436` cells across `5` clusters; `1,500` HVGs, `dims 1:10`, `k_param = 50`)
- **QC & Parameter Assessment**: **No coherent doublet-cluster signature identified in these exports; no reassignment indicated.**
  - **k=50 sensitivity result**: the same `7,436` cells resolve into `5` clusters. C3/C4/C5 reproduce the same broad k=100 cDC2-with-DC3-features, cDC1, and mregDC programs (24/25 top genes shared each; not identical partitions); the k=100 pDC cluster predominantly splits into C1 (`2,841`) + C2 (`1,756`) (`4,597` vs `4,499`, with additional membership changes per the clustree-review barcode comparison). Both C1 and C2 retain the listed pDC-associated program (`TCF4⁺ IRF7⁺ IRF8⁺ BCL11A⁺ TLR7⁺ TLR9⁺ CD4⁺ BLNK⁺ NCR1⁺ CMKLR1⁺`), supporting a shared broad pDC identity. The partition is compatible with within-pDC heterogeneity; these one-versus-rest exports do not establish whether it represents biological states, subtypes, or technical variation. No cluster shows a coherent foreign-lineage (T/B/plasma/monocyte/EC/epithelial) signature in these exports.
- **Master Index & References**: See `../../../NOTES_BARCODE_REASSIGNMENT_AND_CLEANUP.md` (Auray et al. 2016 `PMID: 27837108`; Maisonnasse et al. 2016 `PMID: 26530136`; Mair et al. 2022 `PMID: 35371050`; Baillou et al. 2025 `PMID: 41235223`).

---

## Summary Table (`harmony_snn_res.0.3`)

| Cluster | Cells (%) | Proposed Known Subset | QC / Doublet Flag | Key Named Top-25 Markers | Additional Canonical Markers (`markers_significant.csv`) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **1** | 2,841 (38.21%) | **Plasmacytoid DCs, State A (`TCF4⁺ CD4⁺ TCF4/BLNK/NCR1^{hi}`)** | Keep (Clean) | `NEGR1`, `CLSPN`, `NRXN1`, `CRHBP` (`62%`), `COBLL1` (#18, `80%`), `MGP` | `CD4` (#29, `86.5%`), `TCF4` (#51, `94%`), `NCR1` (#63, `60%`), `TLR3` (#82, `60%`); pDC program shared with C2 (see rationale) |
| **2** | 1,756 (23.61%) | **Plasmacytoid DCs, State B (`VPREB1⁺ MGP⁺ BLNK⁺`)** | Keep (Clean) | `VPREB1` (#1, `78%`), `TNF`, `MGP` (#18, `90%`), `BLNK` (#20, `98%`) | `NCR1` (#39, `58%`), `IRF7` (#66, `47%`), `CD4` (#85, `86%`); same pDC program as C1 at lower enrichment |
| **3** | 1,460 (19.63%) | **cDC2-like DCs with Inflammatory / DC3-like Features** | Keep (Clean) | `FCER1A` (#1, `32%`), `IL1B` (#3, `38%`), `CD1.1 (review)` (#4, `64%`), `CD207` (#5, `25%`), `S100A4` (#6, `70.4%`), `TLR4` (#7, `37%`), `MS4A2` (#9, `88%`), `FCGR3A` (#11, `79%`), `MRC1` (#14, `46%`), `CD163` (#16, `30%`) | `KLF4` (#71, `76%`), `SIRPA` (#184, `57%`), `SLA-DRB1 (review)` (`99%`) |
| **4** | 1,049 (14.11%) | **Classical Type 1 DCs (`cDC1`)** | Keep (Clean) | `XIRP2` (#2, `58%`), `CLNK` (#3, `54.5%`), `ANPEP` (`CD13`, #4, `54%`), `ME1` (#5, `45%`), `GABRA3`, `CADM1` (#10, `83%`), `XCR1` (#8, `87%`), `SNX22` (#12, `37%`), `DNASE1L3` (#18, `68%`) | `ZNF366` (#41, `53%`), `LY75` (#54, `60%`), `BATF3` (`68%`), `WDFY4` (`81%`), `ITGAE` (#98, `42.5%`), `FLT3` (`76%`) |
| **5** | 330 (4.44%) | **Mature / Migratory `mregDC`-like State** | Keep (Clean) | `TMEM150C` (#1, `72%`), `KREMEN1` (#2, `29%`), `IDO1` (#3, `76%`), `SLCO5A1`, `FSCN1` (#6, `81%`), `S1PR1` (#9, `59%`), `IL4I1` (#18, `49%`), `PDCD1LG2` (#19, `48%`) | `STAT4` (#27, `41%`), `CCR7` (#42, `96%`), `CD274` (#51, `57%`), `CD40` (#67, `84%`), `LAMP3` (#173, `93%`), `CD83` (#490, `94%` detection, log2FC `+1.25`) |

---

## Concise Cluster Rationale

* **C1 (`2,841`) + C2 (`1,756`) — pDC split by lower k**: Both retain the porcine pDC program (`TCF4` ~`94%`, `BLNK` ~`95–98%`, `NCR1` ~`60%`, `CD4` `86%`, `CMKLR1`, `BCL11A`, `IRF7`, `TLR7/9/3`) with small quantitative shifts (maximum one-versus-rest log2FC ≈ `1.81`; not C1-vs-C2 effect sizes). No subtype-level distinction is supported; treat as one pDC population unless cell-level or functional evidence says otherwise.
* **C3 (`1,460`) — cDC2 with DC3-like features**: Identical program to k=100 C2 (`FCER1A⁺ CD1.1⁺ CD207⁺ MS4A2⁺ FCGR3A⁺ TLR4⁺ S100A4⁺ MRC1⁺ CD163⁺ IL1B⁺`). Finer cDC2/DC3 identity unresolved.
* **C4 (`1,049`) — `cDC1`**: `XCR1^{hi} CADM1^{hi} ANPEP⁺ CLNK⁺ SNX22⁺ DNASE1L3⁺ BATF3⁺ LY75⁺ ITGAE⁺ ZNF366⁺ FLT3⁺`.
* **C5 (`330`) — `mregDC`-like**: `CCR7^{hi} FSCN1^{hi} IDO1^{hi} IL4I1⁺ PDCD1LG2⁺ CD274⁺ CD40⁺ S1PR1⁺ TMEM150C⁺` with `LAMP3` broadly detected and `CD83` at `94%` detection (log2FC `+1.25`).

> Note: pre-red-palette outputs (viridis dotplot + prior CSVs/heatmaps) are archived in `old/` in this directory; current plots use the red scale.
