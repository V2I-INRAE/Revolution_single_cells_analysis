# Master Subset Annotation Index, Doublet/Reassignment Ledger & Handoff Guide (`job-44542861`)

This document is the central index for the 12 subclustered coarse compartments of `job-44542861` (**238,627 recorded original cells** across 23 porcine Control and EVLP lung samples). It links the original identity reports, designated doublet removals and planned lineage reassignments. **Cleanup/reassignment counts below are expected outcomes, not executed operations.** The reviewed regrouping and parameter policy is in Sections 4–5.

**Current companion analyses:** [nested marker panels](NOTES_REVIEWED_MARKER_PANELS.md) and [sequential regrouping proposal](NOTES_CLUSTER_REGROUPING_PROPOSAL.md). Original identity reports remain unchanged for provenance; their fine labels and earlier parameter recommendations are not automatically endorsed. The annotation record `data/job-44542861/annotated_clustered.INFO` and all 12 cluster/sample summaries support the corrected total.

---

## 1. Directory of Subset Reports & Literature Dossiers

### A. Subset Identity Reports (`12 / 12` Complete)

“Clean Cells Kept” is the expected count after the listed cleanup/transfers, with the NK-derived T destination scenario specified below. Optional mast/ILC2-like exports are not applied in these totals.

| # | Coarse Subset | Subcluster Job | Selected Res | Total Cells | Clean Cells Kept | Doublets Removed | Reassigned Out / In | Subset Identity Report Path |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | **`B_cells`** | `job-44614229` | `0.4` (9 cl.) | `26,790` | `26,790` | `0` | `0` / `0` | `B_cells/job-44614229/markers_res0.4/cluster_identities_res0.4.md` |
| 2 | **`CD4_T_cells`** | `job-44614236` | `0.4` (6 cl.) | `25,888` | `20,774` | `0` | `-6,783` (C3, C4) / `+1,669` | `CD4_T_cells/job-44614236/markers_res0.4/cluster_identities_res0.4.md` |
| 3 | **`CD8_T_cells`** | `job-44614237` | `0.4` (6 cl.) | `20,552` | `38,518`* | `0` | `0` / `+17,966`* | `CD8_T_cells/job-44614237/markers_res0.4/cluster_identities_res0.4.md` |
| 4 | **`DC`** | `job-44614230` | `0.3` (4 cl.) | `7,436` | `7,436` | `0` | `0` / `0` | `DC/job-44614230/markers_res0.3/cluster_identities_res0.3.md` |
| 5 | **`Endothelial`** | `job-44614231` | `0.3` (6 cl.) | `14,154` | `13,836` | `318` (C6) | `0` / `0` | `Endothelial/job-44614231/markers_res0.3/cluster_identities_res0.3.md` |
| 6 | **`Epithelial`** | `job-44614232` | `0.4` (11 cl.) | `10,029` | `9,426` | `603` (C7) | `0` / `0` | `Epithelial/job-44614232/markers_res0.4/cluster_identities_res0.4.md` |
| 7 | **`Granulocytes`** | `job-44614233` | `0.3` (6 cl.) | `23,037` | `23,037` | `0` | `0` (`1,002` Mast opt.) | `Granulocytes/job-44614233/markers_res0.3/cluster_identities_res0.3.md` |
| 8 | **`Macrophage_Monocytes`** | `job-44614234` | `0.4` (10 cl.) | `42,214` | `41,524` | `690` (C9, C10) | `0` / `0` | `Macrophage_Monocytes/job-44614234/markers_res0.4/cluster_identities_res0.4.md` |
| 9 | **`NK_cells`** | `job-44614238` | `0.4` (7 cl.) | `24,342` | `19,389` | `0` | `-4,953` (C3, C4) / `0` | `NK_cells/job-44614238/markers_res0.4/cluster_identities_res0.4.md` |
| 10 | **`Stromal`** | `job-44614235` | `0.4` (7 cl.) | `11,398` | `11,079` | `319` (C6) | `0` / `0` | `Stromal/job-44614235/markers_res0.4/cluster_identities_res0.4.md` |
| 11 | **`unconventional_T_cells`** | `job-44614239` | `0.4` (8 cl.) | `20,852` | `22,374`* | `0` | `0` (`399` ILC2 opt.) / `+1,522` | `unconventional_T_cells/job-44614239/markers_res0.4/cluster_identities_res0.4.md` |
| 12 | **`T_cells_unresolved`** | `job-44614240` | `0.4` (5 cl.) | `11,935` | **`0` (Dissolved)** | `2,514` (C2) | `-9,421` (C1, C3, C4, C5) | `T_cells_unresolved/job-44614240/markers_res0.4/cluster_identities_res0.4.md` |
| **Total** | **12 original → 11 expected containers** (scenario) | — | — | **`238,627`** | **`234,183` (`98.14%`)** | **`4,444` (`1.86%`)** | **`21,157` planned outgoing transfers** | — |

*\*Scenario, not a resolved lineage decision: the displayed CD8 count of 38,518 assumes all 4,953 NK C3/C4 cells enter CD8; unconventional-T then contains 22,374. If all enter unconventional-T instead, CD8=33,565 and unconventional-T=27,327. Before that decision, there are 33,565 specified-recipient CD8 cells, 22,374 unconventional-T cells, and 4,953 T-lineage cells pending placement. Do not call them confirmed NKT/MAIT/gamma-delta cells or count them twice. Eleven final containers assumes this routing is resolved and optional mast/ILC2-like exports are not made; those choices change container counts, not the retained global total.*

### B. Associated Literature & Parent Annotation Dossiers
* **Porcine Lung Macrophage & Monocyte Literature Review**: `results/literature_review/porcine_lung_macrophage_monocyte_subsets/literature_report.md` (with `paper_summary_table.csv` and `search_log.md`)
* **Whole-Lung Parent (`res 0.7`) Literature Scratchpad**: `results/literature_review/SCRATCHPAD.md`
* **Whole-Lung Parent (`res 0.7`) Cluster Identities**: `results/job-44542861/clusters_plots/markers/0.7/cluster_identities_res0.7.md`

---

## 2. All Identified Doublet Clusters (`REMOVE (Doublet)` — `4,444` Cells Total)

*(Note: `B_cells`, `CD4_T_cells`, `CD8_T_cells`, `DC`, `Granulocytes`, `NK_cells`, and `unconventional_T_cells` had **0** doublet clusters.)*

| Subset | Subcluster Job | Resolution & Cluster ID | Cells (% of Subset) | Doublet Composition & Key Markers | Exact Extraction Rule (`subclustered.rds`) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **`Endothelial`** | `job-44614231` | `res 0.3` — **Cluster `6`** | `318` (`2.25%`) | **Endothelial + Leukocyte Doublets** (`PTPRC/CD45` `86.2%`, `LCP1`, `CD53`, `HCST`, `TYROBP`, `LAPTM5`, `CD3E`, `CCL5`, `S100A12` + pan-EC RNA) | `harmony_snn_res.0.3 == "6"` |
| **`Epithelial`** | `job-44614232` | `res 0.4` — **Cluster `7`** *(or `res 0.3` — **Cluster `6`**; identical barcodes)* | `603` (`6.01%`) | **Epithelial + T/NK Leukocyte Doublets** (`PTPRC/CD45` `75.1%`, `CD2` `42.6%`, `CD3D` `38.1%`, `CD3E` `35.0%`, `GZMA` `25.0%`, `CCL5` `44.3%` + diluted `EPCAM`/`SFTPA1`/`SCGB1A1`) | `harmony_snn_res.0.4 == "7"` *(or `res.0.3 == "6"`)* |
| **`Macrophage_Monocytes`** | `job-44614234` | `res 0.4` — **Cluster `10`** | `245` (`0.58%`) | **Multi-Lineage Doublet Soup** (`COL1A1/2/3⁺` Fibroblast + `SFTPC/SCGB1A1⁺` Epithelial + `JCHAIN⁺` Plasma + `GNLY/CD247⁺` T/NK) | `harmony_snn_res.0.4 == "10"` |
| **`Macrophage_Monocytes`** *(Borderline EC-adherent complex)* | `job-44614234` | `res 0.4` — **Cluster `9`** | `445` (`1.05%`) | **Patrolling Non-Classical Monocyte (`FCGR3A^{hi} CD14^{lo} CSF1R⁺ MAFB⁺ POU2F2⁺ KLF2⁺ HES1⁺ SIGLEC5⁺ ADGRG1^{hi}`) + Capillary EC (`PLVAP⁺ KIT⁺ PECAM1⁺ HPGD⁺ CDH5⁺`) Doublet / Complex** (`MARCO⁻ TIMD4⁻ VSIG4⁻ MRC1⁻`) | `harmony_snn_res.0.4 == "9"` |
| **`Stromal`** | `job-44614235` | `res 0.4` — **Cluster `6`** | `319` (`2.80%`) | **Stromal + Leukocyte (T/NK/Myeloid) Doublets** (`PTPRC/CD45` `55.2%`, `SRGN` `77.7%`, `LAPTM5` `59.9%`, `CXCR4` `54.9%`, `RGS1` `53.9%`, `FCGR3A` `38.2%`, `CCL5` `37.0%`, `CD2` `33.5%`, `CD3D` `27.6%` + diluted `COL1A1`/`COL3A1`/`DCN`) | `harmony_snn_res.0.4 == "6"` |
| **`T_cells_unresolved`** | `job-44614240` | `res 0.4` — **Cluster `2`** | `2,514` (`21.06%`) | **T-Cell + Alveolar Macrophage Doublets** (`APOE` `28.0%`, `FCN2` `25.8%`, `BTC` `28.1%`, `CD163` `35.0%`, `LGALS3` `32.9%`, `C3` `26.7%`, `CD9` `33.2%`, `CYBB` `27.0%`, `LGMN` `26.7%`, `ANXA1` `42.6%`, `NCF2` `25.0%` + diluted `CD3D`/`CD3E`/`CD3G`/`CD2`) | `harmony_snn_res.0.4 == "2"` |

---

## 3. Planned Lineage Reassignments (`21,157` Mandatory Outgoing + `1,401` Optional)

Of the 21,157 planned outgoing transfers, **16,204 have specified recipient compartments** and **4,953 NK-derived T cells require recipient-lineage review**. Historical descriptive identities/marker summaries below are retained as working annotations, not receptor-validated subtype claims. CD247 alone does not establish T lineage; cytotoxic/NK-receptor programmes alone do not distinguish NK-like alpha-beta from gamma-delta T cells. Consult the companion marker review for limitations. Optional transfers retain all cells.

| Source Subset | Subcluster Job | Resolution & Cluster ID | Cells (% of Subset) | Working Biological Identity & Key Markers | Target Coarse Subset |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **`CD4_T_cells`** | `job-44614236` | `res 0.4` — **Cluster `3`** | `3,748` (`14.48%`) | **Naive / Central Memory `CD8⁺` T Cells** (`CD8B` `95.0%`, `CD8A` `71.0%`, `CD4` `0.8%`, `SELL` `69.6%`, `CCR7` `53.0%`, `LEF1`, `TCF7`) | **`CD8_T_cells`** |
| **`CD4_T_cells`** | `job-44614236` | `res 0.4` — **Cluster `4`** | `3,035` (`11.72%`) | **`GZMK⁺ CCL5⁺` Effector-Memory `CD8⁺` T Cells** (`CD8A` `74.4%`, `CD8B` `62.2%`, `CD4` `4.1%`, `CCL5` `90.9%`, `GZMK` `61.3%`, `KLRK1` `57.2%`, `CXCR3`) | **`CD8_T_cells`** |
| **`NK_cells`** | `job-44614238` | `res 0.4` — **Cluster `3`** | `2,892` (`11.88%`) | **Activated `CD3⁺ CD8αα⁺ FCGR3A(CD16)⁺ GNLY^{hi} CSF1⁺ MHC-II⁺` Cytotoxic NK-like T / NKT Cells** (`CD3D` `75.1%`, `CD3E` `70.1%`, `CD3G` `63.9%`, `CD247` `65.9%`, `CD8A` `63.3%`, `CD8B⁻`, `TYROBP⁻` `8.3%`) | **`CD8_T_cells`** *(or `unconventional_T_cells`)* |
| **`NK_cells`** | `job-44614238` | `res 0.4` — **Cluster `4`** | `2,061` (`8.47%`) | **Resident-Memory `CD3⁺ CD8αα⁺ TCF7⁺ CRTAM^{hi} ZNF683⁺ GZMM⁺` Cytotoxic T / NKT Cells** (`CD3D` `80.5%`, `CD3E` `75.3%`, `CD3G` `72.1%`, `CD247` `72.3%`, `CD8A` `60.8%`, `CD8B⁻`, `TYROBP⁻` `7.7%`) | **`CD8_T_cells`** *(or `unconventional_T_cells`)* |
| **`T_cells_unresolved`** | `job-44614240` | `res 0.4` — **Cluster `1`** | `4,017` (`33.66%`) | **Parenchymal Tissue-Resident / Effector-Memory `CD8⁺` T Cells ($T_{EM}/T_{RM}$)** (`GZMK` `31.9%`, `GZMA` `57.8%`, `CXCR6` `39.2%`, `CCR5` `46.6%`, `CRTAM` `45.1%`, `CD96` `32.7%`, `JAML` `26.1%`, `KLRK1` `50.9%`, `CD8A` `56.9%`, `CCL5` `86.9%`, `CD69` `62.8%`, `RUNX3` `41.0%`) | **`CD8_T_cells`** |
| **`T_cells_unresolved`** | `job-44614240` | `res 0.4` — **Cluster `3`** | `2,213` (`18.54%`) | **Intravascular Cytolytic Effector `CD8⁺` T Cells ($T_{EMRA}/T_{EFF}$)** (`ADGRG1` `70.2%`, `GZMB` `64.8%`, `GNLY` `61.3%`, `GZMH` `47.9%`, `PRF1` `50.0%`, `NKG7` `63.7%`, `FCGR3A` `87.0%`, `TBX21` `47.7%`, `KLF2` `78.7%`, `KLF3` `37.7%`, `S1PR1` `34.0%`, `IL10` `26.3%`, `CD8A` `54.0%`, `CD2` `82.0%`, `CD3D` `64.0%`) | **`CD8_T_cells`** |
| **`T_cells_unresolved`** | `job-44614240` | `res 0.4` — **Cluster `4`** | `1,669` (`13.98%`) | **Helper / Effector-Memory `CD4⁺` T Cells** (`CD40LG` `33.9%`, `CD28` `43.3%`, `MAF` `47.9%`, `IL7R` `50.4%`, `IL6R` `27.4%`, `IL18R1` `33.1%`, `TCF7` `33.0%`, `CD5` `65.8%`, `CD6` `38.0%`, `GPR183` `68.0%`, `SOCS3` `46.3%`, `CD8A⁻` `18%`, `GNLY⁻` `2%`, `GZMB⁻` `4%`) | **`CD4_T_cells`** |
| **`T_cells_unresolved`** | `job-44614240` | `res 0.4` — **Cluster `5`** | `1,522` (`12.75%`) | **`CD2⁻ WC1⁺` Porcine $\gamma\delta$ T Cells** (`SOX13` `44.5%`, `RHEX` `33.4%`, `BLK` `37.4%`, `GATA3` `56.5%`, `MAF` `46.4%`, `ZBTB16` `28.2%`, `ID2` `86.0%`, `IL2RA` `34.2%`, `IL7R` `56.7%`, `AQP3` `51.9%`, `CGA` `39.4%`, `AHR` `32.0%`, `CD2⁻` `13%`, `LCK⁻` `4%`, `CD6⁻` `6%`) | **`unconventional_T_cells`** |
| **`Granulocytes`** *(Optional)* | `job-44614233` | `res 0.3` — **Cluster `4`** | `1,002` (`4.35%`) | **Resident Pulmonary Mast Cells** (`KIT` `46.9%`, `MS4A2` `38.4%`, `FCER1A` `30.3%`, `GATA2`, `PTGS1`, `HPGD`, `RAB27B`, `MLPH`, `VWA5A`) | Keep in `Granulocytes` (as `Mast_cells`) or promote to coarse `Mast_cells` |
| **`unconventional_T_cells`** *(Optional)* | `job-44614239` | `res 0.4` — **Cluster `8`** | `399` (`1.91%`) | **Pulmonary Group 2 Innate Lymphoid Cells (`ILC2`)** (`CD3D/E/G⁻`, `ID2^{hi}` `90%`, `GATA3^{hi}` `89%`, `KIT^{hi}` `84.7%`, `MAF^{hi}` `83%`, `IL2RA^{hi}` `64.2%`, `HLF^{hi}` `60.4%`, `CD2⁺` `60.2%`, `RORA⁺` `60%`, `ZBTB16⁺` `57%`, `HOPX^{hi}` `54.1%`, `IL7R⁺` `52%`, `IL1RL1/ST2⁺` `42%`, `PTGDR2/CRTH2⁺` `26.1%`) | Keep in `unconventional_T_cells` (as `ILC2`) or promote to coarse `ILC2` |

---

## 4. Reviewed regrouping and graph-parameter policy

**This section supersedes the earlier promises of specific lineage splits with smaller `k_param` and the untested assertion of an established `k=50` optimum.** The evidence supports candidate groupings and sensitivity tests, not guaranteed splits.

Read the [full sequential regrouping proposal](NOTES_CLUSTER_REGROUPING_PROPOSAL.md) and [reviewed marker-panel tables](NOTES_REVIEWED_MARKER_PANELS.md) before adopting the historical identity labels. In particular, broad identity, original fine-cluster provenance, and state/programme should be separate annotations. A clustree plateau or shared parent on one original graph is not donor/resampling validation, subtype equivalence, or a developmental trajectory.

### Baseline versus sensitivity settings

`R/config.R` currently sets `k_param=100`; no subset-optimal k has been demonstrated. Across the original 7,436–42,214-cell subsets, 100 neighbours correspond to roughly **0.24%–1.35%** of cells, but that fraction alone does not determine biological resolution. Smaller k may reveal structure **or fragment it**; larger k may oversmooth. A rare cluster being smaller than twice k is a reason to inspect sensitivity, not proof of mixed lineage.

For graph refits, start from the reported subset HVGs/dimensions and k=100; compare k=50 with dimensions/features fixed and resolutions 0.2/0.3/0.4. Reserve targeted k≈30 tests for specific, cell-level-supported questions. Do not change k merely to obtain fewer reporting groups. Metadata roll-ups do not require new clustering.

| Subset | Reported original settings (HVGs; dims; k; res) | Pragmatic annotation / refit recommendation |
|---|---|---|
| B_cells | 2,000; 1:20; 100; 0.4 | Merge mature-plasma C3+C4 as a candidate curated identity. Broad ASC C3+C4+C7+C8 is a reporting umbrella with secondary states retained. Existing res0.3 is the first comparison; no k change required. |
| CD4_T_cells | 2,000; 1:20; 100; 0.4 | After C3/C4 transfer out and unresolved-T C4 transfer in, refit the 20,774-cell pool before graph-dependent conclusions. C1+C5 can share a broad helper-memory/activated umbrella; C5 is not resolved Treg/Tfh. |
| CD8_T_cells | 2,000; 1:20; 100; 0.4 | Refit after incoming transfers; 33,565 cells before NK-derived T placement, 38,518 if all enter here. No forced original six-group mapping or C3+C4 merge: C3 lineage remains provisional. |
| DC | 1,500; 1:20; 100; 0.3 | Keep four operational groups. Do not force C2 into cDC2/DC3/MoDC subtypes or infer C4 precursor lineage from a graph parent. Targeted sensitivity only if warranted. |
| Endothelial | 2,000; 1:15; 100; 0.3 | Remove C6; retain five vascular identities and separate activation scores. Refit if cleaned graph/neighbourhood analyses are needed; no obligatory venous subdivision. |
| Epithelial | 2,000; 1:20; 100; 0.4 | Remove C7. Broad ciliated C1+C5 and secretory C2+C3 are recommended; tuft/brush-like C8+C9 is an optional umbrella. Preserve secondary programmes and rare lineages. |
| Granulocytes | 500; 1:10; 100; 0.3 | Broad neutrophils C1+C2+C5 with inflammatory/IFN/granule-rich states; minimal merge C1+C2. Keep eosinophils, basophils and mast identity distinct. No new graph required for reporting. |
| Macrophage_Monocytes | 2,000; 1:20; 100; 0.4 | Remove C9/C10. Broad PAM C1+C2+C3 and monocyte/MoMac C5+C6+C8 with state reporting; retain C4 and niche-provisional C7 separately. k=50 is sensitivity, not a proven optimum. |
| NK_cells | 2,000; 1:15; 100; 0.4 | Move C3/C4 out to T-lineage review; refit for new boundary decisions. C5 remains low-information, not a validated subtype or automatically C2. C6/C7 identities are provisional NK/ILC-like states. |
| Stromal | 2,000; 1:15; 100; 0.4 | Remove C6. Default C1+C4 matrix/adventitial-like umbrella; C3+C5 optional alveolar/niche roll-up, with fine distinctions preserved. Keep contractile C2 and mural C7 separate. No guaranteed pericyte/VSMC or myofibroblast/ASM split with smaller k. |
| unconventional_T_cells | 2,000; 1:20; 100; 0.4 | Refit after adding unresolved-T C5. Broad CD2-low gamma-delta-like C1–C6 with state reporting; C7 and ILC2-like C8 distinct. Do not force low-signal C5 into another fine cluster. |
| T_cells_unresolved | 2,000; 1:20; 100; 0.4 | Dissolve after C2 removal and other-cluster transfers; no new graph for this container. |

**Stromal caution:** C1 and C4 coexist in all five baseline-control samples; C1 cannot simply be described as an EVLP-only replacement for C4. Shared matrix markers support a broad umbrella, not a proven single subtype. C5's niche-associated programme warrants preservation. Details and control counts are in the proposal.

**Implementation status:** the inspected `R/clusters/subcluster.R` and `R/utils/cli.R` support `--hvg` and `--dims`, not `--k-param`. A per-run k override would require a separately approved code change. Do not alter the shared configuration globally for one subset. No such implementation or job submission has been performed.

---

## 5. Handoff guidance — cleanup, hierarchical annotation, then selective refit

**Plan only: this document does not certify an executed cleaned RDS.** Original subset objects preserve parent cell barcodes according to the annotation records. Recheck this at execution. The detailed source/job/resolution index and validation criteria are in [the regrouping proposal](NOTES_CLUSTER_REGROUPING_PROPOSAL.md).

### A. Barcode cleanup and reassignment

1. Extract removal and transfer barcode sets from the original subset jobs using their selected resolutions in Sections 2–3. Never match a bare cluster number across jobs or resolutions.
2. Assert source ownership, unique barcode keys, disjoint retained/removed sets, and exactly one curated destination for every retained cell. Removal takes precedence; transfers conserve cells.
3. Remove the designated **4,444** doublet barcodes. Expected retention is **234,183 of 238,627 cells (98.14%)**.
4. Apply the **16,204** transfers with specified recipient compartments; route the additional **4,953** NK-derived T cells through explicit lineage review before choosing CD8 versus unconventional-T. Leaving NK is planned, but the recipient subtype is not established. Preserve a pending-review bucket if necessary rather than force or double-count them.
5. Preserve optional mast (1,002) and ILC2-like (399) cells as distinct identities, whether exported to dedicated compartments or retained in their original containers. Global retention is unchanged by these optional transfers.
6. Save any future cleaned object under a new name, without overwriting the original; verify barcodes, counts, metadata and source provenance after reopening.

### B. Hierarchical annotations without compulsory reclustering

Build a barcode lookup retaining `orig_subset`, `orig_subcluster_job`, `orig_resolution`, `orig_subcluster_id`, `curated_compartment`, broad/fine curated identities, programme/state annotations, confidence, QC status and reassignment status.

Use the proposal rather than blindly copying every historical biological label. Prioritise B C3+C4, epithelial C2+C3 and C1+C5, granulocyte C1+C2, and the broad stromal C1+C4 umbrella. Preserve source labels and secondary states. Broader ASC/PAM/neutrophil/gamma-delta roll-ups are reporting choices, not proof that constituent populations are identical. Purely metadata-level regrouping does not require renormalisation or new UMAPs.

### C. Selective graph refits and validation

For receiving CD4/CD8/unconventional pools, refit features, scaling/PCA, integration, neighbours and clustering before graph-dependent analyses. Refit NK for new boundaries after outgoing transfers; refit other cleaned compartments if their graph is needed downstream. Do not concatenate old embeddings or graphs. Existing retained labels may still support explicitly provisional metadata-level analyses.

Start with the bounded k=100 versus k=50 comparison in Section 4, retaining reported dimensions/features initially. Evaluate closest-cluster multigene separation, QC effects, donor replication, seed/graph sensitivity and held-out donor performance where feasible. Treat repeated samples/timepoints from a donor as dependent observations. Low detection is not biological absence; absence of a five-gene panel is not a removal criterion.

Keep state-specific composition/DE alongside broad counts so merging cannot hide EVLP-associated responses. Recompute marker panels against cleaned curated labels. Use cell-cycle/stress-feature sensitivity analyses when justified; do not automatically regress away condition biology. A clustree is not evidence of temporal cell conversion.

**Expected final global count remains 234,183.** Eleven compartments assumes the unresolved-T container is dissolved, NK-derived T placement is resolved, and mast/ILC2-like populations remain inside existing containers; optional dedicated compartments or a temporary pending-review bucket alter the container count, not cell retention.

---

## 6. Master Literature Reference Index (By Lineage)

**Historical reference index retained for provenance, not newly literature-verified in this documentation update.** The preceding abstract-level review identified attribution/PMID problems and unsupported marker-specific interpretations in some entries, including the “Shin et al.” attribution for PMID 39875369, the Fairbairn journal/DOI, and the Schupp PMID (the verified atlas is **34030460**, not 34097468). See [verified anchors and qualifications](NOTES_REVIEWED_MARKER_PANELS.md). Do not use these historical claims alone to justify a subtype, niche, or merger; original reference entries have not been comprehensively corrected here.

### A. Porcine & Mammalian Monocytes, Alveolar Macrophages (PAMs), IMs, PIMs & DCs
1. **Álvarez et al. (2023)** *Cells* (`PMID: 37626913`, `PMC10453229`, DOI: `10.3390/cells12162103`) — Comprehensive review of porcine monocyte/macrophage markers; establishes PAMs as `CD203a(ENPP1)^{hi} CLEC12B⁺ CD206(MRC1)^{hi} CD169(SIGLEC1)^{hi} CD163^{hi} CD172a(SIRPA)^{hi} FCGR1A^{hi} CD14^{lo} CD16^{lo}` (`Macrophage_Monocytes` C1–C3).
2. **Poderoso et al. (2020)** *Dev Comp Immunol* (`PMID: 32535044`, DOI: `10.1016/j.dci.2020.103767`) — Demonstrates `CLEC12B` is specifically expressed on porcine alveolar macrophages (`Macrophage_Monocytes` C1–C3) and absent on blood monocytes.
3. **Shin et al. (2025)** *Nat Commun* (`PMID: 39875369`, `PMC11775223`, DOI: `10.1038/s41467-024-54676-2`) — `306,052`-cell porcine BALF scRNA-seq atlas; defines steady-state `MRC1-NREP^{high}` (`PON1⁺ PPBP⁺ CD9⁺`) PAMs (`Macrophage_Monocytes` C1), lipid-associated PAMs (`C2`), and `PDCD1LG2^{high}` macrophages (`C7`).
4. **Chen et al. (2025)** *Front Immunol* (`PMID: 40242774`, `PMC12000058`, DOI: `10.3389/fimmu.2025.1574120`) — `114,881`-cell pig macrophage atlas across 49 organs; validates `PLSCR1` + `FCGR1A` as lung-specific PAM markers (`C1–C3`), `FABP5/LPL/GPNMB/SPP1` lipid-remodeling module (`C2`), and `C4BPA` in pig macrophages (`C7`).
5. **Zhang et al. (2026)** *Virulence* (`PMID: 41843736`, `PMC13014559`, DOI: `10.1080/21505594.2026.2646800`) — `59,062`-cell whole piglet lung scRNA-seq atlas; distinguishes `S100A8/9⁺ CD9⁺` AMs (`C1–C3`), `MARCO⁺ CD163⁺ C1QA/B/C⁺ HMOX1⁺` IMs (`C4`, `C7`), and `CD14⁺ FCN1⁺` monocytes (`C5`, `C6`).
6. **Ke et al. (2023)** *Microbiol Spectr* (`PMID: 37078873`, `PMC10269733`, DOI: `10.1128/spectrum.04753-22`) — Shows `MARCO` is constitutively expressed on primary porcine PAMs and restricts PRRSV (`Macrophage_Monocytes` C7).
7. **Zong et al. (2025)** *Virulence* (`PMID: 40944943`, `PMC12459362`, DOI: `10.1080/21505594.2025.2560894`) — Identifies `SLAMF7^{high}` porcine BAL macrophage subpopulation (`Macrophage_Monocytes` C7).
8. **Muir et al. (2024)** *PLoS Pathog* (`PMID: 39024231`, `PMC11257366`, DOI: `10.1371/journal.ppat.1011910`) — Spectral flow and scRNA-seq of Babraham pig BAL vs. PBMC myeloid and lymphoid cells (`Macrophage_Monocytes`, `CD4_T_cells`, `CD8_T_cells`).
9. **Maisonnasse et al. (2016a,b), Bordet et al. (2018), Crisci et al. (2020)** (`PMID: 26530136`, `PMID: 27992536` / `PMC5167224`, `PMID: 29977043` / `PMC6033929`, `PMID: 33365028` / `PMC7750501`) — Dissect the porcine lung DC/macrophage network (`cDC1`, `cDC2`, `moDC`, `moMacro`, `AM`, and `AM-like` `PIM`), showing a `10:1` AM:PIM ratio and near-identical bulk RNA-seq profiles between sorted AMs and PIMs (`DC`, `Macrophage_Monocytes`).
10. **Auray et al. (2016) & Baillou et al. (2025)** (`PMID: 27798170`; *Front Immunol* `PMID: 41235223` / `PMC12606228`, DOI: `10.3389/fimmu.2025.1639553`) — Define porcine `XCR1⁺ CADM1⁺` `cDC1`, `FCER1A⁺ CD1.1⁺` `cDC2`, `CD163⁺` `DC3`, and `CD4^{hi} TCF4⁺ IRF8⁺` `pDC` (`DC` C1–C4).
11. **Fairbairn et al. (2013)** *J Leukoc Biol* (`PMID: 23667115`, DOI: `10.1189/jlb.0113030`) — Establishes the `CSF1`-driven maturation series from `CD14^{hi} CD163^{lo}` (`Macrophage_Monocytes` C6) to `CD14⁺ CD163^{hi}` (`C5`, `C8`) porcine monocytes.
12. **Gouin et al. (2023)** *Front Immunol* (`PMID: 37465668`, `PMC10351384`, DOI: `10.3389/fimmu.2023.1142228`) — Human `0h, 4h, 10h` EVLP lung scRNA-seq; shows `CD14⁺ FCN1⁺ VCAN⁺` `cMos` and `MoMacs` drive `IL1B/IL1R1/IL1R2/AREG/EREG` expression during EVLP (`Macrophage_Monocytes` C5, C6).
13. **Mulder et al. (2021)** *Immunity* (`PMID: 34297922`, DOI: `10.1016/j.immuni.2021.07.007`) — MoMac-VERSE cross-tissue atlas defining conserved `TIMD4⁺ MARCO⁺ VSIG4⁺ C1Q⁺` resident macrophages, `C1Q^{hi} MAF⁺` interstitial macrophages, and `IL1B⁺ AREG⁺ EREG⁺` `MoMacs` (`Macrophage_Monocytes` C4, C5, C7).

### B. Porcine T Cells ($\alpha\beta$, $\gamma\delta$), NK Cells & Innate Lymphoid Cells (`ILCs`)
14. **Herrera-Uribe et al. (2021)** *Front Genet* (`PMID: 34249103`, `PMC8261551`, DOI: `10.3389/fgene.2021.689406`) — Bulk and single-cell reference atlas of porcine PBMCs; defines `CD4⁺ CD8α⁻` naive vs. `CD4⁺ CD8α⁺` activated/memory porcine `CD4⁺` T cells (`CD4_T_cells` C1 vs. C2), `CD2⁻ WC1⁺` vs. `CD2⁺` $\gamma\delta$ T cells (`unconventional_T_cells` C1–C6 vs. C7), and `CD21⁺`/`CD21⁻` B cells and ASCs (`B_cells`).
15. **Wiarda et al. (2022) & Wiarda et al. (2026)** *Life Sci Alliance* (`PMID: 35995567`, `PMC9396248`, DOI: `10.26508/lsa.202201442`) & *Front Immunol* (`PMID: 41859083`, `PMC12995613`, DOI: `10.3389/fimmu.2026.1704257`) — Single-cell atlases of porcine intestinal, primary, and secondary lymphoid organs; define `SELL^{hi} ID3⁺ ITGB7⁺` $\gamma\delta$ T cells (`unconventional_T_cells` C1), `CD2⁺ CD3E⁻ ID2^{hi} KIT^{hi} IL7R⁺` porcine tissue ILCs (`unconventional_T_cells` C8 `ILC2`), and `AICDA⁺ RGS13⁺` GC B cells (`B_cells` C9).
16. **Rodríguez-Gómez et al. (2016, 2019)** *Dev Comp Immunol* (`PMID: 26920461`) & *Front Immunol* (`PMID: 30915070`, `PMC6421308`, DOI: `10.3389/fimmu.2019.00396`) — Establish that `~80–90%` of porcine `CD2⁻` $\gamma\delta$ T cells in blood and lung are `GATA3^{hi} CD27⁺ Perforin⁻` (`unconventional_T_cells` C1–C6), whereas `EOMES⁺` $\gamma\delta$ T cells (`~14%` in lung) are `CD2⁺ CD27⁺` (`unconventional_T_cells` C7).
17. **Gu et al. (2022) & Zhang et al. (2024)** *Cell Rep* (`PMID: 35793622`, `PMC9704770`, DOI: `10.1016/j.celrep.2022.111050`) & *Front Immunol* (`PMID: 38384475`, `PMC10879363`, DOI: `10.3389/fimmu.2024.1339787`) — Single-cell atlases of porcine thymopoiesis and peripheral T cells; confirm `SOX13, BLK, GATA3, ETV5, CD163L1(WC1)` as hub regulators of `CD2⁻` $\gamma\delta$ T cells (`unconventional_T_cells` C1–C6) and `IKZF2` in `CD2⁺` $\gamma\delta$ T cells (`C7`).
18. **Takamatsu et al. (2006) & Le Page et al. (2022)** *Anim Health Res Rev* (`PMID: 17389056`) & *Dev Comp Immunol* (`PMID: 34863955`, DOI: `10.1016/j.dci.2021.104330`) — Reviews of porcine unconventional T cells, documenting `CD80/CD86` and MHC-II upregulation on activated APC-like porcine $\gamma\delta$ T cells (`unconventional_T_cells` C3, C6).
19. **Mair et al. (2014, 2016) & Gerner et al. (2009, 2022)** *Front Immunol* (`PMID: 27471504`, `PMC4943943`), *Dev Comp Immunol* (`PMID: 18601948`), *Annu Rev Anim Biosci* (`PMID: 34724393`) — Characterize porcine `CD3⁻ NKp46(NCR1)^{high}` tissue NK cells (`NK_cells` C1), `CD3⁻ NKp46^{low} CD8αα⁺ CD16^{high}` cytolytic NK cells (`NK_cells` C2), and `CD3⁺ NKp46⁺ CD8⁺` NK-like T cells enriched in pig lungs (`CD8_T_cells` C3, `NK_cells` C3–C4).

### C. Pulmonary Endothelial, Epithelial, Stromal & Granulocyte Lineages
20. **Gillich et al. (2020) & Schupp et al. (2021)** *Nature* (`PMID: 33057196`) & *Circulation* (`PMID: 34097468`) — Define the mammalian alveolar capillary dichotomy between `GPIHBP1⁺ APLNR⁺ HPGD⁺ KIT⁺` General Capillaries (`gCap` / CAP1; `Endothelial` C3) and `APLN⁺ EDNRB⁺ PRX⁺ CA4⁺ KITLG⁺ VWF⁻ TEK⁻` Alveolar Aerocytes (`aCap` / CAP2; `Endothelial` C4), plus `ACKR1⁺` Venous (`C2`), `GJA5⁺ BMX⁺ HEY2⁺` Arterial (`C5`), and `PROX1⁺ CCL21⁺` Lymphatic (`C1`) ECs.
21. **Travaglini et al. (2020) & Sikkema et al. (2023)** *Nature* (`PMID: 33208946`) & *Nat Med* (`PMID: 37291214`, Human Lung Cell Atlas) — Reference atlases for pulmonary epithelial (`AT1`, `AT2`, `Club`, `Multiciliated`, `Basal`, `FOXI1⁺ CFTR⁺` `Ionocytes`, `POU2F3⁺` `Tuft`), stromal (`COL13A1⁺ ITGA8⁺ NPNT⁺` Alveolar Fibroblasts, `PI16⁺ MFAP5⁺` Adventitial Fibroblasts, `WIF1⁺ GREM1⁺ MYH11⁺` Ductal Myofibroblasts/ASM, `NOTCH3⁺ PDGFRB⁺` Pericytes), and immune compartments.
22. **Montoro et al. (2018) & Perniss et al. (2023)** *Nature* (`PMID: 30069044`) & *Eur Respir J* (`PMID: 37531421`) — Define airway `POU2F3⁺` `Tuft-1` (`CHAT^{hi} DCLK2⁺` cholinergic/neuronal; `Epithelial` C9) vs. `Tuft-2` (`SUCNR1⁺ ALOX5AP^{hi} IL17RB^{hi} GFI1B⁺ PTPRC⁺` succinate-sensing immunological; `Epithelial` C8) brush cells and `FOXI1⁺` pulmonary ionocytes (`Epithelial` C10).
23. **Lee et al. (2017) & Narvaez Del Pilar et al. (2022)** *Cell* (`PMID: 28886383`) & *Development* (`PMID: 35302583`) — Define the three-axis lung mesenchymal organization and `LGR5⁺ FGF10^{hi} HGF^{hi} SFRP1^{hi}` epithelial-niche mesenchymal cells (`Stromal` C5).
