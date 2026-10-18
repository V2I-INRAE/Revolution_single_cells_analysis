# Pragmatic cluster regrouping proposal — `job-44542861`

**Status: documentation and analysis proposal only. No barcode cleanup, metadata reassignment, graph refit, clustering, or marker recomputation has been executed by this review.**

Companions:
- [Master cleanup/reassignment ledger](NOTES_BARCODE_REASSIGNMENT_AND_CLEANUP.md): extraction rules, original run identities, and planned transfers.
- [Reviewed nested marker panels](NOTES_REVIEWED_MARKER_PANELS.md): top-5/10/25 candidate panels and specificity limitations.
- [Parent annotation record](../../../../data/job-44542861/annotated_clustered.INFO): recorded cell counts and original annotation validation.

## 1. Decision framework

The aim is fewer defensible **primary identities**, not fewer clusters at any cost. Use three complementary layers:

| Layer | Suggested metadata | Purpose |
|---|---|---|
| Provenance | `orig_subset`, `orig_subcluster_job`, `orig_resolution`, `orig_subcluster_id` | Preserve the exact source assignment, including cells transferred between compartments. |
| Curated identity | `curated_compartment`, `curated_identity_broad`, `curated_identity_fine`, `annotation_confidence` | Broad reporting groups plus defensible secondary distinctions; allow an unresolved value. |
| Programme/state | `state_annotation` and continuous multigene scores | Activation, interferon response, lipid remodelling, ER remodelling, cell cycle, and other programmes without assuming immutable subtypes. |
| QC/decisions | `qc_status`, `reassignment_status` | Separate removal, optional transfer, pending lineage placement, and retained cells. |

A broad roll-up does **not** imply that constituent cells are biologically identical. Preserve original fine labels even where a new curated fine label merges them. Report secondary states alongside broad abundance results when pooling could hide an EVLP-associated change. Programme scores must use available, appropriately annotated features, not a single marker or absence of a sparse transcript.

### Evidence reviewed and limits

- All 12 subset identity reports and the cleanup ledger were read; all 12 `clustree.png` files were inspected sequentially.
- The preceding marker review included significant/prefiltered CSVs and 22 heatmap/dotplot PDFs for retained compartments. Additional unplotted marker candidates were assessed from CSVs, not new plots.
- Cluster and sample count CSVs agree within all 12 subsets. Donor/sample summaries were checked for selected uncertain B, stromal, NK, and unconventional-T groups.
- Clustree records transitions across resolution on one original graph. It is **not** a resampling test, independent donor validation, developmental trajectory, or proof of subtype equivalence. Labels are local to each resolution.
- Parent/child relationships below are visual observations from PNGs. Exact barcode overlap and one-to-one correspondence must be checked on the saved objects before adopting a lower-resolution mapping.
- Marker `pct.2` pools the original other clusters, including subsequently excluded cells. Missing competitor DE rows were not interpreted as zero expression.
- No cleaned pairwise DE, held-out donor classifier, new graph, or new literature search was run for this proposal. Existing literature qualifications are in the marker-panel companion.

## 2. Cleanup and reassignment prerequisites

### 2.1 Reconciled cell accounting

The 12 `cluster_summary.csv` files sum to **238,627** cells. This matches the recorded actual count in `data/job-44542861/annotated_clustered.INFO`. The earlier ledger totals of 228,627 starting / 224,183 retained cells were arithmetic/documentation errors.

| Quantity | Expected cells | Interpretation |
|---|---:|---|
| Original parent / sum of 12 subsets | 238,627 | Recorded original count, not newly filtered here. |
| Designated doublet removals | 4,444 | 1.86% of original cells. |
| Retained after designated removals | **234,183** | 98.14%; expected, not an executed cleaned-object count. |
| Transfers with specified recipient compartments | 16,204 | CD4-to-CD8 and four unresolved-T cluster transfers. |
| Additional cells leaving NK with recipient pending | 4,953 | NK C3+C4; T-lineage placement needs resolution. |
| Total planned mandatory outgoing transfers | 21,157 | Transfers conserve cells; they are not additional removals. |
| Optional mast / ILC2-like compartment transfers | 1,002 / 399 | Conserve cells and keep distinct identities. |

### 2.2 Designated removals

| Source | Selected resolution | Cluster(s) | Cells removed |
|---|---:|---|---:|
| Endothelial | 0.3 | C6 | 318 |
| Epithelial | 0.4 | C7 | 603 |
| Macrophage_Monocytes | 0.4 | C9+C10 | 445+245=690 |
| Stromal | 0.4 | C6 | 319 |
| T_cells_unresolved | 0.4 | C2 | 2,514 |
| **Total** | | | **4,444** |

Use original job/resolution/barcode keys, never a bare cluster number. The epithelial report states that res0.4 C7 and res0.3 C6 have identical barcodes; verify that equivalence before using the alternate rule. Preserve the existing removal decision for macrophage C9 without treating its historical “patrolling nonclassical monocyte” label as a validated porcine identity.

### 2.3 Mandatory transfers and unresolved routing

| Source, at res0.4 | Cells | Recipient / decision |
|---|---:|---|
| CD4 C3+C4 | 3,748+3,035=6,783 | CD8 compartment; fine identity assigned after unified review/refit. |
| Unresolved-T C1+C3 | 4,017+2,213=6,230 | CD8 compartment; preserve source identity as provenance. |
| Unresolved-T C4 | 1,669 | CD4 compartment. |
| Unresolved-T C5 | 1,522 | Unconventional-T / gamma-delta-like compartment, with lineage qualification. |
| NK C3+C4 | 2,892+2,061=4,953 | Remove from NK assignment; final CD8 versus unconventional-T placement remains pending. Do not declare NKT or gamma-delta identity from cytotoxic/NK receptors alone. |

After these actions, the expected CD4 pool is **20,774**, NK pool **19,389**, and original unresolved-T container is dissolved. CD8 contains **33,565** cells before the NK-derived T-cell decision and **38,518** if all 4,953 enter CD8. Unconventional-T contains **22,374**, or **27,327** if all 4,953 enter there; these figures include the 399 ILC2-like cells. A temporary explicit T-lineage review bucket is preferable to a forced recipient label. Do not count it in both recipients.

Mast C4 may stay under Granulocytes or become a separate compartment. ILC2-like C8 may stay in the broad unconventional container or become separate. Neither should be discarded or merged into a T/neutrophil identity because of container membership. Exporting them changes container totals, not the global retained count.

## 3. Sequential subset assessment

All cluster IDs below use each original selected resolution. Counts for proposed merges of original clusters are pre-inflow counts unless stated otherwise. These are annotation proposals, not new clustering results.

### 3.1 B cells — merge the mature-plasma split; preserve ASC states

**Cleanup:** none; 26,790 cells.

| Primary grouping | Original clusters | Cells | Secondary interpretation |
|---|---|---:|---|
| Atypical-memory/MZ-like B | C1 | 6,384 | Provisional identity; retain original label. |
| Activated B | C2 | 5,631 | Activation programme; do not force a pre-GC anatomical interpretation. |
| Costimulatory/B1-like B | C5 | 2,552 | B1-like interpretation provisional. |
| Resting/naive-like B | C6 | 1,145 | Retain pending donor-aware validation. |
| GC B | C9 | 255 | Distinct AICDA/RGS13 programme; do not erase a coherent rare group. |
| **ASC umbrella** | **C3+C4+C7+C8** | **10,823** | Mature core C3+C4=9,708; ER-remodelling C7=625; plasmablast-like C8=490. |

**Strongest immediate merge: C3+C4.** Their LGALS1 detection is 68.6/69.3%, RGS1 98.6/95.6%, and CD69 94.1/92.8%. IGLL1 (review) is a C3 lead but also detected in C8; F8 is a sparse, noncanonical C4 lead. The clustree shows res0.3 C1 dividing into selected C3/C4. A 0.4–0.5 plateau does not validate two mature-plasma cell types.

C7 and C8 branch earlier and have useful programmes; the ASC umbrella must not replace their secondary reporting. All three proposed secondary ASC groups occur across the five donor prefixes; this alone does not prove stable marker differences. C8 is only **plasmablast-like**, not definitively a cycling population: BIRC5 detection is 31%.

**Acceptance/veto checks:** pairwise mature-ASC comparison excluding immunoglobulin-only/clonal and technical differences; donor-replicated C7 ER programme; C8 programme versus both C7 and mature ASCs, with and without cell-cycle features. Treat cell-cycle masking as sensitivity analysis, not a requirement to remove meaningful biology. If an apparent discrete boundary vanishes beyond a continuous programme, keep a score rather than invent a subtype. Do not rely on ERN1 alone: C7/C8 detection is 74.6/75.9%.

**Parameters:** no new graph needed for the C3+C4 metadata merge. Compare existing 0.3 versus 0.4 assignments first. Keep k=100 as the baseline; k=50 is optional sensitivity, not a prerequisite.

### 3.2 CD4 T cells — regroup only after outgoing and incoming transfers

**Cleanup:** transfer C3/C4 out (6,783); add unresolved-T C4 (1,669); expected pool 20,774.

- Candidate broad helper-memory/activated umbrella: **C1+C5=12,776 original cells**. Keep C5 as activated/cycling-enriched secondary state, not a confirmed Treg/Tfh class.
- Retain C2 naive/central-memory and C6 unresolved/quiescence-like provisionally. Do not place C6 in C2 merely because both were called resting; BACH2 does not discriminate C6 from C2.
- Do not assign incoming helper cells a fine cluster solely from their original report; establish their placement in the recipient pool.

**Clustree:** res0.3 C1 divides into selected C1/C5; C2 and C6 have separate earlier branches. Later mixing also cautions against immutable activation subtypes. This supports an optional broad C1+C5 reporting umbrella, not proof that C5 lacks specialised helper populations.

**Parameters:** refit the changed recipient pool before graph-dependent analysis, initially using reported 2,000 HVGs/dims1:20 and k=100. Compare k=50 and res0.2–0.4 without promising that smaller k separates Tregs/Tfh. Require cell-level multigene evidence for any such annotation. Porcine CD8A expression alone does not invalidate memory-CD4 identity.

### 3.3 CD8 T cells — defer hard merging until the enlarged pool is refit

**Cleanup:** add 13,013 cells from CD4/unresolved-T, giving 33,565; final count depends on the NK-derived T-cell decision.

- C1+C4 may share a broad effector/resident-memory reporting umbrella after T-lineage confirmation; preserve GZMK-rich versus integrin-rich secondary programmes.
- Retain C2 cytotoxic effector and C5 type-17-like programmes. MAIT identity remains unconfirmed without appropriate receptor/orthogonal evidence.
- Keep C3 innate-like/lineage-provisional and C6 unresolved state. **Do not merge C3+C4 solely because of their shared graph parent.**
- Incoming CD4 C3 may supply a clearer naive population missing from the original CD8 grouping. Do not force the recipient data into the original six labels.

**Clustree:** selected C3/C4 share res0.3 C3; selected C2 remains a separate branch from res0.1. Selected C1 contributes cells to the C4-derived branch at a later resolution. Topology does not resolve C3's lineage.

**Parameters:** unified clean-pool feature/PCA/integration refit; baseline k=100, sensitivity k=50, dims1:20, res0.2–0.4. Evaluate marker programmes and source mixing without requiring source mixing when it would erase real biological differences. Do not claim the resulting graph establishes a developmental trajectory.

### 3.4 DC — retain four operational groups, without forced subdivision

**Cleanup:** none; 7,436 cells.

Retain C1 pDC, C2 cDC2/inflammatory-DC-like unresolved group, C3 cDC1, and C4 mature/migratory DC. C4 is a state-defined group, but its precursor lineage is not reliably assigned by this review. An optional very coarse cDC umbrella can include C2+C3+C4, without erasing their labels.

**Clustree:** res0.2 C2 divides into selected C2/C4. Do not infer that all mature C4 cells are biologically cDC2-derived. Higher resolutions mainly subdivide pDC rather than demonstrate a clean C2 decomposition.

**Parameters:** k=100/res0.3 is an adequate baseline for the four operational groups. Only pursue targeted C2 k=30/50 tests if multigene, cell-level heterogeneity and the scientific question justify them. The previous claim that lower k will cleanly separate cDC2 from DC3/MoDC is untested.

### 3.5 Endothelial — keep the five identities

**Cleanup:** remove C6=318; expected 13,836.

Keep C1 lymphatic, C2 venous, C3 general-capillary, C4 aerocyte, and C5 arterial. Use inflammatory/adhesion scores across identities instead of creating obligatory activated-EC subtypes. SELE/SELP are shared between C2/C5 and cannot establish a venous identity or its subdivision alone.

**Clustree:** the five identities are recognisable at lower resolution, but selected gCap C3 and removed C6 share a res0.2 parent. Higher-resolution subdivisions are chiefly lymphatic. Selecting lower resolution does not clean the object.

**Parameters:** no graph rerun required for retained metadata labels. If new neighbourhood/trajectory analyses are needed, refit after removal, retaining k=100/dims1:15 baseline and optionally comparing k=50. No required resting/activated/systemic-venous split.

### 3.6 Epithelial — seven broad groups, with preserved programmes

**Cleanup:** remove selected C7=603; expected 9,426.

| Broad grouping | Original clusters | Cells | Preserve as secondary programmes |
|---|---|---:|---|
| Ciliated | C1+C5 | 2,958 | Apical-tip/radial-spoke-rich versus dynein-rich. |
| Secretory | C2+C3 | 3,201 | ATP13A5/distal-like versus canonical secretory intensity. |
| Tuft/brush-like, provisional | C8+C9 | 768 | Immune-like versus neuronal/cholinergic-like. |
| AT2 | C4 | 1,287 | Keep distinct. |
| AT1 | C6 | 662 | Keep distinct. |
| Ionocyte | C10 | 292 | Keep distinct. |
| Basal | C11 | 258 | Keep distinct. |

**Clustree:** res0.3 C1 gives selected C2/C3; res0.3 C2 gives C1/C5; res0.3 C4 gives C8/C9. This offers a practical coarse-resolution starting point, subject to exact barcode cross-tabulation.

**Strongest merge:** C2+C3; SCGB1A1 detection is 88.3/86.6% and TFF3 95.9/95.5%, so C3 has no defensible five-gene positive discriminator panel. Distal/proximal anatomy remains provisional. C1/C5 share ciliated identity, with a quantitative programme split. C8/C9 pooling is an optional reporting umbrella, not proof of identical states or mature tuft identities. C8 requires coherent cell-level lineage/QC review; PTPRC alone neither proves a doublet nor validates the tuft label.

**Parameters:** existing res0.3 is the first comparison, after original-barcode exclusion. No k change is needed to obtain broader annotations. Preserve rare basal/ionocyte groups; do not tune parameters to a preselected group count.

### 3.7 Granulocytes — one broad neutrophil identity, several states

**Cleanup:** no designated doublets. Mast C4=1,002 is optional compartment export, not removal.

| Broad identity | Original clusters | Cells |
|---|---|---:|
| Neutrophils | C1+C2+C5 | 17,950 |
| Eosinophils | C3 | 3,814 |
| Basophils | C6 | 271 |
| Mast cells | C4 | 1,002 |

Preserve C1 inflammatory, C2 interferon-high, and C5 granule-rich/maturation-like programmes. The minimal merge is **C1+C2=17,232**; C5 remains worth secondary analysis and should not be called definitively immature from RNA alone.

**Clustree:** selected C1/C2 share res0.2 C1. C5 already has a distinct branch at res0.1, as do eosinophils/mast/basophils; higher resolution mainly fragments inflammatory neutrophils.

**Parameters:** retain the reported 500-HVG/dims1:10/k=100 baseline. No new graph is necessary for broad reporting. Do not lose basophils because they are rare, or merge them with mast cells based on shared FCER1A/MS4A2.

### 3.8 Macrophage/monocytes — four broad groups, state-aware analyses

**Cleanup:** remove C9=445 and C10=245; expected 41,524. Keep C7 with its explicit endothelial-signal caveat.

| Broad grouping | Original clusters | Cells | Preserve |
|---|---|---:|---|
| PAM-like | C1+C2+C3 | 27,735 | C2 lipid-remodelling; original C1/C3 provenance, without asserting homeostatic versus quiescent subtype identity. |
| Interstitial-like macrophages, niche provisional | C4 | 5,434 | CCL8/FCGR2B/MAF-associated programme. |
| Monocyte/MoMac continuum | C5+C6+C8 | 7,217 | C5 repair-associated, C6 inflammatory, C8 differentiating-like programmes. |
| MARCO/TIMD4 macrophages, niche provisional | C7 | 1,138 | Distinct programme and EC-signal QC flag. |

These are coarse reporting groups. Classical-monocyte versus MoMac comparisons should continue to stratify the C5/C6/C8 provenance/programmes rather than use only the pooled total.

**Clustree:** C1/C2 share res0.2 C1, whereas C3 is already separate at 0.2; all three share a res0.1 ancestor. C5/C6 share res0.3 C2, and C8 separates earlier. C7 shares res0.3 C6 with removed C10. Thus no single lower resolution implements all desired biological merges and exclusions.

C3's weak positive signature is not a deletion criterion. Do not merge C7 into C4 by cross-species niche analogy, and do not use PLVAP/CDH5 as macrophage identity markers. NREP is a dataset-level C1 candidate, not an established PAM-lineage marker; neither its presence nor lack of retrieved supporting literature alone establishes contamination.

**Parameters:** annotation grouping first. Keep k=100 as baseline; k=50 is an optional post-cleanup sensitivity, not a demonstrated optimum.

### 3.9 NK — resolve transfers, then avoid forcing the low-information cluster

**Cleanup:** transfer C3+C4=4,953 out of NK; expected remaining 19,389.

Retain C1 cytokine-responsive-like, C2 cytolytic, C6 quiescent NK/ILC-like, and C7 TCF7-high NK/ILC-like programmes provisionally. C5 should be labelled low-information/unclassified state rather than a new subtype. BTC/ANXA1/TXNRD1 are sparse leads, not a validated classifier.

**Clustree:** C1/C7 share res0.2 C1; C6 has a separate earlier branch. C5 receives contributions from both major NK branches and shows later mixing. Consequently, do not automatically map C5 to C2. An optional C1+C7 reporting umbrella contains 9,785 original cells, conditional on compatible lineage evidence; it does not demonstrate anatomical ILC1 or progenitor function.

**Parameters:** refit the clean pool for new graph-dependent analysis; k=100 versus 50, dims1:15, res0.2–0.4. Review QC, donor effects, T-lineage evidence, and neighbour assignment with an explicit rejection/unresolved option. Do not remove C5 simply because it lacks a five-gene panel.

### 3.10 Stromal — default five coarse groups; optional four-group roll-up

**Cleanup:** remove C6=319; expected 11,079.

| Default broad grouping | Original clusters | Cells | Qualification |
|---|---|---:|---|
| Matrix/adventitial-like fibroblast umbrella | C1+C4 | 5,994 | Preserve both original clusters and C1-associated activation score; not a proven single subtype or anatomical niche. |
| Alveolar-like fibroblasts | C3 | 1,617 | Retain NPNT/BMPER/TF programme. |
| Niche-like fibroblasts | C5 | 529 | Retain LGR5/CRABP1/SFRP1/HGF/FST programme. |
| WIF1/GREM1 contractile group | C2 | 2,751 | Myofibroblast/airway-contractile composition unresolved. |
| Vascular mural group | C7 | 188 | Pericyte/VSMC composition unresolved; retain rare coherent group. |

**Optional roll-up:** C3+C5=2,146 as a broad alveolar/niche axis, yielding four groups. Keep C3/C5 separate by default for niche/signalling analysis. Never merge C2/C7 simply because both express ACTA2/MYH11.

**Clustree:** selected C1/C4 share res0.2 C1 and split at 0.3. C3/C5 share res0.3 C3 and split at 0.4. Contractile C2 and mural C7 already have distinct res0.1 branches. Split order is not an independent biological ranking.

**Important correction:** C1 activation does not establish C1/C4 subtype equivalence. Both occur in all five baseline controls:

| Baseline sample | C1 | C4 |
|---|---:|---:|
| REVO26-C | 144 | 69 |
| REVO27-C | 47 | 16 |
| REVO29-C | 71 | 23 |
| REVO30-C | 182 | 95 |
| REVO31-C | 68 | 33 |

Their coexistence rules out the simple story that C1 exists only as an EVLP-induced replacement for C4; it does not by itself prove distinct immutable subtypes. Shared COL14A1 detection (97.1/99.2%) and LUM (93.0/94.1%) support a coarse matrix umbrella, while C4's ADAMTS16/DKK2/ACKR2/MFAP5 programme warrants preservation. C1 IL6/PTGS2/IGF1 detection is sparse; zero counts do not mean absence of activation.

**Acceptance/veto checks:** compare C1/C4 lineage programmes across donors and along a multigene activation score. If C4 markers persist after accounting for state/QC, retain a fine identity distinction. Assess C5's programme against C3 across donors, especially before pooling niche/signalling results. No physical localisation or causal C4-to-C1 conversion has been demonstrated.

**Parameters:** no graph needed for coarse metadata roll-up. If graphs are refit after cleanup, k=100 baseline versus 50 with dims1:15 is a bounded first comparison. Targeted k≈30 experiments on the contractile/mural questions require cell-level evidence; they cannot guarantee a pericyte/VSMC or myofibroblast/ASM split. C7's size relative to k is a sensitivity concern, not proof that two lineages were merged.

### 3.11 Unconventional T — broad lineage plus activation states

**Cleanup:** add unresolved-T C5=1,522; keep ILC2-like C8 distinct or export it without loss. NK-derived T placement remains separate and unresolved.

- Broad CD2-low gamma-delta-like umbrella: **C1–C6=19,998 original cells**, provisionally **21,520** with the 1,522 incoming cells.
- Keep C7 CD2-high/EOMES-high gamma-delta-like group=455 separate.
- Keep C8 ILC2-like=399 separate; it is not a gamma-delta activation state.
- Preserve C1 resting/circulating-like, C2 stress/early-activation-like, C3 effector-like, C4 poised, and C6 costimulatory programmes. C5 is low-information/unclassified state, not a validated subtype.

**Clustree:** C2/C5 share res0.3 C1; C1 joins them at a res0.2 parent; C6 shares an earlier ancestor. C3/C7 also share an earlier parent despite different CD2/EOMES programmes, while C4 separates early. A common parent is not a reason to merge C7 with the CD2-low group or with incoming NK-like T cells.

C1 is 55.4% control, C2 only 1.5%, and C6 4.1%, from the sample-count tables. These are cluster compositions, not donor-normalised enrichment estimates or evidence of a temporal C1-to-C2/C6 transition. A broad lineage count must be accompanied by state analysis or it could hide the relevant EVLP response. Gamma-delta assignments remain qualified where direct TCR evidence is lacking; CD247 alone is not T-lineage proof.

**Parameters:** refit after incoming-cell integration; k=100 versus 50, dims1:20, res0.2–0.4. Retain ILC2-like cells separately for any gamma-delta-specific graph. No obligation to reproduce six CD2-low states or to infer a trajectory.

### 3.12 T cells unresolved — dissolve, do not recluster

Remove C2=2,514 and transfer C1/C3 to CD8, C4 to CD4, C5 to unconventional-T. No cells remain in the original container.

**Clustree counterexample:** selected C1 and C4 share a res0.2 parent despite their differing recipient-lineage interpretations. The later 0.3–0.6 plateau therefore cannot establish that the unresolved container is one biologically valid cell type.

No new k/resolution search is warranted here. Preserve provenance and review transferred cells in the recipient pools; do not paste together saved graphs or equate source-cluster numbers across objects.

## 4. Bounded parameter and validation plan

### What to do before any new graph

1. Extract original barcodes using `(subset, subcluster job, selected resolution, cluster)` keys.
2. Assert unique source ownership, disjoint removal/retention sets, and exactly one curated destination per retained barcode. A pending lineage bucket is a valid explicit destination.
3. Apply designated removals first; route transfers next; preserve mast/ILC2-like identities and optional handling.
4. Reconcile expected global retention **234,183** and scenario-specific recipient counts. New barcode-level results take precedence if a discrepancy is found; do not silently force totals.
5. Create broad/fine/state annotations without overwriting source labels. Marker-poor cells may remain unresolved rather than be deleted or confidently assigned by one marker.

### Graph policy

| Situation | Recommended first action |
|---|---|
| Coarse reporting in unchanged cell pools | Metadata roll-up; existing graph need not change. |
| CD4/CD8/unconventional pools after incoming transfers | Refit subset features, scaling/PCA, integration, neighbours and clustering before graph-dependent conclusions. Reuse appropriate normalised RNA under the existing pipeline; do not concatenate old embeddings/graphs. |
| NK after outgoing transfers | Refit before deciding ambiguous boundaries; old labels remain useful provenance. |
| Other compartments after removals | Metadata retained labels may be used without a new graph; refit if using cleaned neighbourhoods/trajectories or testing boundary stability. |
| First graph sensitivity test | **k=100 baseline versus k=50**, reported features/dimensions held fixed, resolutions 0.2/0.3/0.4. Start with the recorded seed; check leading choices with another seed if needed. |
| Rare-lineage question supported by cell-level data | Optional targeted k≈30 experiment, not a universal prescription. |

Use the reported original subset dimensions (1:20 for B/CD4/CD8/DC/epithelial/macrophage/unconventional/unresolved; 1:15 endothelial/stromal/NK; 1:10 granulocytes), not today's global defaults as a substitute for per-run provenance. Keep HVG/dimension changes separate from the first k comparison to avoid confounding sensitivity sources. Smaller k may fragment data; larger k may oversmooth. Neither is a direct “number of biological types” dial, and no k optimum has been demonstrated here. UMAP `n.neighbors` is not clustering `k.param`.

Implementation note: the inspected `R/clusters/subcluster.R` / `R/utils/cli.R` accept `--hvg` and `--dims`, **not `--k-param`**. A future per-run k override needs an explicitly approved code change with recorded settings. Do not change the shared configuration globally just to force one subset's outcome. No such change or job submission is part of this document.

### Evidence needed before adopting fine merges/splits

- Inspect **cell-level coexpression** and coherent multigene programmes, rather than cluster averages alone.
- Evaluate the closest retained competitor with cleaned pairwise DE and coverage, accounting for donor, QC and relevant condition/state effects. Treat donors/biological experimental units—not cells—as independent replication; timepoints from the same donor are repeated observations.
- Check both baseline and EVLP samples. Coexistence at baseline limits a simple activation-only explanation but does not establish fixed subtypes. Condition enrichment can be meaningful state biology and is not automatically a batch artifact.
- For important boundaries, test donor-held-out programme/classifier generalisation with a rejection option and control for class imbalance. There are five donor prefixes, not 23 independent donors. Do not overpromise validation for a 188-cell mural group or other rare populations.
- Recalculate marker panels against curated labels. Old pooled-rest panels are candidates, not validated post-regrouping classifiers.
- Compare broad and secondary-state composition results so that pooling does not conceal a response. Use continuous state scores as complementary—not replacement—analyses.
- Do not regress away activation, cell cycle or condition by default. Use justified sensitivity analyses; the original pipeline already includes cell-cycle handling/regression, so “cycling” labels need particular care.

## 5. Constrained scientific challenge and resulting decisions

Exactly **one** sci-oracle consultation was used for this regrouping proposal, restricted to the B-cell ASC and stromal merger uncertainties. It was instructed to reason only from supplied local evidence, with no tools, file inspection, literature search, delegation, or modification. No further oracle consultation was needed for the other subsets.

| Question | Reviewer challenge | Adopted interpretation |
|---|---|---|
| Are several plasma/ASC clusters necessary as primary identities? | Endorsed C3+C4 mature merge and broad four-cluster ASC umbrella only with secondary reporting. | Mature C3+C4 merge is the strongest recommendation; C7/C8 remain available for state-specific analyses. |
| Does sparse BIRC5 prove a discrete cycling C8 type? | No; require a reproducible programme beyond one cell-cycle marker. | “Plasmablast-like”; cell-cycle-feature sensitivity rather than mandatory regression or a definitive cycling identity. |
| Does C1 activation prove stromal C1/C4 are one subtype? | No; broad roll-up is different from subtype equivalence. | Default matrix/adventitial-like umbrella; preserve C1/C4 fine provenance and programme differences. |
| Should stromal C3/C5 be merged by default? | Specialized C5 niche factors warrant caution. | Optional alveolar/niche roll-up; keep C3/C5 separate by default in secondary/niche analyses. |

The reviewer suggested activation-marker-negative comparisons; these should be implemented as activation-score-aware analyses, not literal zero-count gating, because dropout precludes equating non-detection with inactivity. Its recommendation to retain C7/C8 secondary analyses does not prohibit valid broad ASC abundance reporting; the two reporting levels address different questions.

## 6. Source index and decision status

All paths below are relative to this folder. Each marker directory contains the identity report, `cluster_summary.csv`, `sample_cluster_counts.csv`, marker CSVs, and the previously reviewed heatmap/dotplot PDFs.

| Sequential order | Subset / job | Selected resolution | Identity report | Inspected clustree |
|---|---|---:|---|---|
| 1 | B_cells / job-44614229 | 0.4 | [Report](B_cells/job-44614229/markers_res0.4/cluster_identities_res0.4.md) | [PNG](B_cells/job-44614229/clustree.png) |
| 2 | CD4_T_cells / job-44614236 | 0.4 | [Report](CD4_T_cells/job-44614236/markers_res0.4/cluster_identities_res0.4.md) | [PNG](CD4_T_cells/job-44614236/clustree.png) |
| 3 | CD8_T_cells / job-44614237 | 0.4 | [Report](CD8_T_cells/job-44614237/markers_res0.4/cluster_identities_res0.4.md) | [PNG](CD8_T_cells/job-44614237/clustree.png) |
| 4 | DC / job-44614230 | 0.3 | [Report](DC/job-44614230/markers_res0.3/cluster_identities_res0.3.md) | [PNG](DC/job-44614230/clustree.png) |
| 5 | Endothelial / job-44614231 | 0.3 | [Report](Endothelial/job-44614231/markers_res0.3/cluster_identities_res0.3.md) | [PNG](Endothelial/job-44614231/clustree.png) |
| 6 | Epithelial / job-44614232 | 0.4 | [Report](Epithelial/job-44614232/markers_res0.4/cluster_identities_res0.4.md) | [PNG](Epithelial/job-44614232/clustree.png) |
| 7 | Granulocytes / job-44614233 | 0.3 | [Report](Granulocytes/job-44614233/markers_res0.3/cluster_identities_res0.3.md) | [PNG](Granulocytes/job-44614233/clustree.png) |
| 8 | Macrophage_Monocytes / job-44614234 | 0.4 | [Report](Macrophage_Monocytes/job-44614234/markers_res0.4/cluster_identities_res0.4.md) | [PNG](Macrophage_Monocytes/job-44614234/clustree.png) |
| 9 | NK_cells / job-44614238 | 0.4 | [Report](NK_cells/job-44614238/markers_res0.4/cluster_identities_res0.4.md) | [PNG](NK_cells/job-44614238/clustree.png) |
| 10 | Stromal / job-44614235 | 0.4 | [Report](Stromal/job-44614235/markers_res0.4/cluster_identities_res0.4.md) | [PNG](Stromal/job-44614235/clustree.png) |
| 11 | unconventional_T_cells / job-44614239 | 0.4 | [Report](unconventional_T_cells/job-44614239/markers_res0.4/cluster_identities_res0.4.md) | [PNG](unconventional_T_cells/job-44614239/clustree.png) |
| 12 | T_cells_unresolved / job-44614240 | 0.4 | [Report](T_cells_unresolved/job-44614240/markers_res0.4/cluster_identities_res0.4.md) | [PNG](T_cells_unresolved/job-44614240/clustree.png) |

**Recommended first annotation changes, after cleanup:** B C3+C4; epithelial C2+C3 and C1+C5; granulocyte C1+C2 as a broad state-aware group; stromal C1+C4 as a coarse umbrella. Broader ASC/PAM/neutrophil/gamma-delta roll-ups are reporting hierarchies, not mandates to erase fine biology. All new regroupings remain proposed until explicitly implemented and validated; original identity reports are unchanged.
