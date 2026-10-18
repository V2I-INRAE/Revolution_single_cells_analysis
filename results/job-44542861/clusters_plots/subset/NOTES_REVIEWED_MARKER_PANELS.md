# Reviewed nested marker panels — job-44542861

## Scope and interpretation

Saved from the completed conversational review at the owner's request. Companion to [the cleanup ledger](NOTES_BARCODE_REASSIGNMENT_AND_CLEANUP.md). Existing identity reports and ledger are unchanged.

All 12 identity reports were reviewed, with marker CSVs and 22 heatmap/dotplot PDFs for the retained compartments. Three literature scouts checked biological interpretation at abstract level. No new literature search was conducted to save these tables.

**Nested columns:** top 5; add five for top 10; add fifteen for top 25. **These are exploratory within-original-subset panels, not validated whole-lung or post-cleanup classifiers.** Expanded panels contain supporting state/identity genes, not 25 individually exclusive discriminators. Withhold means insufficient discrimination, not missing work. Annotation warnings (†) preserve the literal source labels `CD1.1 (review)` and `IGLL1 (review)`.

Selection prioritised detection coverage, enrichment, nearest-cluster contrasts where available, and coherent expression patterns. Every proposed gene has a positive significant marker row and positive pooled detection difference; this does not establish nearest-neighbour specificity. Missing competitor DE rows were not treated as zero expression. `pct.2` pools the original other clusters, including later exclusions.

Heatmaps use scaled RNA expression and at most 500 sampled cells per cluster. Dotplot size is detection and colour is scaled average expression. Unplotted additions were assessed from CSVs, not visually validated in new plots.

### Excluded/deferred clusters

- Macrophage/monocyte C9/C10; endothelial C6; epithelial C7; stromal C6.
- CD4 C3/C4; NK C3/C4; all T_cells_unresolved.
- Optional reassignment candidates conservatively deferred: granulocyte C4 (mast) and unconventional-T C8 (ILC2).
- All B and DC clusters retained; CD8 identities remain provisional pending integration/reclustering.

## Macrophages and monocytes

| Cluster / working identity | Top 5 | Add 5 → 10 | Add 15 → 25 |
|---|---|---|---|
| C1 — NREP-high PAM-like; state label provisional | NREP, CYP7A1, GPR82, RFLNA, TMEM236 | NAAA, PPBP, SLC2A12, XDH, FCAR | BHLHE41, FCGR1A, ADGRE3, UPP1, PDLIM1, PIK3IP1, WDR93, DOK6, PPM1J, P2RY14, CLEC7A, TNFSF13, SLC2A9, SORT1, PADI2 |
| C2 — Lipid-remodelling PAM | OCSTAMP, TUBA8, S100A2, FABP5, LPL | CD24, GPD1, SLC45A3, SPTBN1, SPRY2 | ABHD8, FABP3, AVPI1, FJX1, GCNT1, PM20D2, GALNT4, TCN1, LRIG3, HAPLN4, CKB, TASL, SH3BGRL2, TYSND1, ARL11 |
| C3 — basal PAM; weak separation | **Withhold** | **Withhold** | **Withhold** |
| C4 — Interstitial/perivascular-like macrophage; niche provisional | FCGR2B, CCL2, MAF, RNF128, CCL8 | DNASE1L3, IGF1, ABCA9, JAML, HTRA1 | PID1, ACP5, CADM1, GPR183, CCL14, FGL2, ADORA3, MS4A7, STAB1, HAVCR2, GPR31, AGAP1, CCND2, ST3GAL5, C3AR1 |
| C5 — Repair-associated monocyte/early MoMac | PLBD1, EREG, S100A4, CCR2, PYCARD | AREG, IL1R1, CHI3L2, TIMP1, SULT1A3 | CCR7, CD52, ITGB7, COL18A1, MILR1, TNIP3, EBI3, AHRR, SATB1, RGS18, VCAN, F13A1, IL1B, IL1R2, DGAT2 |
| C6 — Inflammatory classical monocyte | RETN, KCNJ15, GPR84, CHI3L1, TREM1 | SDS, IL1R2, STEAP4, NRG1, G0S2 | IL1B, THBS1, TGM1, TGM3, DGAT2, NFE2, SELL, BATF3, NLRP12, OSCAR, MCEMP1, CSF3R, ADAM8, CYP1B1, ACOD1 |
| C7 — MARCO/TIMD4 macrophage; niche provisional | MARCO, IGFBP2, TIMD4, IDO1, C4BPA | CCL19, SLAMF7, PDCD1LG2, TRPC6, DSE | KIAA1217, GPRIN3, DYSF, ITGA9, ACER2, SDC3, GHSR, NEBL, TSPAN33, CD180, ADGRG6, SAA3, CXCL10, C1R, NUPR1 |
| C8 — Differentiating MoMac | LRMDA, PRKG1, CHST11, LRCH1, RAP1GAP2 | ZFYVE9, MAML3, ABHD2, SH3PXD2B, IGF2BP3 | HDAC9, GFOD1, HDAC4, OSBPL3, SIPA1L1, DENND5A, TBXAS1, SCARB1, PCNX1, WWOX, SIK3, TNS3, XYLT1, JAK2, GAB2 |

C3 has no convincing positive signature versus C1/C2. MYO1D detection is 21.3% in C3 versus 52.2% in C1 and 57.4% in C2. C5/C6/C8 share monocyte genes; FCN1/VCAN/F13A1/IL1B alone do not distinguish them. C7's PLVAP/CDH5 signal is a QC caveat, not macrophage identity evidence.

**Literature qualification:** NREP is a dataset-level C1 discriminator, not a verified porcine PAM-lineage marker; NREP alone does not prove contamination either. C1's earlier “homeostatic” interpretation and C7's anatomical niche are provisional. MARCO supports porcine PAM identity as well, so it cannot establish an interstitial niche. ENPP1/CLEC12B support PAM enrichment but not C1–C3 separation.


---

## B cells, DCs and granulocytes

**Top 10 = top 5 + next 5; top 25 = top 10 + next 15.** These are exploratory **within-original-subset** candidate panels, not validated classifiers. Full panels were checked against `markers_significant.csv`: every listed gene has positive detection enrichment (`pct.1 > pct.2`). Missing competitor DE rows were **not** interpreted as zero expression.

**“Withhold” means insufficient evidence for the requested highly discriminating panel—not missing work.** `†` preserves an annotation warning: the CSV literally labels these genes `CD1.1 (review)` or `IGLL1 (review)`.

### B cells

| Cluster / identity | Preferred top 5 | Next 5 → top 10 | Next 15 → top 25 |
|---|---|---|---|
| **C1 atypical-memory/MZ-like** | ZBTB32, FCRL5, CD72, PTPN22, SERPINB9 | ITGAD, HCK, FCGR2B, CD1.1†, TBX21 | SCIMP, CNR2, YBX3, C4BPB, LMO2, SPIB, SLA, PLEK, ITGB5, ST3GAL2, ADCY7, LY86, CLEC17A, GNS, CD302 |
| **C2 activated B — provisional** | SATB1, IRS2, CCR7, DUSP2, MAFF | CD83, NRIP1, ID3, MYBL2, PHACTR1 — **provisional combination** | **Withhold:** insufficient demonstrated discrimination from other activated/GC B cells |
| **C3 plasma A — insufficient separation** | **No supported top 5.** IGLL1† is the strongest lead | **Withhold** | **Withhold** |
| **C4 plasma B — insufficient separation** | **No supported top 5.** F8 is a tentative, noncanonical lead | **Withhold** | **Withhold** |
| **C5 costimulatory/B1-like — provisional** | BHLHE41, TCF7, CD80, FGD2, DUSP4 | TNFRSF1B, SPN, CD86, SEMA7A, FCRLA — **provisional combination** | **Withhold:** extended genes mostly provide shared activation/state support |
| **C6 resting/naive-like** | FHIT, KCNQ5, TMCC3, SMYD3, EBF1 | AFF3, ARHGAP24, BACH2, PREX1, TOX | TTC28, RAD51B, WWOX, BTBD9, AHCYL2, UBE2E2, CAMK1D, PITPNC1, CCND3, ANKRD44, ATXN7L1, MAML2, TBC1D5, AKT3, PRKCE |
| **C7 ER-remodelling ASC** | FNDC3B, ABTB3, GAB2, CADPS2, OSBPL10 | ERN1, SCFD2, ABTB2, PITPNM2, ENTPD1 | CCSER1, KIAA1217, TNS3, PLPP1, MED12L, HDLBP, SEC24A, CREB3L2, EIF2AK3, PTPRE, RBM47, TMCC1, BCAR3, LSS, AFG2A |
| **C8 plasmablast-like/cycling-enriched** | TOM1L1, RFC3, IL6R, PDIA5, LAYN | PADI2, CLDN4, BIRC5, CRELD2, PLAAT3 | ERP29, AQP3, SDF2L1, DUT, RAB37, TPST2, CNPY2, ADAM19, ATOX1, DNAJC3, SLAMF1, NHERF1, CLIC5, SLC38A12, GMPPB |
| **C9 germinal-centre B — strongest B panel** | RGS13, AICDA, NUGGC, S1PR2, GCSAM | BCL7A, MYBL1, MEF2B, P2RY8, BCL6 | AFF2, NUAK1, RRM2B, NEIL1, FANCI, ELL3, OSBPL3, LPP, RAB3IP, DUSP6, SIAH1, MTMR12, NDC1, CCDC181, NREP |

### Dendritic cells

| Cluster / identity | Preferred top 5 | Next 5 → top 10 | Next 15 → top 25 |
|---|---|---|---|
| **C1 pDC** | VPREB1, CD4, SMIM5, BLNK, ATP13A2 | TCF4, NCR1, TDRD1, COBLL1, CMKLR1 | BCL11A, IRF7, TLR7, TLR9, NEXN, IGLL1†, CD8B, CD93, NIBAN3, DPYS, BHLHA15, GAPT, AOAH, WASF3, FGFR1 |
| **C2 cDC2/DC3–MoDC mixture** | CD1.1†, FCER1A, MS4A2, S100A4, CD207 | FCGR3A, MRC1, IL1B, CD163, TLR4 | CXCL2, AREG, CD300C, AXL, IL1A, PTAFR, HBEGF, CFP, CEBPB, FCGR2B, KLF4, SIRPA, CHPT1, AMPD3, IGSF6 |
| **C3 cDC1** | XCR1, CADM1, ANPEP, CLNK, DNASE1L3 | BATF3, LY75, ZNF366, WDFY4, SNX22 | XIRP2, ME1, ATP9A, GNA14, SLAMF8, ARHGAP18, RAB7B, RAB32, SCIMP, ZBED2, NEBL, CLEC7A, ITGAE, TOX, RIMS2 |
| **C4 mature/migratory DC** | IDO1, FSCN1, CCR7, TNFRSF11A, TMEM150C | IL4I1, PDCD1LG2, CD274, S1PR1, CCL22 | SLCO5A1, S1PR5, CD200, ANKRD33B, NXN, GPR157, STAT4, IL21R, KREMEN1, CACNB3, PLEKHG1, ETV7, WNT10B, ACHE, LAMP3 |

**C2’s panel identifies the existing mixed cluster; it does not separate cDC2 from DC3/MoDC.** pDC genes such as `VPREB1/IGLL1/CD4` discriminate inside this DC subset, not against every lung lineage.

### Granulocytes

**C4 mast cells deferred** because the parent requested conservative exclusion of optional reassignment candidates.

| Cluster / identity | Preferred top 5 | Next 5 → top 10 | Next 15 → top 25 |
|---|---|---|---|
| **C1 inflammatory neutrophils — state panel** | CSF1, EDN1, IL1A, ICAM1, CXCL2 | OLR1, NFKB2, CCL5, CCL4, NLRP3 — **provisional combination** | **Withhold:** insufficient demonstrated nearest-cluster specificity |
| **C2 ISG-high neutrophils — state panel** | IFIT2, IFIT3, VNN1, IFIT1, STEAP4 | DHRS9, OAS2, FCGR3A, ISG15, HERC6 — **provisional combination** | **Withhold:** shared interferon/maturation programme |
| **C3 eosinophils** | EPX, CCR3, ZBTB16, FAR1, IGF1R | OR51E1, ALDH3A2, CLNK, ALOX5, PYGL | STAP1, PPA2, GRK3, PAFAH2, MFN2, TMEM30A, ST3GAL1, MIDEAS, PGD, MCTP1, ADGRG3, STXBP5, CEMIP2, HES1, DACH1 |
| **C5 granule-rich/immature-like neutrophils** | LTF, LCN2, MMP8, S100A8, S100A12 | FCN1, S100A9, CD52, S100A11, MEGF9 | NRG1, IL1R2, CD24, UPP1, OSCAR, TNFAIP6, MMP25, FFAR2, TSPO, CDKN2D, FBXO9, EHD1, RAB10, IGSF6, SLC26A6 |
| **C6 basophils** | TMPRSS11E, FCER1A, MS4A2, GCSAML, CA8 | RHEX, ITGA2, HRH4, GATA2, GPR82 | SYTL3, RNF150, PLBD1, HDAC9, MMP9, XK, RUNX2, IKZF2, AP1S3, PADI2, SC5D, GGA2, MAML3, PDE7B, ETS1 |

`TMPRSS11E` is an excellent **dataset-level candidate**, not an established canonical basophil marker.

## Main numerical caveats

| Issue | Direct evidence | Consequence |
|---|---|---|
| **B C3/C4 are poorly separated** | C3/C4 detection: LGALS1 **68.6/69.3%**, RGS1 **98.6/95.6%**, CD69 **94.1/92.8%** | These cannot support the proposed C4-specific signature. Shared JCHAIN/PRDM1/TNFRSF17/XBP1/MZB1 identify plasma cells, not the split |
| **B C3 IGLL1† is not exclusive** | **80.4% C3; 57.8% C8** | Useful lead, insufficient for a five-gene discriminatory panel |
| **B C4 F8 needs caution** | **25.2% C4 versus 5.8% pooled rest**, log2FC **2.27** | Sparse, noncanonical lead; inspect annotation/cell-level pattern rather than promote it to a signature |
| **B C2 shares activation markers** | CD83 **82.0% C2 versus 87.1% C9**; CCR7 **73.8% C2 versus 61.6% C5** | Use combinations; withhold a confident top 25 |
| **B C5 SPN is not C5-specific** | **42.5% C5 versus 61.8% C8** | Supporting gene only, removed from core five |
| **B C7 ERN1 is not C7-specific** | **74.6% C7 versus 75.9% C8** | Supporting UPR marker, not the core discriminator |
| **DC C4 LAMP3 is shared with pDCs** | **94.0% C4 versus 90.6% C1** | Prefer IDO1/FSCN1/TMEM150C; retain LAMP3 only as supporting context |
| **Granulocyte C5 MMP9 is stronger in basophils** | **38.2% C5 versus 73.1% C6** | Excluded from the C5 panel |
| **Sparse lineage markers need coverage partners** | EPX **30.9%**; LTF **30.1%**; MMP8 **28.3%** in target clusters | Pair with ZBTB16/FAR1 or S100A8/S100A12 rather than relying on the sparse marker alone |

### Strong positive anchors

“Rest” means the **original pooled subset**, not cleaned pairwise comparisons.

| Cluster | Marker | Target / rest detection | log2FC |
|---|---|---:|---:|
| B C1 | ZBTB32 | 44.0 / 4.8% | 3.51 |
| B C6 | KCNQ5 | 72.6 / 21.2% | 3.69 |
| B C8 | TOM1L1 | 55.5 / 1.2% | 5.47 |
| B C9 | RGS13 | 55.3 / 0.2% | 8.83 |
| B C9 | AICDA | 63.9 / 0.8% | 6.77 |
| DC C1 | ATP13A2 | 94.5 / 4.5% | 5.34 |
| DC C2 | MS4A2 | 88.2 / 12.3% | 4.26 |
| DC C3 | XCR1 | 79.4 / 12.9% | 4.68 |
| DC C4 | FSCN1 | 83.6 / 4.3% | 5.76 |
| Gran C1 | ICAM1 | 75.2 / 40.5% | 1.85 |
| Gran C2 | IFIT2 | 48.0 / 16.5% | 2.26 |
| Gran C3 | EPX | 30.9 / 0.4% | 6.85 |
| Gran C5 | LTF | 30.1 / 0.9% | 8.15 |
| Gran C6 | TMPRSS11E | 60.5 / 0.5% | 8.83 |

## What visual inspection contributed

- **B C3/C4:** the heatmap shows a large shared plasma-cell block. C4’s existing dotplot top five is dominated by unannotated identifiers, with F8 the remaining named gene. This is not evidence for 25 named discriminators.
- **B C9:** a coherent GC block and restricted RGS13/AICDA/NUGGC/S1PR2 dots provide substantially cleaner discrimination.
- **B C7/C8:** FNDC3B extends into C8; C8’s TOM1L1/LAYN/PADI2/CLDN4 dots are more restricted than shared IGHM.
- **DC:** broad blocks separate well; C2 shows heterogeneous coverage of cDC2 versus inflammatory markers. FCER1A/CD207 need higher-coverage partners.
- **Granulocyte C1/C2:** inflammatory and interferon programmes separate in aggregate but overlap; extended panels should not be advertised as exclusive subtypes.
- **Granulocyte C5:** restricted but sparse LTF/LCN2/MMP8 signals are complemented by broader S100A8/S100A12 detection.

Heatmaps represent scaled RNA expression in **at most 500 sampled cells per cluster**; dotplots show detection prevalence and scaled average expression. **Unplotted additions are CSV-supported, not visually validated.** This is established in `R/plots/subset_markers_plots.R:9–57`.

## Evidence inspected

Under `results/job-44542861/clusters_plots/subset/`:

- `B_cells/job-44614229/markers_res0.4/cluster_identities_res0.4.md:14–22`
- `DC/job-44614230/markers_res0.3/cluster_identities_res0.3.md:15–18`
- `Granulocytes/job-44614233/markers_res0.3/cluster_identities_res0.3.md:14–19`

Also inspected each corresponding `markers_significant.csv`, `markers_prefiltered.csv`, `heatmap_top25.pdf`, and `dotplot_top5.pdf`, plus the inherited cleanup ledger. **All six PDFs were inspected visually. No project documents were changed.** Literature validation is handled by the parent’s literature agents.

---

## CD4 and CD8 T cells

**Nested sets:** top10 = core5 + next5; top25 adds the final column. These are exploratory within-subset panels, not validated post-cleanup classifiers. **CD8 C1, C3 and C4 have explicit specificity shortfalls; do not force 25.**

| Subset / cluster | Provisional identity | Core5 | Next5 → top10 | Next15 → top25 |
|---|---|---|---|---|
| CD4 C1 | Effector/helper-memory | PRDM1, CCR5, CXCR3, BHLHE40, ITGB1 | AHR, FCGR3A, GPR183, IL2RB, TNF | CD40LG, MAF, S100A4, CAB39L, CPD, CD86, ANXA2, AHNAK, AQP3, NIBAN1, TNFRSF1B, ACTN4, GATA3, ST8SIA4, EMB |
| CD4 C2 | Naive/central-memory | LEF1, SELL, SATB1, SLC1A4, ENPP1 | CCR7, TCF7, GCNT4, IFNGR2, ICAM2 | SPTBN1, RFLNB, NAV3, SH3BGRL2, RASGRP2, NCK1, DYRK2, ZBTB10, FGD3, TXK, NISCH, LDLR, ZNF395, CDCA7L, CYBRD1 |
| CD4 C5 | Mixed activated/cycling; Treg/Tfh-enriched | TIFA, TOP2A, LGMN, DAPP1, CD80 | TNFRSF18, IRF4, SORL1, CNDP2, TACC3 | B3GLCT, HACE1, NET1, MDFIC, ACSL1, MAP3K5, TBC1D4, CDR2, ARID5A, RABEP1, PRDM2, ICOS, BCL6, ZC3H12D, SESN3 |
| CD4 C6 | Distinct state; “quiescent” provisional | MAML2, SMYD3, WWOX, CHST11, RUNX1 | THSD7B, PREX1, ATXN1, XYLT1, BCL11B | SLC9A9, PDE8A, LPP, SIK3, MSI2, KLF12, AFF3, PRKCB, SKAP1, BCAS3, BTBD9, ATP8A1, ELMO1, ANKRD44, ZFAND3 |
| CD8 C1 | GZMK effector-memory | **Four selective candidates:** GZMK, TNFSF9, CXCR3, CD8B. **GZMA only as shared fifth/context feature.** | CXCR6, CCR5, FRMD4B, PCED1B, JAML — shared contextual features | **Not supported.** Even top10 is a combination panel, not ten strongly specific markers. |
| CD8 C2 | Cytotoxic/vascular-associated effector | CX3CR1, ADGRG1, GZMB, KLF2, RASA3 | GZMH, KLF3, IL10, S1PR1, TBX21 | ADD3, PTPN4, TMEM71, FLNA, VCL, ITGB2, AHNAK, CDC42EP3, KLHL8, RIT1, RHOQ, CSGALNACT2, LBH, GNLY, PLEK |
| CD8 C3 | NCR1/KIT innate-like; lineage unresolved | NCR1, KIT, P2RY14, SPRY2, FGL2 | FCER1G, STAP1, TXK, CXXC5, HOPX | **Not supported as a high-discrimination panel.** Additional prominent markers overlap C4/C5. |
| CD8 C4 | ITGAE/ITGA1 resident-memory-like | ITGAE, ITGA1, GCSAM, CD96, EGLN3 | CCDC50, FUT8, PDE7A, DUSP4, CRTAM | **Not supported.** Further expansion mainly introduces C3-shared or generic activation genes. |
| CD8 C5 | Type-17 / MAIT-like | IL23R, RORC, IL7R, AQP3, IL18R1 | KLRB1, TNF, EOMES, SH2D1A, MAF | NFIL3, CMC1, PDE4D, GPR183, LTB, ZBTB16, PHACTR2, SMARCE1, THEMIS, PTPN6, BLTP1, SCML4, MYLIP, MCOLN2, AUH |
| CD8 C6 | Distinct state; “quiescent” provisional | THSD7B, BCL11B, MAML2, PREX1, XYLT1 | WWOX, SMYD3, CHST11, RUNX1, MSI2 | IMMP2L, ATXN1, TSPAN5, AFF3, TRPS1, SKAP1, PDE8A, ATP8A1, CDK6, BCAS3, BTBD9, SLC9A9, ANKRD44, SIK3, ZFAND3 |

**Interpretation safeguards**

- Expanded top25 lists are **candidate feature sets**, not evidence that every included gene is highly specific. CD40LG/MAF, CCR7/TCF7 and GNLY provide shared biological context.
- **Do not use** IL6ST to separate CD4 C2/C5; BACH2 to identify CD4 C6; ZNF683/XCL1/IKZF2 to distinguish CD8 C3/C4.
- CD4 C3/C4 were excluded. **CD8A alone does not invalidate porcine memory-CD4 identity.**
- NK receptors can occur on cytotoxic γδ T cells; **CD247 alone cannot establish T lineage**. CD8 C3 remains lineage-provisional; CD8 C5 remains MAIT-like, not confirmed MAIT.
- Both C6 panels identify a shared transcriptional state across CD4/CD8, not lineage-exclusive identity.

All four plots were visually inspected; unplotted additions were assessed from marker CSVs only. No project documents changed.

---

## NK and unconventional T cells

**Read-only review completed:** identity reports, cleanup ledger, marker CSVs, and visual inspection of NK/unconventional-T heatmaps and dotplots. **No project documents modified.**

**Excluded/deferred:** NK C3/C4; unconventional C8 (optional ILC2 reassignment); **all** T_cells_unresolved clusters.

Panels are **within-original-subset candidates**, not validated cleaned-subset classifiers. Columns are nested: **core 5 + next 5 = 10; + next 15 = 25**. Full rows contain 25 unique genes, checked against positive significant marker rows with positive detection differences. Supporting extensions are **not individually exclusive markers**.

### NK cells

| Cluster | Working identity / confidence | Core 5 | Next 5 → 10 | Next 15 → 25 |
|---|---|---|---|---|
| C1 | Resident/cytokine-responsive NK; good separation | GZMK, CXCR6, EOMES, KLRB1, TNFSF9 | IL18RAP, SH2D1A, SLAMF7, IL18R1, NCR1 | CLEC1B, RAB19, ERO1B, ADAM19, PLXNC1, NRP1, COTL1, SERPINF1, BCL2A1, ARHGAP18, LYST, SPINT2, STAP1, CCR5, CLEC2B |
| C2 | Mature cytolytic NK; good separation | CX3CR1, S1PR5, GZMH, GZMB, GNLY | KLRG1, PTH1R, PTAFR, ITGAM, KLF2 | RASA3, VCL, ARHGAP31, S100A4, S100A1, ADD3, CDC42EP3, PTPN4, LBH, SEMA4D, CD226, KLF3, GPR18, RASGRP2, ADGRG1 |
| C5 | Basal/low-complexity NK; **weak** | **BTC, ANXA1, TXNRD1**; AGPS and RGCC are weaker additions, not equally convincing core markers | **Not supported as a highly discriminating panel** | **Not supported; do not pad** |
| C6 | Quiescent NK/ILC-like; good state separation, identity provisional | AFF3, ZEB1, PREX1, CHST11, SMYD3 | BACH2, THSD7B, WWOX, KCNQ5, GPR39 | PIK3CD, IMMP2L, BTBD9, BCAS3, ZFAND3, ZNF609, MAML2, KLF12, SRPK2, ATP8A1, SIK3, ANKRD44, XYLT1, FRYL, TCF12 |
| C7 | TCF7-high NK/ILC1-like; good state separation, identity provisional | TCF7, SPOCK2, SLC38A1, SERPINB1, EEIG1 | P2RY14, RBPJ, NCOA7, GATA3, CD96 | BEND5, CYTH4, NDFIP2, ATP2B4, PAFAH2, ACSL4, GRK3, RPS6KA3, NIBAN1, NCK2, ERN1, ADSL, SOCS3, ANXA6, RALGDS |

### Unconventional T cells

**C1/C2 are overlapping state signatures:** their extensions are exploratory candidates, **not evidence that 25 highly specific discriminators exist**.

| Cluster | Working identity / confidence | Core 5 | Next 5 → 10 | Next 15 → 25 |
|---|---|---|---|---|
| C1 | Resting/circulating CD2− γδ T; moderate | SELL, CAPG, ID3, LY9, ITGB7 | CD27, S1PR1, RHEX, IL6R, CD52 | ANXA1, RASA3, BIN2, ADD3, TRAF3IP3, SIT1, AHNAK, TLE5, PIK3C2B, RASGRP1, SPTBN1, GRK6, RCAN3, TBC1D10C, KCTD12 |
| C2 | Early-activated/stress-adapted γδ T; moderate/weak specificity | SGK1, BMF, ZRANB1, LIX1L, SATB1 | BCL2, PDE4D, SH3BP4, CCR7, IL18R1 | ARID5B, CTSC, TLE1, MTMR12, TP53I11, WHRN, KLF9, IPCEF1, TP53INP1, PLEKHA1, NET1, TRAF4, GNA13, PIM1, IL1RL1 |
| C3 | Resident effector-memory γδ T; strong separation | TIGIT, CXCR6, CGA, CD96, CCL5 | TNF, ICOS, AQP3, PRDM1, FCGR3A | JAML, DUSP4, BHLHE40, GRK3, HLA-DRA, SLA-DQB1, CRTAM, ITGB1, IL2RB, NIBAN1, IL2RA, ID2, IL7R, ITGA4, CD74 |
| C4 | Transcriptionally poised γδ T; good state separation | TCF12, IMMP2L, WWOX, SMYD3, BCL11B | BACH2, AFF3, MAML3, PDE7B, ZEB1 | KCNQ5, MAML2, SLC9A9, FRYL, SSH2, AHCYL2, PITPNC1, XYLT1, KLF12, ZFAND3, ZNF609, PREX1, ANKRD44, CHST11, SIK3 |
| C5 | Basal/low-signal γδ T; **insufficient** | **AGPS only, provisional; no defensible core 5** | **Not supported** | **Not supported** |
| C6 | Co-stimulatory/metabolically activated γδ T; strong separation | TNFRSF9, TNFSF11, IRF4, BCL2L1, NPM3 | TNFRSF18, IL2RA, CD80, CD83, SLC1A5 | NME1, NEK6, BCL2A1, SLC29A1, SLC7A5, SRM, C1QBP, PPA1, NR4A3, NOLC1, PLAGL2, SLC35F2, RANBP1, IFRD2, TCP1 |
| C7 | CD2+ EOMES-high γδ T; strong within-compartment separation | CD2, EOMES, LCK, FCRL6, HPSE | NPY, ADGRG1, TRAT1, IKZF2, TOX | CD6, CAMK4, IL2RB, CXXC5, CRTAM, CD2AP, DGKA, GPR18, EPHB6, TGFBR3, PRKCH, SPRED2, HDAC7, RAP1GAP2, CD27 |

`HLA-DRA` is the exact supplied feature label—not evidence of human cells. The ambiguous literal feature `SLA-DRB1 (review)` was omitted.

## Numerical anchors and caveats

**“Remainder” is pooled `pct.2` from the original subset**, including subsequently excluded clusters. It is **not nearest-retained-cluster specificity**.

| Target / marker | Target detection | Pooled remainder | log2FC | Interpretation |
|---|---:|---:|---:|---|
| NK C1 / GZMK | 56.2% | 9.7% | 2.74 | Strong; NK C2 detection is 1.1% |
| NK C1 / CXCR6 | 60.1% | 14.4% | 2.27 | C2 3.5%; C7 26.3% |
| NK C2 / CX3CR1 | 35.7% | 4.9% | 2.95 | Specific but incomplete coverage |
| NK C2 / GNLY | 94.3% | 40.4% | 2.34 | Broad coverage, not exclusive |
| NK C5 / BTC | 28.1% | 3.7% | 3.46 | Best C5 candidate, low sensitivity |
| NK C5 / ANXA1 | 34.0% | 11.4% | 2.44 | Sparse enriched signal |
| NK C5 / TXNRD1 | 26.3% | 10.0% | 2.10 | Sparse enriched signal |
| NK C6 / AFF3 | 68.9% | 30.9% | 2.51 | C5 28.7%; C2 23.5% |
| NK C7 / TCF7 | 67.8% | 19.5% | 2.26 | Also expressed in excluded NK C4 |
| γδ C1 / SELL | 38.8% | 16.1% | 1.67 | Moderate; C2 18.1%, C3 9.4% |
| γδ C2 / SGK1 | 90.2% | 60.2% | 1.18 | State enrichment, not exclusivity |
| γδ C3 / TIGIT | 27.2% | 0.5% | 5.98 | Highly specific, low coverage |
| γδ C3 / CXCR6 | 30.9% | 1.0% | 5.78 | Highly specific, low coverage |
| γδ C5 / AGPS | 30.1% | 19.4% | 1.07 | Modest candidate only |
| γδ C6 / TNFRSF9 | 49.2% | 1.7% | 5.28 | Strong activation-state discriminator |
| γδ C6 / NPM3 | 69.5% | 7.1% | 4.19 | Strong within-subset blasting signature |
| γδ C7 / CD2 | 83.5% | 4.7% | 4.92 | Excellent within-subset division |
| γδ C7 / EOMES | 40.2% | 0.4% | 7.53 | Restricted, incomplete coverage |

### Main limits on specificity

- **NK C5:** BTC/ANXA1/TXNRD1 are the convincing leading candidates. AGPS/RGCC are weaker. ERC1 and GRK5 would be poor filler because detection is higher in retained C6.
- **Unconventional C5:** only **18 named positive significant genes** remained after excluding ENSSSCG/LOC features. All except AGPS had **log2FC ≤0.51**. MAP4K4, QKI, ANXA1 and other candidates are shared; five genuinely discriminating markers are unsupported.
- **Unconventional C2:** BCL2 detection is **72.7% versus 70.1% in C7**; CCR7 is **51.7% versus 67.8% in C6**. They support state interpretation, not unique identification.
- **Unconventional C1:** CD27 and S1PR1 are more prevalent in C7. ANXA1 is **54.0% in C1 versus 44.2% in C5**. These are not core-exclusive markers.
- **Unconventional C4:** BACH2 is **67.3% in C4 versus 61.8% in C6**; therefore retained as contextual support rather than a core discriminator.
- **NK C1:** NCR1 is **89.9% in C1 versus 83.0% in C7**, and CD69 is also substantial in C7. Neither distinguishes those neighbors well.

## Visual interpretation

1. **NK C1/C2:** clear opposing heatmap blocks and dotplot patterns—GZMK/CXCR6/EOMES versus CX3CR1/GZMH/S1PR5.
2. **NK C5:** a few enriched but small dots; much of its broader heatmap block overlaps other clusters. This supports the shortfall rather than a forced 25-gene list.
3. **NK C6/C7:** coherent quiescence-associated versus TCF7/SPOCK2-associated blocks, but these do not establish anatomical ILC1 identity or progenitor function.
4. **Unconventional C3/C6/C7:** strongest distinct blocks—TIGIT/CXCR6/CGA; TNFRSF9/NPM3/NME1; and EOMES/FCRL6/HPSE/LCK, respectively.
5. **Unconventional C1/C2/C5:** substantial overlap. ANXA1, SGK1 and MAP4K4 illustrate why high relative expression is not equivalent to cluster exclusivity.
6. **Unconventional C5’s exported heatmap includes negative-marker rows**, including ZFP36L1. A “top25” export is not necessarily 25 positive discriminators.
7. **NK C6 and unconventional C4 share SMYD3/WWOX/BACH2-related programs.** These distinguish states within their compartments, not NK versus γδ T globally.

Heatmaps use scaled RNA expression and at most 500 cells per cluster; dotplot color represents scaled average expression and size represents detection. **Only exported plotted features were visually assessed**; unplotted additions are CSV-supported candidates. Missing DE rows were **not** treated as zero expression.

**Lineage caution:** CD247 alone is not T-lineage evidence; NK cells also use it. Cytotoxic γδ cells can share NK receptors. NK-like T/MAIT-like labels require provisional wording; porcine memory CD4 cells can express CD8A.

### Inspected sources

Under `results/job-44542861/clusters_plots/subset/`:

- `NK_cells/job-44614238/markers_res0.4/`
- `unconventional_T_cells/job-44614239/markers_res0.4/`

For both: identity Markdown, `markers_significant.csv`, `markers_prefiltered.csv`, `heatmap_top25.pdf`, `dotplot_top5.pdf`.

Also:

- `T_cells_unresolved/job-44614240/markers_res0.4/cluster_identities_res0.4.md`
- `NOTES_BARCODE_REASSIGNMENT_AND_CLEANUP.md`

---

## Endothelial, epithelial and stromal cells

**Reviewed all three complete identity reports, master exclusion ledger, significant/prefiltered marker CSVs, and visually inspected all six heatmap/dotplot PDFs.** Excluded Endothelial C6, Epithelial C7, Stromal C6. No project files changed; rasterizations only in `/tmp`.

**Interpretation:** these are **nested candidate panels**, not validated classifiers. Top 5 emphasizes discrimination; genes 11–25 include supporting identity/state markers that are not individually specific. All listed genes were programmatically checked against their matching significant-marker CSV: present, positive log2FC, positive pooled detection difference, unique within panel.

### Endothelial

| Cluster / identity | Top 5 | Add 5 → 10 | Add 15 → 25 |
|---|---|---|---|
| C1 lymphatic | CCL21, MMRN1, PROX1, PDPN, RELN | FLT4, SEMA3D, NTN1, HS3ST1, PARD6G | TRPC6, TNFAIP8L3, DTX1, PTX3, HOXD8, HOXD9, TBX1, ACKR2, GJA1, CCN3, ABI3BP, OSBPL6, CD24, ANGPT2, LYVE1 |
| C2 activated venous | ACKR1, MATN2, LHX6, OLR1, ADGRG6 | IL6, OMD, F3, CYP2S1, PCDH17 | VCAM1, SELE, SELP, CXCL9, CXCL10, ADAMTS4, CCL2, TRAF1, PLA1A, TSPAN11, MME, ART3, GPX3, VWF, SERPINE1 |
| C3 general capillary | GPIHBP1, ANK3, HPGD, APLNR, KIT | MCAM, CD93, CCSER1, CLEC1B, TNC | TSPAN15, CDK6, ST6GALNAC3, TMEM100, CLIC5, CADM1, PLCB1, LRP11, CHST9, CLEC3B, HMCN1, NBEA, FRYL, SLC39A8, ITGA1 |
| C4 aerocyte | EDNRB, TRNP1, APLN, MFGE8, SLC16A9 | PRX, DPP4, KITLG, ITGA4, SEMA7A | CCDC68, PDZRN3, PCDH1, LMO7, ITGA3, PDGFD, RDH10, ADGRL2, USP54, NCALD, GJC1, CA4, RGCC, SIPA1L2, STXBP6 |
| C5 arterial | GJA5, BMX, SEMA3G, MFAP5, HEY2 | FN1, PDE3A, SOX17, PTGIS, EFNB2 | CXCL12, FBLN1, FBLN5, IGFBP2, ADIRF, BMPER, SSTR1, CX3CL1, COL14A1, FGF2, HTR2B, HEY1, KCNE3, CPXM2, SLC52A3 |

**C3 caution:** strongest specific markers have limited coverage: GPIHBP1 32.8%, HPGD 29.0%, APLNR 26.5%, KIT 40.9%. Combining them with higher-coverage support genes is preferable to expecting universal positivity.

### Epithelial

| Cluster / identity | Top 5 | Add 5 → 10 | Add 15 → 25 |
|---|---|---|---|
| C1 ciliated, apical-tip enriched **(weak quantitative separation)** | CFAP126, TPPP3, CFAP90, CIMIP2B, CAPSL | SNTN, CAPS, PIERCE1, RSPH14, PIH1D2 | **Not recommended:** further genes mainly reinforce shared ciliated identity, not C1 versus C5 |
| C2 distal-like secretory | ATP13A5, PREX2, ALPL, ENDOD1, C6 | PGC, ADAM28, BPIFB1, CXCL12, SFTPB | ST6GALNAC1, ARSJ, KDR, MEGF6, EDN1, RND1, LBP, ARHGAP42, P2RY2, SLC4A4, CLIC6, ESYT3, NCKAP5, ABCA5, RHOV |
| C3 canonical secretory **(insufficient discrimination)** | **No defensible highly discriminative 5-gene panel versus C2.** SCGB1A1 and SCGB3A1 are expression-intensity leads only | Not supported | Not supported |
| C4 AT2 | SFTPC, LAMP3, ABCA3, LPCAT1, SLC34A2 | NAPSA, C5, FABP3, IYD, SHISA2 | NPC2, PAH, GJB1, HHIP, DEFB1, IL13RA2, KCNE5, ETV5, DNASE2, LHFPL6, CHIA, SFTPD, KCNJ15, AQP9, LGI3 |
| C5 dynein-enriched ciliated | SLC13A1, ANKFN1, DNAH2, DNAH3, CFAP77 | PACRG, CFAP54, DNAH7, RFX3, RFX2 | HHATL, FRMPD2, SPATA18, FAM227A, NEK5, VWA3A, LRRC9, CACHD1, CFAP46, CFAP61, CFAP91, CCDC175, LRGUK, KATNAL2, MAPK10 |
| C6 AT1 | AGER, RTKN2, COL4A3, COL4A4, CAV1 | SPOCK2, VSIG2, CAV2, PALM2AKAP2, AQP5 | HOPX, PDPN, EMP2, COL4A1, COL12A1, CLDN18, PHACTR1, ARHGEF26, DLC1, UPK3B, ANKRD1, DAAM2, FLRT3, TGFB2, COL5A2 |
| C8 immune-like tuft | HMX3, SPIB, IL17RB, ALOX5AP, SOX9 | SUCNR1, NPY1R, RFLNA, GNG13, LRRC25 | CYP4A90, PTN, KIRREL3, PRSS16, IHH, HMX2, SOX8, GFI1B, POU2AF2, DNAJC12, TFPI2, CD79B, CTSK, SLC1A4, PHLDA1 |
| C9 neuronal/cholinergic-like tuft | NAV2, PLEKHA6, PRKG1, TRPM3, CHAT | DCLK2, HTR4, ERC2, SGIP1, FRMD4A | KCNB2, PAQR5, CLSTN2, RASGRF2, FSTL4, FAM184B, DGKB, THSD4, FMN1, ABHD2, RGS21, USH1C, ZSWIM5, SLIT3, TMEM163 |
| C10 ionocyte | FOXI1, ATP6V1G3, ATP6V1C2, ADGRF1, PDE1C | PGF, SLC9A2, SLC4A11, ATP7B, CFTR | TMEM61, STAP1, NCALD, SLC31A2, SLC35F3, ESRRG, PTGER3, TFCP2L1, TGFB3, LRATD1, SCNN1G, KCNMA1, SGPP2, ZACN, VEGFC |
| C11 basal | COL17A1, TP63, KRT5, KRT17, HAS3 | SFN, CLCA2, TNS4, PKP1, ITGA2 | DLK2, SNAI2, COL15A1, KCNN4, TRIM29, THBD, PLA2G7, PTPRZ1, TINAGL1, KCNQ5, SLC46A2, LTBP1, PTPRS, HS3ST1, CLDN1 |

**C3 is the most important failure case.** Its exported top-five dotplot markers are also broadly expressed in C2. Detection in **C3 versus C2**:

- SCGB1A1: **86.6% versus 88.3%**
- SCGB3A1: **60.7% versus 57.7%**
- TFF3: **95.5% versus 95.9%**
- CRISP3: **25.5% versus 27.3%**
- PTGFR: **40.4% versus 53.7%**
- S100A2: **73.9% versus 85.4%**
- WFDC2: **85.0% versus 94.9%**

A candidate **contrast signature**, rather than positive-marker panel, would be **SCGB1A1/SCGB3A1-high with ATP13A5/PREX2/ALPL/ENDOD1-low**, pending direct pairwise testing. “Distal versus proximal” remains an anatomical interpretation, not established by these plots alone.

### Stromal

| Cluster / identity | Top 5 | Add 5 → 10 | Add 15 → 25 |
|---|---|---|---|
| C1 activated matrix fibroblast **(weak/state-dependent)** | IGF1, IL6, PTGS2, VIPR2, COL6A6 | C7, CXCL12, SPON1, ADAMTS4, ADAMTSL1 | **Not recommended:** additional genes mainly shared ECM/activation support |
| C2 ductal-myofibroblast/airway-contractile mixture | WIF1, GREM1, CLDN11, NDNF, TRPC4 | BDNF, CADM1, GPX3, ADGRG6, LRP4 | WNT11, PCDH17, MYH11, ACTA2, LMOD1, CNN1, PDE5A, SYNPO2, COL4A5, PTGS1, MAOB, ROBO2, P3H2, TNC, MET |
| C3 alveolar fibroblast | NPNT, BMPER, TF, FMO2, GYG1 | A2M, COL13A1, ITGA8, CMKLR1, COL8A1 | ADRA1A, LIMCH1, SHISA2, CHI3L1, NAV2, UNC5C, OLFML2A, RELL1, MYO1B, RSPO2, ENPEP, FIBIN, TRABD2B, LPL, PRKG1 |
| C4 adventitial-like fibroblast | ADAMTS16, DKK2, ACKR2, MFAP5, CLEC3B | PROCR, C1QTNF3, ACKR4, ACKR3, DPP4 | WNT10B, SLC7A10, DUOX2, DUOXA2, IGFBP5, ADAMTS5, PCOLCE2, POSTN, IL33, HTRA3, CADM3, HMCN2, SEMA3B, ENPP1, ITGA11 |
| C5 epithelial-niche-like fibroblast | LGR5, CRABP1, SFRP1, FST, HGF | FGF10, SULF1, MAMDC2, FBN2, FGFR2 | SFRP2, COL15A1, ANGPT4, MECOM, ABCA8, AQP3, USP54, PLA2R1, MATN2, TGFBI, CXCL12, HEYL, PTN, DTX4, PRRX1 |
| C7 vascular mural mixture | TINAGL1, NOTCH3, GUCY1B1, MCAM, EDNRB | PDE3A, ITGA7, FRZB, TRPC6, NAALAD2 | FHL5, SHISAL1, PDGFA, ANGPT2, NEURL1B, ITGA1, RCAN2, BCAM, KCNAB1, DES, ABCC9, KCNJ8, COX4I2, EDNRA, RASL12 |

### Eight plot-driven findings

1. **Endothelial heatmap:** strong C1 and C5 blocks; capillary C3/C4 share more signal. C3’s exported high-fold markers have small dots—specificity does not imply coverage.
2. **Endothelial venous versus arterial:** heatmap shows an inflammatory/adhesion band shared by C2 and C5. SELE detection is **59.1% versus 62.4%**, SELP **65.0% versus 63.2%**, VCAM1 **83.3% versus 67.7%**. These were deliberately excluded from C2’s top five.
3. **Epithelial ciliated split:** heatmap and dotplot show a strong shared C1/C5 ciliated program. CFAP126 separates only quantitatively (**80.1% versus 59.6%**), whereas C5 ANKFN1 (**89.8% versus 54.1%**) and DNAH2 (**76.8% versus 40.6%**) offer larger nearest-cluster detection contrasts.
4. **Epithelial secretory split:** C2/C3 top markers share bright blocks. C3 lacks five convincing positive discriminators; C2 has the more useful asymmetric signature.
5. **Tuft split:** shared POU2F3/SH2D7/ALOX5 identity is not discriminatory. C8 SPIB (**91.1% versus 41.2% C9**) and IL17RB (**92.9% versus 53.7%**) help. C9 NAV2 (**98.8% versus 51.6% C8**) helps. **RGS21 does not discriminate by detection** (**86.3% versus 87.5%**); its separation is intensity-based, hence only in the expanded panel.
6. **Stromal C1/C4:** top C1 heatmap rows extend into C4; C1 state genes are relatively sparse. COL14A1 (**97.1% versus 99.2%**) and LUM (**93.0% versus 94.1%**) are poor discriminators despite their annotation prominence.
7. **Stromal C3/C5:** heatmap shows shared COL13A1/ITGA8 modules. These identify the alveolar fibroblast axis, not uniquely C3. Also, A2M has a conspicuous large dot in C7, so it was removed from C3 core five. **HEYL is stronger in C7 than C5** (**66.5% versus 50.7%**), so it was removed from C5 core five.
8. **Stromal contractile split:** C2 and C7 share MYH11/contractile signal, but WIF1/GREM1/CLDN11 distinguish C2 and TINAGL1/NOTCH3/GUCY1B1 distinguish C7. The panels distinguish the **current mixed clusters**; they cannot resolve myofibroblast versus airway smooth muscle or pericyte versus vascular smooth muscle within them.

### Numerical anchors and source locations

Here `pct.2` is the **original pooled other clusters**, including excluded clusters—not cleaned nearest-cluster specificity.

| Subset / cluster / gene | pct.1 / pct.2 | log2FC | CSV line |
|---|---:|---:|---:|
| Endothelial C1 MMRN1 | 93.3% / 4.6% | 6.324 | 9 |
| Endothelial C2 ACKR1 | 89.8% / 1.7% | 6.790 | 2456 |
| Endothelial C3 GPIHBP1 | 32.8% / 2.4% | 4.725 | 3854 |
| Endothelial C4 TRNP1 | 74.1% / 1.6% | 6.742 | 5742 |
| Endothelial C5 GJA5 | 85.4% / 4.5% | 5.198 | 7766 |
| Epithelial C2 ATP13A5 | 60.9% / 5.0% | 3.538 | 2912 |
| Epithelial C4 SFTPC | 99.5% / 7.1% | 8.141 | 7290 |
| Epithelial C8 SPIB | 91.1% / 2.4% | 5.845 | 16788 |
| Epithelial C9 NAV2 | 98.8% / 11.6% | 6.503 | 19156 |
| Epithelial C10 FOXI1 | 91.8% / 0.2% | 9.745 | 21761 |
| Stromal C2 WIF1 | 81.1% / 3.3% | 7.041 | 1404 |
| Stromal C4 ADAMTS16 | 62.7% / 2.0% | 6.137 | 4567 |
| Stromal C5 CRABP1 | 57.1% / 6.8% | 4.126 | 5753 |
| Stromal C7 TINAGL1 | 67.0% / 1.0% | 7.444 | 7890 |

Sources, all under `results/job-44542861/clusters_plots/subset/`:

- `Endothelial/job-44614231/markers_res0.3/`
- `Epithelial/job-44614232/markers_res0.4/`
- `Stromal/job-44614235/markers_res0.4/`

Within each, inspected `cluster_identities_res*.md`, `markers_significant.csv`, `markers_prefiltered.csv`, `heatmap_top25.pdf`, and `dotplot_top5.pdf`. Nearest-cluster percentages above come from actual CSV rows; missing rows were **not** treated as zero.

**Limitations:** only exported genes can be visually assessed in these existing plots. Added panel genes were assessed from CSV statistics, not visually verified in an ungenerated new plot. No pairwise cleaned-cluster DE or donor-level validation was run.

---

## Verified literature anchors and limits

These previously verified abstract-level references support interpretation, not all individual panel genes.

| PMID | Scope | Narrow implication |
|---|---|---|
| [32535044](https://pubmed.ncbi.nlm.nih.gov/32535044/) | Pig, CLEC12B phenotyping | PAM enrichment; lower DC expression prevents absolute exclusivity. |
| [25896759](https://pubmed.ncbi.nlm.nih.gov/25896759/) | Pig, ENPP1/CD203a antibody identification | Strong PAM enrichment relative to monocytes/granulocytes in that assay. |
| [23667115](https://pubmed.ncbi.nlm.nih.gov/23667115/) | Pig monocyte phenotyping/microarrays | Maturation-series interpretation; CD163 alone is not a MoMac/PAM discriminator. |
| [29977043](https://pubmed.ncbi.nlm.nih.gov/29977043/) | Pig AM/PIM comparison | Related AM/PIM programmes limit transcript-only anatomical assignment. |
| [30158925](https://pubmed.ncbi.nlm.nih.gov/30158925/) | Pig DC phenotype/function | XCR1/CADM1-high supports cDC1. |
| [36203585](https://pubmed.ncbi.nlm.nih.gov/36203585/) | Pig multicolour flow | Memory CD4 cells can express CD8alpha. |
| [22539298](https://pubmed.ncbi.nlm.nih.gov/22539298/) | Pig NK phenotype/function | NKp46-negative functional NK cells exist. |
| [39144143](https://pubmed.ncbi.nlm.nih.gov/39144143/) | Pig gamma-delta phenotyping | CD2-positive gamma-delta cells can share NK receptors/cytotoxic programmes. |
| [34030460](https://pubmed.ncbi.nlm.nih.gov/34030460/) | Human/mouse endothelial atlas | Supports aerocyte/vascular framework; cross-species extrapolation is qualified. |
| [40467553](https://pubmed.ncbi.nlm.nih.gov/40467553/) | Human airway atlas/differentiation | Tuft-like states can include tuft-ionocyte progenitors; porcine tuft subtypes remain provisional. |
| [39806101](https://pubmed.ncbi.nlm.nih.gov/39806101/) | Mouse/human pericytes | Some pericytes share smooth-muscle genes; contractile genes alone do not resolve stromal C7. |
| [37465668](https://pubmed.ncbi.nlm.nih.gov/37465668/) | Human EVLP scRNA-seq | Inflammation changes during EVLP; separate state from identity. Pathway evidence does not validate specific genes. |

### Citation/annotation cautions carried forward

The existing ledger's “Shin et al. 2025” attribution for PMID 39875369 does not match the retrieved record (Lim et al.), and its claimed NREP/PON1 PAM subtype definitions were not established by the reviewed abstract. The Schupp endothelial atlas is PMID **34030460**, not 34097468. Existing source documents were not corrected in this task.

### Next analytical step

Pairwise marker evaluation against the closest retained cluster, donor-aware replication and donor-held-out classification after cleanup. Do not force a five-gene panel for macrophage C3, B C3/C4, epithelial C3 or unconventional-T C5; NK C5 has only three convincing leads.
