# Proposed cluster identities — refit job-44890936, resolution 0.4

Proposal derived from `markers_top25.csv` per cluster (FindAllMarkers Wilcoxon,
`only.pos=FALSE`, min.pct 0.25, logfc.threshold 0.25, padj < 0.05; 50,512
significant rows over 22 clusters), checked against the full
`markers_significant.csv` for lineage contrasts (positive and negative
enrichments; all top-25 rows are positive by descending-avg_log2FC rank),
the lineage panels in `R/clusters/markers.R`, and the pig lung marker
literature (`results/literature_review/SCRATCHPAD.md`). Reviewed with the
oracle before writing; six wordings were revised per review (clusters 1, 3, 7,
14, 17, 20) and the methods line was corrected from `only.pos` to the
run-recorded `only.pos=FALSE` (see `markers_res0.4/.INFO`).

Parent run: job-44890936 (refit of `clean_annotated_object.rds`; 234,183 cells;
HVG 3000; dims 1–30; Harmony max_iter 30 with early stopping; Leiden 0.4 gives
22 clusters). Cluster numbers here are this run's own labels and do not
correspond to job-44542861 res 0.7 clusters. "Not returned" below means not a
significant positive marker under these filters; it is not absence of
expression. TRDC is absent from the BD panel, so gamma-delta receptor identity
is untestable. CD5/CD6 non-return does not exclude T identity in pig
(Mair 2016, CD3+NKp46+ cells can lack CD5/CD6).

| Cl | Cells | Proposed identity | Key evidence |
|---|---|---|---|
| 1 | 56,057 | T cells, CD8A-enriched, naive/central-memory-like program; subtype unresolved | CD3D 65%, CD3E 64%, CD3G 49%, CD2 74%, CD5 59%, CD6 38%, CD8A 41%, CD28 39%, IL7R 46%, TCF7 43%, CCR7 26%, LCK/LAT/ITK, CD247 61%; GZMA 31%, CCL5 48%, CRTAM 36%, RUNX3 33%; CD4/CD8B not returned positive; CD40LG/ICOS/LEF1/SELL not returned |
| 2 | 27,426 | Macrophages, FCN2⁺ CLEC12B⁺ S100-high (alveolar-like, provisional) | FCN2 97%, CLEC12B 53%, SIGLEC1 77%, S100A12 98%, LYZ 76%, LGALS3 99%, APOE 99%, FCGR1A 85%, CD163 97%, MRC1 94%, PPARG, GPNMB; FCN1 83%, MAFB, CD68, TLR4 |
| 3 | 27,219 | NK-like cytotoxic lymphocytes; NK/T composition unresolved | NCR1 52%, KLRB1 55%, XCL1 60%, EOMES 27%, KLRK1 91%, GZMA 81%, GZMB 35%, GNLY 52%, PRF1 68%, IL2RB 94%, CCL5 98%, CRTAM 51%, TBX21 51%, CXCR6 30%, NKG7 74%, CD244 31%; CD3D 65%, CD2 85%, CD3G 37%, ITK 52%, LCK 45%, CD8A 37%; CD5/CD6 not returned positive |
| 4 | 18,414 | Neutrophils, with interferon-responsive program | CSF3R 90%, CXCL8 90%, G0S2 66%, RETN 67%, TREM1 60%, GPR84 69%, MCEMP1 69%, SELL 68%, CD14 56%, C5AR1 76%, TLR4 61%; IFIT1/IFIT2/OAS1/OAS2/MX1/ISG15 enriched within this cluster |
| 5 | 17,894 | GATA3⁺ SOX13⁺ T cells (γδ-like hypothesis, receptor subtype unresolved) | CD3D 74%, CD3E 70%, CD247 84%, TCF7 79%, GATA3 84%, SOX13 82%, BLK 59%, RHEX 84%, LEF1 36%, CCR7 38%, CD27 27%; CD2 negatively enriched (log2FC −4.15, 3.7% vs 34.5%); no CD4/CD8A/CD8B/CD5/CD6 returned positive |
| 6 | 16,021 | B cells | MS4A1 81%, CD19 49%, CLEC17A 36%, PAX5 64%, CD79B 89%, CD79A 47%, BANK1 72%, CXCR5 39%, TNFRSF13C 37%, FCRL3 36%, EBF1 59% |
| 7 | 10,787 | Plasma cells, IgM-enriched | JCHAIN 98%, TNFRSF17 73%, MZB1 71%, IGHM 81%, POU2AF1 77%, DERL3 46%, TXNDC5 88%, PRDM1 88%, XBP1 82%; no IGHA/IGHG/IGK/IGL returned positive |
| 8 | 8,231 | Fibroblasts, PI16⁺ adventitial-like | PI16 76%, DPT 82%, LUM 88%, CHRDL1 84%, FBLN1 96%, COL1A1 100%, COL3A1 99%, TCF21 78%, CD34 57%, PDPN 69% |
| 9 | 7,772 | Blood vascular endothelium; vascular subtype unresolved | ACKR1 33%, SELE 34%, SELP 35%, VWF 66%, PECAM1 99%, CLDN5 85%, ESAM 69%, CDH5 87%, CD34 69%, KDR 38%, FLT1 49%, TEK 64%, NOS3 40%, SOX17 37%, DLL4 36%, NOTCH4 57%, EPHB4 34%, CXCL9 39%, IL33 43% |
| 10 | 7,601 | IL1B-high CD14⁺ FCN1⁺ inflammatory monocytes | IL1B 74%, FCN1 88%, CD14 83%, CSF1R 66%, VCAN 77%, MAFB 82%, F13A1 93%, AREG 68%, EREG 52%, PLBD1 76%, MS4A7 27%, C5AR1 82%, S100A8 48% |
| 11 | 6,064 | Lymphatic endothelium | CCL21 86%, MMRN1 95%, FLT4 69%, RELN 81%, LYVE1 77%, PROX1 54%, PDPN 90%, PECAM1 98%, CDH5 84% |
| 12 | 5,352 | Macrophages, C1Q-high CCL8/CCL14⁺ (interstitial-like, provisional) | CCL8 76%, CCL14 87%, C1QA 96%, C1QB/C 95%, MS4A7 58%, FCGR2B 86%, STAB1 73%, MSR1 83%, CD163 96%, MRC1 95%, CSF1R 69%, MAFB 88%, CCL2 70%, DNASE1L3 37%, PLA2G2D 25% |
| 13 | 4,616 | Eosinophils (moderate confidence) | EPX log2FC 11.01, 27%, CCR3 32%, ALOX5 72%, CLNK, ADGRG3, EMB, MCTP1, RGS18; large fold change does not establish whole-cluster homogeneity |
| 14 | 4,535 | pDC; tDC contribution unresolved (NCR1 inspection trigger, see caveat) | TCF4 95%, IRF8 99%, PLD4 84%, SPIB 75%, BLNK 98%, TLR7 67%, TLR9 31%, TLR3 59%, CD4 88%, VPREB1 77%, CLEC12A 68%, FLT3 68%, BCL11A 93%, IRF7 48%, RUNX2 81%, LAMP3 91%, IRF4 71%, CD8B 52%; NCR1 60% (see caveat) |
| 15 | 3,695 | Secretory/club-like epithelium with surfactant-associated transcripts | SCGB1A1 79%, TFF3 91%, SCGB3A1 52%, SCGB3A2 78%, BPIFB1 49%, PIGR 82%, WFDC2 90%, KRT7 65%, KRT19 78%, SFTPB 66%, SFTPD 70%, SFTPA1 79%, NAPSA 45%; co-expression does not establish transition |
| 16 | 3,051 | Ciliated epithelium | CIMIP1 91%, SNTN 75%, RSPH4A 87%, FOXJ1 77%, CAPS 77%; PDPN 69%, HOPX 83%, AGER 38% (substantial co-detection, not "low-level"; ciliated identity kept) |
| 17 | 2,880 | cDC1, XCR1⁺ with MS4A2 enrichment | XCR1 40%, FLT3 75%, CADM1 49%, BATF3 48%, IRF8 72%, WDFY4 50%, CD83 89%, CCR7 38%, MS4A2 65%, KIT 27%; FCER1A/CD207/CLEC10A/CD1C/SIRPA/ITGAM not returned positive; IRF4 weak (27%); no distinct cDC2 at this resolution |
| 18 | 2,832 | Smooth muscle cells; pericyte contribution unresolved | MYH11 84%, CNN1 45%, ACTG2 40%, ACTA2 77%, TAGLN 89%, MYLK 73%, GREM1 52%, WIF1 80% (identity weighted on contractile program, not GREM1/WIF1 alone) |
| 19 | 1,291 | AT2 cells | SFTPC 100%, SLC34A2 90%, SFTPA1 100%, SFTPD 96%, NAPSA 93%, LAMP3 96%, CHIA 47%; AGER 80% recorded without inferring AT2→AT1 transition |
| 20 | 1,064 | MARCO⁺ TIMD4/VSIG4-enriched macrophages (provisional; endothelial-signal and homogeneity unresolved) | MARCO 88%, TIMD4 48%, C4BPA 92%, C1QA 92%, C1QB 93%, CD163 85%, CSF1R 48%, VSIG4 64%, SIGLEC1 40%, GPNMB 77%, HMOX1 63%; PLVAP 82%, CDH5 40% (see caveat); FCN2/S100/LYZ/APOE not returned |
| 21 | 714 | Tuft cells (POU2F3⁺ chemosensory epithelium) | POU2F3 95%, GNG13 72%, RGS21 88%, SH2D7 96%, IL17RB 76%, CHAT 34%, POU2AF2 49%, AVIL 84%, SOX9 66%, RGS13 84%, HPGD 99%, LTC4S 74%, EPCAM 78%, KRT8 90%, KRT18 85%; KIT 83% compatible with tuft identity (see note); mast proteases/HDC not returned; MS4A2 not positively enriched here |
| 22 | 667 | AT1 cells | AGER 98%, HOPX 96%, PDPN 90%, COL4A3 74%, COL4A4 77%, CLDN18 83%, RTKN2 84%, AQP5 78% |

## Caveats (oracle-agreed wording)

- Cluster 1: T-cell identity is strongly supported, but CD4/CD8 composition,
  receptor subtype and differentiation state remain unresolved. CD8A enrichment
  does not establish a CD8-only population in pigs. IL7R/TCF7/CCR7/CD28 support
  a naive/central-memory-like program, not a uniformly naive cluster. CD4 and
  CD8B were not returned as significant positive markers; this does not
  establish absence of their expression or of CD4 T cells.
- Cluster 3: the cluster has a strong NK-like cytotoxic program, but significant
  CD3D/CD3G and T-cell signaling transcripts leave conventional NK versus
  NK-like T-cell composition unresolved. CD5/CD6 were not returned as
  significant positive markers, which neither demonstrates their absence nor
  excludes T identity. Porcine CD3+NKp46+ T cells can lack CD5/CD6 (Mair 2016).
  Neither NKT nor MAIT identity is established by these markers.
- Cluster 14: the broad TCF4/IRF8/PLD4/SPIB/BLNK/TLR7/9/CD4 program supports
  pDC identity. IRF4 enrichment alone does not establish transitional DC
  identity; a tDC contribution remains unresolved. NCR1 is detected in 60.3% of
  cells and requires cell-level inspection of expression strength and
  accompanying lineage programs. This signal does not, by itself, justify
  NK/NKT reassignment, and cluster homogeneity has not been established.
- Cluster 17: XCR1 with FLT3/CADM1/IRF8/BATF3/WDFY4 supports cDC1 identity.
  MS4A2 and KIT enrichment do not establish mast-cell identity, and MS4A2 must
  not be substituted for FCER1A as cDC2 evidence. CCR7/CD83 are consistent with
  a maturation-associated program but do not demonstrate migration. No distinct
  cDC2 population is identified at this resolution; its presence remains
  unresolved.
- Cluster 20: a coherent MARCO/TIMD4/VSIG4/C1Q/CD163/CSF1R program supports a
  provisional macrophage assignment. PLVAP (82.4%) and CDH5 (40.2%) provide
  substantial endothelial-associated signal whose cellular origin is unresolved.
  Prior doublet filtering does not exclude residual macrophage–endothelial
  doublets; ambient RNA, shared expression and cluster heterogeneity also remain
  possible. Homogeneity has not been established, and neither a doublet identity
  nor endothelial transition is assigned.
- Cluster 21 KIT note: KIT enrichment is compatible with tuft identity and is
  not treated as evidence of mast contamination. The assignment rests primarily
  on POU2F3/GNG13/RGS21/SH2D7/IL17RB/CHAT/POU2AF2 and epithelial expression
  (EPCAM/KRT8/KRT18). Mast proteases CPA3/TPSAB1/TPSB2 and HDC were not
  returned as significant markers, and MS4A2 was not positively enriched in
  this cluster; these results do not prove absence of a minor mast population.
  Tuft-intrinsic KIT expression in this pig dataset remains an interpretation
  (Lara et al., Sci Adv 2026, mouse tuft KIT, PMID 41739917) rather than an
  independently validated finding.
- No distinct CD4+ T cluster is identified at 0.4 (CD4 transcript not returned
  positive in any T cluster). CD4/CD8 deferral to subclustering is provisional:
  anchor T lineage with CD3D/E/G, LCK, LAT, ITK; examine CD4 vs CD8A/CD8B within
  that compartment (porcine CD4 cells can express CD8A); treat
  TCF7/LEF1/CCR7/SELL/IL7R/CD28 as differentiation-state, not CD4, evidence.
  A parent-to-refit barcode contingency table is still the smallest check of
  whether parent CD4/CD8 populations merged. Subclustering cannot recover the
  unmeasured TRDC receptor.
- No distinct cDC2 population is identified (FCER1A/CD207/CLEC10A/CD1C not
  returned positive anywhere). Inspect FLT3/MHC-II with SIRPA, FCER1A, CD207,
  IRF4, FCGR2B against XCR1/CADM1/BATF3/IRF8 and the pDC program, including
  neighboring monocyte/DC populations, rather than presuming cDC2 sits inside
  cluster 17. Do not use human CD1C/CLEC10A non-return as the decisive
  exclusion.
- No separate IFN-high mixed-lineage cluster is assigned at 0.4;
  IFN-responsive transcripts (IFIT1/IFIT2/OAS/MX1/ISG15) are enriched within
  neutrophils (cluster 4) and endothelium (cluster 9). Anatomical qualifiers
  (alveolar-like, interstitial-like, adventitial-like) are hypotheses, not
  established locations.

## Proposed pre-subclustering groups

Groups below are working compartments for subclustering, not final annotations.
They roughly recover the 11 cleaned compartments (see
`data/job-44542861/clean_annotated_object.INFO`); exact correspondence needs a
barcode contingency table, not assumed from similar markers or sizes.

- Lymphoid T/NK/ILC: clusters 1, 3, 5 (CD8A-enriched T; NK-like cytotoxic;
  GATA3+SOX13+ T). Priority: CD4 vs CD8, NK vs NK-like T, TRDC still
  untestable; inspect CD4/CD8A/CD8B, CD40LG/ICOS, RUNX3/cytotoxic program.
- B/plasma: clusters 6, 7.
- Myeloid phagocytes: clusters 2, 10, 12, 20 (FCN2-S100 macrophages;
  inflammatory monocytes; C1Q-CCL macrophages; MARCO macrophages). Priority:
  cluster 20 homogeneity/endothelial signal; monocyte–macrophage boundaries
  (Baillou 2025 DC3 in-blood caveat: CD14-low DC3 can fall in cDC2 gates).
- Granulocytes: clusters 4 (neutrophils), 13 (eosinophils).
- Dendritic cells: clusters 14 (pDC, NCR1 inspection), 17 (cDC1; cDC2
  unresolved, inspect neighbors too).
- Endothelial: clusters 9 (blood vascular), 11 (lymphatic).
- Epithelial: clusters 15 (secretory/club-like), 16 (ciliated), 19 (AT2),
  21 (tuft), 22 (AT1). Priority: surfactant/AGER overlaps recorded above are
  not transition evidence.
- Stromal: clusters 8 (fibroblasts), 18 (smooth muscle; pericyte unresolved).
