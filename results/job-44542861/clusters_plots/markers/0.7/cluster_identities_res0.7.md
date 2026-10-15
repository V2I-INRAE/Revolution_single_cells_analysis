# Proposed cluster identities — job-44542861, resolution 0.7

Proposal derived from `markers_top25.csv` per cluster (FindAllMarkers Wilcoxon, only.pos, min.pct 0.25, logfc 0.25, padj < 0.05; 29,660 significant rows over 24 clusters), checked against the full `markers_significant.csv` for lineage contrasts, the lineage panels in `R/clusters/markers.R`, and the pig lung marker literature (`results/literature_review/SCRATCHPAD.md`, section "Marker -> publication map (job-44493014, res 0.7)"). Reviewed with the oracle; three labels were downgraded per review (clusters 8, 12, 21) and "absent/negative" wording replaced with "not returned as significant positive markers".

Parent run: job-44542861 (scDblFinder doublet detection; 238,627 cells; 21,123 genes; Harmony dims 1-30; Leiden 0.7 gives 24 clusters). Cluster numbers here are this run's own labels and do not correspond to job-44493014 clusters. "Not returned" below means not a significant positive marker under these filters; it is not absence of expression. TRDC/TRGC1 are absent from the BD panel, so gamma-delta identity is untestable.

| Cl | Cells | Proposed identity | Key evidence |
|---|---|---|---|
| 1 | 27,448 | Macrophages, FCN2⁺ CLEC12B⁺ S100-high (alveolar-like, provisional) | FCN2 97%, CLEC12B 53%, SIGLEC1 77%, S100A12 98%, LYZ, LGALS3 99%, APOE 99%, FCGR1A, CD163, MRC1 |
| 2 | 25,888 | CD4⁺ T cells, naive/central-memory-like | CD28, CD40LG, ICOS, LEF1, IL7R 68%, CD4 27%, CD6, TCF7 64%, CCR7 51%, ITK; naive vs central-memory unresolved |
| 3 | 24,342 | NK-like cytotoxic lymphocytes (T contribution unresolved) | NCR1 55%, KLRB1, XCL1 61%, EOMES, KLRK1 92%, GZMA 85%, PRF1, IL2RB 95%, CD244 |
| 4 | 20,852 | GATA3⁺ SOX13⁺ T cells (γδ-like hypothesis, receptor subtype unresolved) | CD3D 73%, CD3E 70%, CD247 83%; RHEX 79%, SOX13 78%, BLK, GATA3 82%; no CD4/CD8A/CD8B/CD2/CD5/CD6 returned; TRDC absent from panel |
| 5 | 20,552 | CD8A⁺ cytotoxic T cells (αβ/γδ unresolved) | CD8A 59%, GZMB/H/K, GZMA 61%, PRF1, GNLY, CCL5 90%, CRTAM, CXCR6, RUNX3, TBX21; KLRK1 58% |
| 6 | 18,411 | Neutrophils | CSF3R 90%, CXCL8 90%, G0S2, RETN, TREM1, GPR84, MCEMP1; IFIT1/IFIT2/OAS2 enriched within this cluster |
| 7 | 16,010 | B cells | MS4A1, CD19, CLEC17A, PAX5, CD79B 89%, CD79A, BANK1 72%, CXCR5, TNFRSF13C, FCRL3/5, EBF1 |
| 8 | 11,935 | T cells, CD8A-enriched; subtype unresolved (lowest confidence) | CD3D/E/G 41-53%, CD2 60%, CD5, CD6, LCK, LAT; CD8A 39% (no CD8B, no CD4); GZMA/PRF1 only, CCL5 64%, CCR5, RUNX3, CRTAM, KLRK1 32%; max top-25 log2FC 2.01; no CCR7/TCF7/LEF1/CD27/CD28, exhaustion, MAIT/NKT, ASC or B program returned |
| 9 | 10,780 | Plasma cells (IgM-dominant) | JCHAIN 98%, TNFRSF17 73%, MZB1 71%, IGHM 81%, POU2AF1, DERL3, TXNDC5, IGLL1 (review), CLGN |
| 10 | 8,509 | Fibroblasts (PI16⁺ adventitial-like) | PI16 74%, DPT, LUM 86%, CHRDL1, FBLN1 96%, COL1A1/2, COL3A1 98%, TCF21, AGTR1, C7, COL14A1 |
| 11 | 8,009 | IL1B-high CD14⁺ FCN1⁺ inflammatory monocytes | IL1B 71%, PLBD1, AREG, F13A1 89%, EREG, VCAN, CSF1R, FCN1 86%, CD14 80%, MAFB, MS4A7 |
| 12 | 7,986 | Blood vascular endothelium; vascular subtype unresolved (provisional) | ACKR1 32%, SELE 34%, SELP, VWF, PECAM1 98%, CLDN5, ESAM, CDH5, CD34, KDR/FLT1/TEK, EMCN, EFNB2, SEMA3G, SOX17, NOS3, DLL4/NOTCH4/EPHB4, EDNRB/AQP1/PRX, CXCL9/10, IL33, CX3CL1; GJA5/HEY1/HEY2 not returned; mixed arterial/venous/capillary signals at this resolution |
| 13 | 6,168 | Lymphatic endothelium | CCL21 86%, MMRN1 95%, FLT4 68%, RELN 80%, LYVE1 77%, TBX1, HS3ST1, PROX1 54%, PDPN |
| 14 | 5,378 | Macrophages, C1Q-high CCL8/CCL14⁺ (interstitial-like, provisional) | CCL8 76%, CCL14 87%, C1QA 96%, MS4A7, FCGR2B 86%, STAB1, MSR1, DNASE1L3, PLA2G2D |
| 15 | 4,626 | Eosinophils | EPX (log2FC 10.8) 27%, CCR3 32%, ALOX5 72%, CLNK, ADGRG3, EMB, MCTP1, RGS18 — moderate confidence |
| 16 | 4,542 | pDC (IRF4⁺; tDC contribution unresolved) | VPREB1 77%, BLNK 98%, TCF4 95%, IRF8 99%, CD4 88%, CD8B 52%, TLR3/7/9, PLD4, CLEC12A, SMIM5, LAMP3 91%, BCL11A, FLT3, RUNX2, IRF4 |
| 17 | 3,845 | Secretory/club-like epithelium with surfactant-associated transcripts | SCGB1A1 80%, TFF3 91%, SCGB3A1, BPIFB1, PIGR, WFDC2, KRT7/19, SCNN1B/G + SFTPB 64%, SFTPD 69%, NAPSA 43%; co-expression does not establish transition |
| 18 | 3,351 | Ciliated epithelium | CIMIP1 86%, SNTN, RSPH4A 81%, RSPH1, TEKT2, CAPS 71%, FOXJ1 72%; AGER/PDPN low-level co-detection caveated |
| 19 | 2,889 | cDC1, XCR1⁺ (MS4A2⁺) | XCR1 40%, TIFAB, FLT3 75%, MS4A2 65%, CD1.1 (review), HLA-DRA/SLA-DQB1, CD83, NOTCH4, CADM1, IRF8, BATF3 |
| 20 | 2,894 | Smooth muscle cells | GREM1 51%, WIF1 79%, CNN1, ACTG2, MYH11 83%, ACTA2 76%, LMOD1, NDNF, TAGLN 89%, MYLK; PDGFRB present, pericyte split unresolved |
| 21 | 1,379 | MARCO⁺ macrophages, TIMD4/VSIG4-enriched (provisional) | MARCO 71%, TIMD4 38%, C4BPA, SAA3, IDO1, CCL19, CD180, C1QA/B/C 76-80%, VSIG4 51%, HMOX1, GPNMB, SIGLEC1 32%, CSF1R, CD14, CD163 78%; FCN2/S100/LYZ/APOE not returned; PLVAP 69% + CDH5 33% endothelial-signal caveat, homogeneity unchecked |
| 22 | 1,335 | AT2 cells | SFTPC 100%, SLC34A2 89%, SFTPA1 100%, SFTPD 95%, NAPSA 92%, LAMP3 95%, CHIA, IYD |
| 23 | 753 | AT1 cells | AGER 95%, COL4A3/4, CLDN18 77%, RTKN2 79%, AQP5 73%, HOPX 94%, PDPN 85%, SCEL, ARHGEF26 |
| 24 | 745 | Tuft cells (POU2F3⁺ chemosensory epithelium) | POU2F3 93%, GNG13 70%, RGS21 85%, SH2D7 94%, IL17RB 74%, CHAT 33%, POU2AF2, HMX2/3, SOX9 |

Notes: no separate IFN-high mixed-lineage cluster was identified among the 24 clusters at resolution 0.7 in this run; IFN-responsive transcripts (IFIT1/IFIT2/OAS2) are enriched within the neutrophil and endothelial clusters. This does not establish the fate of the previous run's mixed cluster (no barcode correspondence checked). Anatomical qualifiers (alveolar-like, interstitial-like, adventitial-like, resident) are hypotheses, not established locations. EPHB4 groups with venous-associated evidence (EFNB2 is the arterial contrast); BMX is significant in cluster 24, not absent dataset-wide.
