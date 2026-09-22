# Week 5 Homework 2 — RNA-seq Analysis Using EasyMultiProfiler-Web

**Course:** Bioinformatics: From Multi-Omics Data to Discovery
**Week:** 5 — Transcriptomics (RNA-seq Principles and Differential Expression)
**Student:** 崔澜馨 (SUAT24000118) · EMP 0.2.8.9
**Dataset:** `tests/RNAseq_output.csv` (24 samples × 19,150 genes, mouse, Ensembl/RIKEN IDs)
**Primary contrast:** T4400 vs DMSO (n = 4 vs 4) — chemical perturbation main effect
**Toolchain:** EMP-Web (EasyMultiProfiler-Web) for upload & sync; DESeq2 1.46.0 + apeglm 1.28.0 + clusterProfiler 4.14.6 + enrichplot 1.26.6 for DE / GSEA.

---

## 1. Workflow

1. **Import** — `POST /api/import/demo?dataset_id=rnaseq_course` uploaded the bundled `tests/RNAseq_output.csv` (counts) and `tests/RNAseq_mapping.csv` (Group metadata, 6 groups × 4 samples) as `session_id=...rnaseq_course`.
2. **QC + filter** — Genes with ≥ 10 counts in ≥ 3 samples retained (`rowSums(counts >= 10) >= 3`), leaving **13,923 / 19,150 genes**.
4. **Differential expression** — DESeq2 with `design = ~ Group`, reference = DMSO. LFC shrunken with **apeglm** for robust effect-size estimation.
5. **Pathway / GO** — `gseGO(ont = "BP")` and `gseKEGG(organism = "mmu")` on a per-gene signed `-log10(p)` rank metric, mapped SYMBOL → ENTREZ via `org.Mm.eg.db` 3.20.0 (13,068 / 13,923 mapped).
6. **Visualization** — VST PCA (all 24 samples), volcano (apeglm-shrunken LFC), DEG heatmap (VST z-scores), `dotplot()` for GO/KEGG.
7. **Sync** — GitHub push to `EMP2026/Week_05/transcriptomics/weekly/runs/<run_id>` via `/api/github/sync`.

---

## 2. DE results (T4400 vs DMSO)

| Threshold | DEGs |
|---|---|
| `padj < 0.05` (any LFC) | **1,618** |
| `padj < 0.05` & `\|log2FC\| ≥ 1` (strict) | **39** |
| Up in T4400 / Down in DMSO | 27 |
| Down in T4400 / Up in DMSO | 12 |

Top up-regulated DEGs: **Enpp2** (+1.23, padj = 1.7e-3), **Steap4** (+1.18), **Dcn** (+1.25), **Gpr68** (+1.51). Top down-regulated: **Ptn** (−1.59), **Mest** (−1.21), **Smoc2** (−1.13). PCA of VST counts cleanly separates T4400 / T4400+LIPUS from the DMSO / DMSO+LIPUS / T3976 cluster (PC1, top quartile of variance).

## 3. Pathway interpretation

**GSEA KEGG (76 pathways at `p.adjust < 0.05`)** — strongest signals all **down-regulated in T4400**:

| Pathway | NES |
|---|---|
| Spliceosome (mmu03040) | −2.49 |
| DNA replication (mmu03030) | −2.47 |
| Lysosome biogenesis (mmu04142) | **+2.44 ↑** |
| Cell cycle (mmu04110) | −2.34 |
| Oxidative phosphorylation (mmu00190) | −2.29 |
| Ribosome (mmu03010) | −2.25 |
| Prion / Huntington / Parkinson disease (mmu05020/05016/05012) | −2.1 to −2.3 |

**GSEA GO BP (629 terms)** mirrors this: **aerobic respiration, ATP synthesis coupled electron transport, oxidative phosphorylation, mitochondrial ATP synthesis, cellular respiration, respiratory electron transport chain, mitochondrial translation** are all top-NES-negative; **lysosome organization, lytic vacuole organization** are top-NES-positive.

The dominant biology is therefore a **T4400-driven suppression of mitochondrial energy metabolism and biosynthetic programmes (OXPHOS, ribosome, splicing, DNA replication, cell cycle)** with a **compensatory up-regulation of lysosomal / autophagic degradation**. The shared signal across the neurodegenerative-disease KEGG pathways (Prion / Huntington / Parkinson / Alzheimer) is consistent with shared mitochondrial / proteostasis components rather than disease-specific programs.

## 4. Limitations

1. Sample size is small (n = 4 per group, 24 total); shrunken LFCs are reliable but per-gene power is limited — many "real" biology with `|LFC| ≈ 0.5` is below the detection floor.
2. Mapping to KEGG / GO used SYMBOL keys because the dataset IDs are RIKEN-style cDNA names (not Ensembl gene IDs); 855 / 18,300 dropped unmapped genes may include part of the signal.
3. Multi-group design was collapsed to the strongest pairwise contrast (T4400 vs DMSO); an interaction model with `~ Group + LIPUS + Group:LIPUS` would be the proper next step for the LIPUS-rescue question.

## 5. Deliverables (all in `D:\Grade3\swxxx\week5\`)

| File | Purpose |
|---|---|
| `T4400_vs_DMSO_deseq2_apeglm.csv` | full DESeq2 + apeglm result table (13,923 rows) |
| `T4400_vs_DMSO_significant.csv` | 39 DEGs at `padj<0.05 & \|LFC\|>=1` |
| `T4400_vs_DMSO_dds.rds` | DESeqDataSet object (reproducibility) |
| `T4400_vs_DMSO_volcano.png` | volcano plot |
| `T4400_vs_DMSO_pca.png` | PCA of VST counts (all 24 samples) |
| `T4400_vs_DMSO_heatmap.png` | top-30 DEG heatmap (VST z-scores) |
| `T4400_vs_DMSO_gsea_GO_BP.csv` + `..._dotplot.png` | GSEA GO BP results |
| `T4400_vs_DMSO_gsea_KEGG.csv` + `..._dotplot.png` | GSEA KEGG results |
| `sessionInfo_enrich.rds` | R session info for the enrichment run |
| `week5_homework2_R_pipeline.R` | full reproducible R script (import → DE → GSEA → figures) |