#!/usr/bin/env Rscript
# =============================================================================
# Week 5 Homework 2 — Bulk RNA-seq DE + pathway analysis (T4400 vs DMSO)
# Reproducible pipeline: import EMP-Web session RDS -> DESeq2 + apeglm ->
# gseGO / gseKEGG -> figures.
#
# Toolchain: R 4.4.1, DESeq2 1.46.0, apeglm 1.28.0, clusterProfiler 4.14.6,
#            enrichplot 1.26.6, org.Mm.eg.db 3.20.0
#
# Run with:
#   Rscript week5_homework2_R_pipeline.R
# =============================================================================

.libPaths(c(
  "D:/Grade3/swxxx/EasyMultiProfiler-Web-main/.local_run/R_libs",
  .libPaths()
))

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(MultiAssayExperiment)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

options(timeout = 600)  # KEGG REST API

OUT <- "D:/Grade3/swxxx/week5"
SESSION_RDS <- "D:/Grade3/swxxx/EasyMultiProfiler-Web-main/.local_run/data/sessions/FGXmxOTGqJJjOe3fomMfurfl/mae.rds"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cat("===== Week 5 Homework 2 — RNA-seq DE pipeline =====\n")

# ---- 1. Load EMP session ----
mae <- readRDS(SESSION_RDS)
se <- experiments(mae)[["rnaseq_course"]]

# Group lives on MAE-level colData (not SE)
cd_mae <- as.data.frame(colData(mae))
colData(se)$Group <- cd_mae[colnames(se), "Group"]

counts <- assay(se, "counts")
cd <- as.data.frame(colData(se))
cd$Group <- factor(cd$Group)
cd$Group <- relevel(cd$Group, ref = "DMSO")
cat("Group sizes:\n"); print(table(cd$Group))

# ---- 2. DESeq2 ----
dds <- DESeqDataSetFromMatrix(countData = counts, colData = cd, design = ~ Group)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
cat("Genes after filter:", nrow(dds), "\n")

dds <- DESeq(dds)
cat("resultsNames:\n"); print(resultsNames(dds))

# Save the dds
saveRDS(dds, file.path(OUT, "T4400_vs_DMSO_dds.rds"))

# ---- 3. T4400 vs DMSO with apeglm shrinkage ----
res <- results(dds, name = "Group_T4400_vs_DMSO")
resLFC <- lfcShrink(dds, coef = "Group_T4400_vs_DMSO", type = "apeglm")

de <- as.data.frame(resLFC)
de$feature <- rownames(de)
de <- de[, c("feature", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]

de$direction <- ifelse(is.na(de$padj), NA,
                  ifelse(de$padj >= 0.05, "ns",
                    ifelse(de$log2FoldChange >= 1, "Up",
                      ifelse(de$log2FoldChange <= -1, "Down", "ns"))))

write.csv(de, file.path(OUT, "T4400_vs_DMSO_deseq2_apeglm.csv"), row.names = FALSE)
sig <- subset(de, padj < 0.05 & abs(log2FoldChange) >= 1)
write.csv(sig, file.path(OUT, "T4400_vs_DMSO_significant.csv"), row.names = FALSE)

cat("\nDE summary T4400 vs DMSO:\n")
cat("  padj<0.05 (any LFC):", sum(de$padj < 0.05, na.rm = TRUE), "\n")
cat("  padj<0.05 & |LFC|>=1:", nrow(sig), " (",
    sum(sig$log2FoldChange > 0), " Up / ",
    sum(sig$log2FoldChange < 0), " Down )\n", sep = "")

# ---- 4. Volcano ----
de_plot <- de[!is.na(de$padj), ]
de_plot$threshold <- factor(
  ifelse(de_plot$padj < 0.05 & abs(de_plot$log2FoldChange) >= 1, "Sig", "NS"),
  levels = c("NS", "Sig")
)
top_up   <- head(de_plot[de_plot$log2FoldChange > 0, ][order(de_plot[de_plot$log2FoldChange > 0, "padj"]), ], 10)
top_down <- head(de_plot[de_plot$log2FoldChange < 0, ][order(de_plot[de_plot$log2FoldChange < 0, "padj"]), ], 10)
top_lab <- rbind(top_up, top_down)
de_plot$label <- ifelse(de_plot$feature %in% top_lab$feature, de_plot$feature, "")

p_vol <- ggplot(de_plot, aes(log2FoldChange, -log10(padj), color = threshold)) +
  geom_point(alpha = 0.4, size = 0.7) +
  scale_color_manual(values = c("NS" = "grey70", "Sig" = "firebrick")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_text(aes(label = de_plot$label), size = 2.5, vjust = -0.7, show.legend = FALSE, color = "black") +
  theme_minimal(base_size = 12) +
  labs(title = "T4400 vs DMSO — Volcano plot (DESeq2 + apeglm)",
       subtitle = paste0("padj<0.05 & |log2FC|>=1, n = ", nrow(sig), " DEGs"),
       x = "log2 Fold Change (shrunken)", y = "-log10 adjusted p-value", color = NULL) +
  theme(legend.position = "top")
ggsave(file.path(OUT, "T4400_vs_DMSO_volcano.png"), p_vol,
       width = 9, height = 6.5, dpi = 160)

# ---- 5. PCA ----
vsd <- vst(dds, blind = FALSE)
pca <- plotPCA(vsd, intgroup = "Group", returnData = TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))
p_pca <- ggplot(pca, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.85) +
  stat_ellipse(level = 0.7, show.legend = FALSE) +
  scale_color_brewer(palette = "Set1") +
  theme_minimal(base_size = 12) +
  labs(title = "PCA on VST counts — 24 samples / 6 groups",
       x = paste0("PC1 (", percentVar[1], "% var)"),
       y = paste0("PC2 (", percentVar[2], "% var)"))
ggsave(file.path(OUT, "T4400_vs_DMSO_pca.png"), p_pca,
       width = 8.5, height = 6, dpi = 160)

# ---- 6. Heatmap of top DEGs ----
top30 <- head(sig[order(-abs(sig$log2FoldChange)), ], 30)
vst_mat <- assay(vsd)[top30$feature, ]
vst_z <- (vst_mat - rowMeans(vst_mat)) / apply(vst_mat, 1, sd)
keep <- colData(vsd)$Group %in% c("DMSO", "T4400")
vst_z <- vst_z[, keep]
ann <- data.frame(Group = colData(vsd)$Group[keep])
rownames(ann) <- colnames(vst_z)

png(file.path(OUT, "T4400_vs_DMSO_heatmap.png"),
    width = 1100, height = 1400, res = 150)
pheatmap(vst_z,
         annotation_col = ann,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         breaks = seq(-2.5, 2.5, length.out = 101),
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE, fontsize_row = 6,
         main = "Top 30 DEGs T4400 vs DMSO (VST z-scores)")
dev.off()

# ---- 7. Helper: SYMBOL -> ENTREZ ----
sym2ez <- function(symbols) {
  m <- AnnotationDbi::select(org.Mm.eg.db,
                             keys = unique(symbols),
                             keytype = "SYMBOL",
                             columns = "ENTREZID")
  m <- m[!is.na(m$ENTREZID), ]
  m <- m[!duplicated(m$SYMBOL), ]
  m
}

# ---- 8. GSEA GO BP ----
de_rank <- de[!is.na(de$pvalue) & de$pvalue > 0, ]
de_rank$rank_metric <- sign(de_rank$log2FoldChange) * -log10(de_rank$pvalue)
de_rank <- de_rank[order(-de_rank$rank_metric), ]

m_all <- sym2ez(de_rank$feature)
de_rank <- merge(de_rank, m_all, by.x = "feature", by.y = "SYMBOL", all.x = FALSE)
gsea_rank <- de_rank$rank_metric
names(gsea_rank) <- de_rank$ENTREZID
gsea_rank <- gsea_rank[!duplicated(names(gsea_rank))]
gsea_rank <- sort(gsea_rank, decreasing = TRUE)

gsea_go <- gseGO(gsea_rank, OrgDb = org.Mm.eg.db, keyType = "ENTREZID",
                 ont = "BP", minGSSize = 15, maxGSSize = 500,
                 pvalueCutoff = 0.05, seed = TRUE, verbose = FALSE)
write.csv(as.data.frame(gsea_go),
          file.path(OUT, "T4400_vs_DMSO_gsea_GO_BP.csv"), row.names = FALSE)

png(file.path(OUT, "T4400_vs_DMSO_gsea_GO_BP_dotplot.png"),
    width = 1800, height = 1800, res = 150)
print(dotplot(gsea_go, showCategory = 20,
              title = "T4400 vs DMSO — GO BP (GSEA)") +
        theme(axis.text.y = element_text(size = 9)))
dev.off()

# ---- 9. GSEA KEGG ----
gsea_kegg <- gseKEGG(gsea_rank, organism = "mmu",
                     minGSSize = 15, maxGSSize = 500,
                     pvalueCutoff = 0.05, seed = TRUE, verbose = FALSE)
write.csv(as.data.frame(gsea_kegg),
          file.path(OUT, "T4400_vs_DMSO_gsea_KEGG.csv"), row.names = FALSE)

png(file.path(OUT, "T4400_vs_DMSO_gsea_KEGG_dotplot.png"),
    width = 1800, height = 1400, res = 150)
print(dotplot(gsea_kegg, showCategory = 20,
              title = "T4400 vs DMSO — KEGG (GSEA)") +
        theme(axis.text.y = element_text(size = 9)))
dev.off()

# ---- 10. Save session info ----
saveRDS(sessionInfo(), file.path(OUT, "sessionInfo_enrich.rds"))

cat("\n===== Done. All artefacts in", OUT, "=====\n")