res_ftsvd <- readRDS("Covid_Su/covid_su_res.rds")
# revise the signs of sample loadings and gene loadings Comp 3
res_ftsvd$A.hat[, 3] <- res_ftsvd$A.hat[, 3] * (-1)
res_ftsvd$B.hat[, 3] <- res_ftsvd$B.hat[, 3] * (-1)

gene.embedding <- sweep(res_ftsvd$B.hat, 2, res_ftsvd$Lambda, "*")

hl <- hclust(as.dist(1 - cor(Matrix::t(gene.embedding))), method = "average")
gene.clu <- cutree(hl, k = 5)
table(gene.clu)

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(uwot))
gene.umap <- umap(X = gene.embedding, metric = "correlation", seed = 2024)
gene.umap.df <- data.frame(gene.umap)
colnames(gene.umap.df) <- c("UMAP1", "UMAP2")
gene.umap.df$Module <- as.character(gene.clu)

gene.umap.df$Gene <- NA
gene_marked <- c("RPS12", "PRF1", "NKG7", "TCF7", "IFI44", "OAS2", "IRF7", "MX2",
                 "GZMK", "DUSP2", "CCL5", "JUN", "MT-CO1")
gene.umap.df[rownames(gene.umap.df) %in% gene_marked, "Gene"] <- rownames(gene.umap.df)[rownames(gene.umap.df) %in% gene_marked]

ggplot(data = gene.umap.df, aes(x = UMAP1, y = UMAP2, label = Gene)) +
  geom_point(size = 0.5, alpha = 1, aes(color = Module)) + scale_color_brewer(palette = "Set2") + theme_classic() +
  geom_text_repel(seed = 2024, aes(fontface = "italic", segment.size = 0.4), max.overlaps = Inf, min.segment.length = 0, size = 4)
ggsave("Covid_Su/figure/metagene/umap_gene.pdf", width = 4.5, height = 4)
