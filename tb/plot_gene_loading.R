res_ftsvd <- readRDS("tb/tb_res.rds")

gene.loading.df <- data.frame(res_ftsvd$B.hat)
gene.loading.df$Gene <- NA
gene_marked <- rownames(res_ftsvd$B.hat)[order(-abs(res_ftsvd$B.hat[, 9]))[1:20]]
gene.loading.df[rownames(gene.loading.df) %in% gene_marked, "Gene"] <- rownames(gene.loading.df)[rownames(gene.loading.df) %in% gene_marked]
gene.loading.df$Chromosome <- "Autosome"
gene.loading.df[(rownames(gene.loading.df) %in% gene_marked) & gene.loading.df$Component9 > 0, "Chromosome"] <- "X"
gene.loading.df[(rownames(gene.loading.df) %in% gene_marked) & gene.loading.df$Component9 < 0, "Chromosome"] <- "Y"

# gene_marked.pc1 <- rownames(res_ftsvd$B.hat)[order(-abs(res_ftsvd$B.hat[, 1]))[1:20]]
gene_marked.pc1 <- c("MYO1F", "ZEB2", "GZMA", "NKG7", "CCL5", "LTB", "JUNB", "RPS3A", "IL7R", "TCF7")
gene.loading.df[rownames(gene.loading.df) %in% gene_marked.pc1, "Gene"] <- rownames(gene.loading.df)[rownames(gene.loading.df) %in% gene_marked.pc1]
gene.loading.df[(rownames(gene.loading.df) %in% gene_marked.pc1) & gene.loading.df$Component1 > 0, "Chromosome"] <- "Up"
gene.loading.df[(rownames(gene.loading.df) %in% gene_marked.pc1) & gene.loading.df$Component1 < 0, "Chromosome"] <- "Down"

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(RColorBrewer))
ggplot(data = gene.loading.df, aes(x = Component1, y = Component9, label = Gene)) +
  geom_point(size = 0.1, alpha = 0.5, aes(color = Chromosome)) + theme_classic() +
  labs(x = "Component 1", y = "Component 9", title = NULL) +
  scale_color_manual(values = c("X" = brewer.pal(8, "Set1")[1], "Y" = brewer.pal(8, "Set1")[2], "Autosome" = "black",
                                "Up" = "#A55194", "Down" = "#43AA8B")) +
  geom_text_repel(seed = 2024, aes(fontface = "italic", color = Chromosome, segment.size = 0.1),
                  max.overlaps = Inf, min.segment.length = 0, size = 2.5) +
  theme(legend.position = "none")
ggsave("tb/figure/loading/gene/loading_gene_pc1&9.pdf", width = 3.5, height = 3)
