res_ftsvd <- readRDS("Covid_Su/covid_su_res.rds")
# revise the signs of sample loadings and gene loadings Comp 3
res_ftsvd$A.hat[, 3] <- res_ftsvd$A.hat[, 3] * (-1)
res_ftsvd$B.hat[, 3] <- res_ftsvd$B.hat[, 3] * (-1)

# revise the signs of sample loadings and gene loadings Comp 4
res_ftsvd$A.hat[, 4] <- res_ftsvd$A.hat[, 4] * (-1)
res_ftsvd$B.hat[, 4] <- res_ftsvd$B.hat[, 4] * (-1)

suppressMessages(source("function/plot.R"))
plot_time_loading(res = res_ftsvd, linewidth = 0.8) + scale_color_manual(values = c("#00AFBB", "#DE9ED6", "#E7B800", "#FC4E07", "#98DF8A"))
ggsave("Covid_Su/figure/loading/temporal/loading_pseudotime.pdf", width = 4.5, height = 3)

meta <- readRDS('Covid_Su/data/meta_sample.rds')
plot_sample_loading(res = res_ftsvd, meta = meta, group = "Symptom_binary", i = 2, j = 3, size = 0.8) +
  labs(color = "Symptom", shape = "Symptom") +
  scale_color_brewer(palette = "Set1")
ggsave("Covid_Su/figure/loading/sample/loading_sample_pc2&3_simple.pdf", width = 4, height = 3)

suppressPackageStartupMessages(library(ggExtra))
pdf("Covid_Su/figure/loading/sample/loading_sample_pc2&3.pdf", width = 3.5, height = 3)
ggMarginal(plot_sample_loading(res = res_ftsvd, meta = meta, group = "Symptom_binary", i = 2, j = 3, size = 0.8) +
             labs(color = "Symptom", shape = "Symptom") + scale_color_brewer(palette = "Set1") +
             theme(legend.position = "none", legend.box.margin = margin(-10, -10, -10, -10),
                   plot.title = element_text(margin = margin(t = 0.1, b = 0.1, unit = "cm"))),
           type = "boxplot", groupColour = T, size = 10, outlier.shape = NA)
dev.off()

plot_gene_loading_bar(res = res_ftsvd, i = 1, numtop = 15) + scale_fill_manual(values = c("#A55194", "#43AA8B")) + labs(title = NULL)
ggsave("Covid_Su/figure/loading/gene/loading_gene_pc1_top15.pdf", width = 3, height = 2.4)
plot_gene_loading_bar(res = res_ftsvd, i = 2, numtop = 15) + scale_fill_manual(values = c("#A55194", "#43AA8B")) + labs(title = NULL)
ggsave("Covid_Su/figure/loading/gene/loading_gene_pc2_top15.pdf", width = 3, height = 2.4)
plot_gene_loading_bar(res = res_ftsvd, i = 3, numtop = 15) + scale_fill_manual(values = c("#A55194", "#43AA8B")) + labs(title = NULL)
ggsave("Covid_Su/figure/loading/gene/loading_gene_pc3_top15.pdf", width = 3, height = 2.4)
plot_gene_loading_bar(res = res_ftsvd, i = 4, numtop = 15) + scale_fill_manual(values = c("#A55194", "#43AA8B")) + labs(title = NULL)
ggsave("Covid_Su/figure/loading/gene/loading_gene_pc4_top15.pdf", width = 3, height = 2.4)
plot_gene_loading_bar(res = res_ftsvd, i = 5, numtop = 15) + scale_fill_manual(values = c("#A55194", "#43AA8B")) + labs(title = NULL)
ggsave("Covid_Su/figure/loading/gene/loading_gene_pc5_top15.pdf", width = 3, height = 2.4)

wilcox.test(res_ftsvd$A.hat[, 1] ~ meta[rownames(res_ftsvd$A.hat), ]$Symptom_binary)[["p.value"]] # 0.04988095
wilcox.test(res_ftsvd$A.hat[, 2] ~ meta[rownames(res_ftsvd$A.hat), ]$Symptom_binary)[["p.value"]] # 0.002036896
wilcox.test(res_ftsvd$A.hat[, 3] ~ meta[rownames(res_ftsvd$A.hat), ]$Symptom_binary)[["p.value"]] # 1.784042e-05
wilcox.test(res_ftsvd$A.hat[, 4] ~ meta[rownames(res_ftsvd$A.hat), ]$Symptom_binary)[["p.value"]] # 0.02218994
wilcox.test(res_ftsvd$A.hat[, 5] ~ meta[rownames(res_ftsvd$A.hat), ]$Symptom_binary)[["p.value"]] # 0.227113
