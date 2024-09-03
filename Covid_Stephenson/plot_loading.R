res_ftsvd <- readRDS("Covid_Stephenson/covid_stephenson_res.rds")
# revise the signs of sample loadings and gene loadings Comp 4
res_ftsvd$A.hat[, 4] <- res_ftsvd$A.hat[, 4] * (-1)
res_ftsvd$B.hat[, 4] <- res_ftsvd$B.hat[, 4] * (-1)

suppressMessages(source("function/plot.R"))
plot_time_loading(res = res_ftsvd, r = 1:5, linewidth = 0.8) + scale_color_manual(values = c("#00AFBB", "#DE9ED6", "#E7B800", "#FC4E07", "#98DF8A"))
ggsave("Covid_Stephenson/figure/loading/temporal/loading_pseudotime.pdf", width = 4, height = 2)

meta <- readRDS("Covid_Stephenson/data/meta_sample.rds")
plot_sample_loading(res = res_ftsvd, meta = meta, group = "Site", i = 2, j = 3, size = 0.8) +
  scale_color_manual(values = c("Cambridge" = "#00AFBB", "Ncl" = "#E7B800", "Sanger" = "#FC4E07"))
ggsave("Covid_Stephenson/figure/loading/sample/loading_sample_site_pc2&3.pdf", width = 4, height = 2)

wilcox.test(res_ftsvd$A.hat[, 4] ~ meta[rownames(res_ftsvd$A.hat), ]$Symptom_binary)[["p.value"]] # 2.13912e-05
plot_sample_loading(res = res_ftsvd, meta = meta, group = "Symptom_binary", i = 1, j = 4, size = 0.8) +
  labs(color = "Symptom", shape = "Symptom") +
  scale_color_brewer(palette = "Set1")
ggsave("Covid_Stephenson/figure/loading/sample/loading_sample_symptom_pc1&4.pdf", width = 4, height = 2)

# suppressPackageStartupMessages(library(ggExtra))
# pdf("Covid_Stephenson/figure/loading/sample/loading_sample_symptom_pc1&4.pdf", width = 4, height = 2)
# ggMarginal(plot_sample_loading(res = res_ftsvd, meta = meta, group = "Symptom_binary", i = 1, j = 4, size = 0.8) +
#              labs(color = "Symptom", shape = "Symptom") + scale_color_brewer(palette = "Set1") +
#              theme(legend.position = "bottom", legend.box.margin = margin(-10, -10, -10, -10)),
#            type = "boxplot", groupColour = T, size = 10, outlier.shape = NA)
# dev.off()

plot_gene_loading_bar(res = res_ftsvd, i = 2, numtop = 5) + scale_fill_manual(values = c("#A55194", "#43AA8B")) + labs(title = NULL)
ggsave("Covid_Stephenson/figure/loading/gene/loading_gene_pc2_top5.pdf", width = 3, height = 1.2)
plot_gene_loading_bar(res = res_ftsvd, i = 3, numtop = 5) + scale_fill_manual(values = c("#A55194", "#43AA8B")) + labs(title = NULL)
ggsave("Covid_Stephenson/figure/loading/gene/loading_gene_pc3_top5.pdf", width = 3, height = 1.2)

plot_gene_loading_bar(res = res_ftsvd, i = 1, numtop = 20) + scale_fill_manual(values = c("#A55194", "#43AA8B")) + labs(title = NULL)
ggsave("Covid_Stephenson/figure/loading/gene/loading_gene_pc1_top20.pdf", width = 3, height = 3)
plot_gene_loading_bar(res = res_ftsvd, i = 4, numtop = 20) + scale_fill_manual(values = c("#A55194", "#43AA8B")) + labs(title = NULL)
ggsave("Covid_Stephenson/figure/loading/gene/loading_gene_pc4_top20.pdf", width = 3, height = 3)
