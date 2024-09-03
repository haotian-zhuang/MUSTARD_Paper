res_ftsvd <- readRDS("simulation/v3/noise_0.1/sim_res.rds")
res_ftsvd$r.square
samp_group <- readRDS("simulation/v3/metasamp.rds")

# revise the signs of sample loadings and temporal loadings Comp 3
res_ftsvd$A.hat[, 3] <- res_ftsvd$A.hat[, 3] * (-1)
res_ftsvd$Phi.hat[, 3] <- res_ftsvd$Phi.hat[, 3] * (-1)

library(ggplot2)
library(RColorBrewer)
# temporal loading----
Phi.data <- res_ftsvd$Phi.hat[, 1:3]
Phi.data <- data.frame(time = res_ftsvd$time, value = as.vector(Phi.data), 
                       component = as.factor(as.vector(t(matrix(rep(1:3, nrow(Phi.data)), 3,)))))
ggplot(data = Phi.data, aes(x = time, y = value, color = component)) +
  geom_line(linewidth = 1.5) + theme_classic() +
  scale_color_manual(values = c(brewer.pal(11, "PiYG")[1], brewer.pal(11, "PuOr")[3], brewer.pal(11, "BrBG")[10])) +
  labs(x = "Pseudotime", y = "Temporal loading", color = "Component") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom")
ggsave("simulation/v3/noise_0.1/temporal_loading.pdf", width = 4, height = 4)

# sample loading----
sample.loading <- reshape2::melt(res_ftsvd$A.hat[, 1:3])
colnames(sample.loading) <- c("Sample", "Comp", "Loading")
sample.loading$Comp <- gsub("(Component)(\\d+)", "\\1 \\2", sample.loading$Comp)
sample.loading$Group <- samp_group[sample.loading$Sample]

sample.loading <- sample.loading[order(sample.loading$Group, sample.loading$Sample, sample.loading$Comp), ]
sample.loading$Sample <- factor(sample.loading$Sample, levels = unique(sample.loading$Sample))

ggplot(sample.loading, aes(x = Sample, y = Loading, fill = Group)) +
  geom_bar(stat = "identity") + labs(x = "Sample", y = paste("Sample loading"), fill = "Sample") + theme_classic() +
  scale_fill_manual(values = c(brewer.pal(9, "Set1")[3], brewer.pal(11, "Spectral")[c(11, 10)])) +
  facet_grid(Comp ~ ., scales = "free") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom")
ggsave("simulation/v3/noise_0.1/sample_loading.pdf", width = 4, height = 4)

# gene loading----
gene.loading <- reshape2::melt(res_ftsvd$B.hat[, 1:3])
colnames(gene.loading) <- c("Gene", "Comp", "Loading")
gene.loading$Comp <- gsub("(Component)(\\d+)", "\\1 \\2", gene.loading$Comp)
gene.loading$Set <- rep(rep(c("Set 1", "Set 2", "Set 3"), times = c(10, 10, 5)), 3)

ggplot(gene.loading, aes(x = Gene, y = Loading, fill = Set)) +
  geom_bar(stat = "identity") + labs(x = "Gene", y = paste("Gene loading"), fill = "Gene") + theme_classic() +
  scale_fill_manual(values = brewer.pal(8, "Set2")[1:3]) +
  facet_grid(Comp ~ ., scales = "free") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom")
ggsave("simulation/v3/noise_0.1/gene_loading.pdf", width = 4, height = 4)
