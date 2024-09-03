res_stephenson <- readRDS("Covid_Stephenson/covid_stephenson_res.rds")
res_su <- readRDS("Covid_Su/covid_su_res.rds")

gene_loading_stephenson <- res_stephenson$B.hat
gene_loading_su <- res_su$B.hat

genelist <- intersect(rownames(gene_loading_stephenson), rownames(gene_loading_su))
gene_loading_stephenson <- gene_loading_stephenson[genelist, ]
gene_loading_su <- gene_loading_su[genelist, ]

cor(gene_loading_su[, 3], gene_loading_stephenson[, 4], method = "spearman")
cor(gene_loading_su[, 3], gene_loading_stephenson[, 4])

top.overlap <- rbind(
  data.frame(numtop = seq(from = 10, to = 100, by = 10),
             propOverlap = sapply(seq(from = 10, to = 100, by = 10), function(i) {
               length(intersect(order(-abs(gene_loading_su[, 1]))[1:i], order(-abs(gene_loading_stephenson[, 1]))[1:i]))/i }),
             Su = "Component 1", Stephenson = "Component 1"),
  data.frame(numtop = seq(from = 10, to = 100, by = 10),
             propOverlap = sapply(seq(from = 10, to = 100, by = 10), function(i) {
               length(intersect(order(-abs(gene_loading_su[, 2]))[1:i], order(-abs(gene_loading_stephenson[, 1]))[1:i]))/i }),
             Su = "Component 2", Stephenson = "Component 1"),
  data.frame(numtop = seq(from = 10, to = 100, by = 10),
             propOverlap = sapply(seq(from = 10, to = 100, by = 10), function(i) {
               length(intersect(order(-abs(gene_loading_su[, 1]))[1:i], order(-abs(gene_loading_stephenson[, 4]))[1:i]))/i }),
             Su = "Component 1", Stephenson = "Component 4"),
  data.frame(numtop = seq(from = 10, to = 100, by = 10),
             propOverlap = sapply(seq(from = 10, to = 100, by = 10), function(i) {
               length(intersect(order(-abs(gene_loading_su[, 2]))[1:i], order(-abs(gene_loading_stephenson[, 4]))[1:i]))/i }),
             Su = "Component 2", Stephenson = "Component 4"),
  data.frame(numtop = seq(from = 10, to = 100, by = 10),
             propOverlap = sapply(seq(from = 10, to = 100, by = 10), function(i) {
               length(intersect(order(-abs(gene_loading_su[, 3]))[1:i], order(-abs(gene_loading_stephenson[, 4]))[1:i]))/i }),
             Su = "Component 3", Stephenson = "Component 4"),
  data.frame(numtop = seq(from = 10, to = 100, by = 10),
             propOverlap = sapply(seq(from = 10, to = 100, by = 10), function(i) {
               length(intersect( union(order(-abs(gene_loading_su[, 2]))[1:i], order(-abs(gene_loading_su[, 3]))[1:i]),
                                 order(-abs(gene_loading_stephenson[, 4]))[1:i]))/i }),
             Su = "Component 2 or 3", Stephenson = "Component 4")
)
top.overlap$Su <- factor(top.overlap$Su, levels = c("Component 1", "Component 2",
                                                    "Component 3", "Component 2 or 3"))
top.overlap$Stephenson <- factor(top.overlap$Stephenson, levels = c("Component 1", "Component 4"))

library(ggplot2)
ggplot(data = top.overlap, aes(x = numtop, y = propOverlap, color = Su)) +
  geom_point() + geom_line() +
  theme_classic() + theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(from = 0, to = 100, by = 20)) +
  scale_y_continuous(limits = c(0, 1)) + scale_color_brewer(palette = "Set1") +
  labs(x = "Number of top genes", y = "Overlap proportion", title = "COVID-Stephenson", color = "COVID-Su") +
  facet_wrap(~Stephenson)
ggsave("Covid_Stephenson/figure/topgene/topgene_overlap_pc1&4.pdf", width = 5, height = 2)
