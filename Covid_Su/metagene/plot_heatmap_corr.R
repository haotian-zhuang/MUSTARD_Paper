res_ftsvd <- readRDS("Covid_Su/covid_su_res.rds")
# revise the signs of sample loadings and gene loadings Comp 3
res_ftsvd$A.hat[, 3] <- res_ftsvd$A.hat[, 3] * (-1)
res_ftsvd$B.hat[, 3] <- res_ftsvd$B.hat[, 3] * (-1)

gene.embedding <- sweep(res_ftsvd$B.hat, 2, res_ftsvd$Lambda, "*")
gene.cor = cor(Matrix::t(gene.embedding))

hl <- hclust(as.dist(1 - gene.cor), method = "average")
gene.clu <- cutree(hl, k = 5)
table(gene.clu)

pl <- gene.cor[hl$order, hl$order]

suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(RColorBrewer))

gene_marked <- c("RPS12", "PRF1", "NKG7", "TCF7", "IFI44", "OAS2", "IRF7", "MX2",
                 "GZMK", "DUSP2", "CCL5", "JUN", "MT-CO1")

Gene.cluster <- factor(gene.clu[hl$order])
Gene.clv <- brewer.pal(8, "Set2")[1:length(unique(Gene.cluster))]
# set.seed(2024)
# Gene.clv <- sample(rainbow(length(unique(Gene.cluster)), alpha = 0.5))
names(Gene.clv) <- levels(Gene.cluster)

pdf(file = "Covid_Su/figure/metagene/heatmap_module.pdf", width = 9.2, height = 7.5)
ht <- Heatmap(matrix = pl, name = "Pearson correlation",
              col = colorRamp2(c(quantile(pl, 0.05), 0, quantile(pl, 0.95)), brewer.pal(7, "RdYlBu")[c(7, 5, 1)]),
              # clustering_distance_rows = function(m) {as.dist(1 - m)}, clustering_distance_columns = function(m) {as.dist(1 - m)},
              # clustering_method_rows = "average", column_title = "linkage = average",
              # clustering_method_columns = "average", show_column_dend = F,
              cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F,
              top_annotation = columnAnnotation(`Component 1` = res_ftsvd$B.hat[hl$order, 1],
                                                `Component 2` = res_ftsvd$B.hat[hl$order, 2],
                                                `Component 3` = res_ftsvd$B.hat[hl$order, 3],
                                                `Component 4` = res_ftsvd$B.hat[hl$order, 4],
                                                `Component 5` = res_ftsvd$B.hat[hl$order, 5],
                                                col = list(`Component 1` = colorRamp2(c(quantile(res_ftsvd$B.hat[, 1], 0.05), 0, quantile(res_ftsvd$B.hat[, 1], 0.95)), c("green", "white", "purple")),
                                                           `Component 2` = colorRamp2(c(quantile(res_ftsvd$B.hat[, 2], 0.05), 0, quantile(res_ftsvd$B.hat[, 2], 0.95)), c("green", "white", "purple")),
                                                           `Component 3` = colorRamp2(c(quantile(res_ftsvd$B.hat[, 3], 0.05), 0, quantile(res_ftsvd$B.hat[, 3], 0.95)), c("green", "white", "purple")),
                                                           `Component 4` = colorRamp2(c(quantile(res_ftsvd$B.hat[, 4], 0.05), 0, quantile(res_ftsvd$B.hat[, 4], 0.95)), c("green", "white", "purple")),
                                                           `Component 5` = colorRamp2(c(quantile(res_ftsvd$B.hat[, 5], 0.05), 0, quantile(res_ftsvd$B.hat[, 5], 0.95)), c("green", "white", "purple"))),
                                                annotation_name_side = "right", show_legend = F),
              left_annotation = rowAnnotation(Gene = anno_mark(at = which(rownames(pl) %in% gene_marked),
                                                               labels = rownames(pl)[rownames(pl) %in% gene_marked],
                                                               labels_gp = gpar(fontface = "italic", fontsize = 10),
                                                               padding = unit(1.5, "mm"), side = "left"),
                                              Module = Gene.cluster, col = list(Module = Gene.clv), show_annotation_name = list(Module = F)),
              row_split = Gene.cluster, column_split = Gene.cluster, row_title = NULL, column_title = NULL,
              use_raster = T, heatmap_legend_param = list(legend_direction = "vertical"))
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = T,
     annotation_legend_list = list(Legend(col_fun = colorRamp2(c(-1, 0, 1), c("green", "white", "purple")), title = "Gene loading",
                                          at = c(-1, 0, 1), labels = c("Negative", "Zero", "Positive"))))
dev.off()
