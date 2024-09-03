res_ftsvd <- readRDS("Covid_Su/covid_su_res.rds")
# revise the signs of sample loadings and gene loadings Comp 3
res_ftsvd$A.hat[, 3] <- res_ftsvd$A.hat[, 3] * (-1)
res_ftsvd$B.hat[, 3] <- res_ftsvd$B.hat[, 3] * (-1)

gene.embedding <- sweep(res_ftsvd$B.hat, 2, res_ftsvd$Lambda, "*")

hl <- hclust(as.dist(1 - cor(Matrix::t(gene.embedding))), method = "average")
gene.clu <- cutree(hl, k = 5)
table(gene.clu)

gene.clu.weight <- data.frame(row.names = names(gene.clu))
for (i in sort(unique(gene.clu))) {
  gene.clu.weight[, paste("Module", i)] <- (gene.clu == i) * 1 / sum(gene.clu == i)
}

expr = readRDS("Covid_Su/data/saver_log2norm_sub.rds")
pseudotime = readRDS("Covid_Su/data/tActivate_pseudotime.rds") # 55,953
all.equal(colnames(expr), names(pseudotime))
meta <- readRDS('Covid_Su/data/meta_sample.rds')
cellanno = sub(':.*', '', names(pseudotime))
samp = rownames(res_ftsvd$A.hat)
meta_cell <- data.frame(row.names = names(pseudotime), pseudotime = pseudotime, sample = cellanno,
                        symptom = meta[cellanno, "Symptom_binary"], time = meta[cellanno, "Time"])

expr = expr[rownames(res_ftsvd$B.hat), ]
sl <- matrixStats::rowSds(expr)
names(sl) <- rownames(expr)
ml <- Matrix::rowMeans(expr)
expr_sc <- (expr - ml)/sl
expr_sc[expr_sc > 10] = 10
expr_sc[expr_sc < -10] = -10
str(expr_sc)

expr.agg <- t(gene.clu.weight) %*% expr_sc[, cellanno %in% samp] # rank (+contrast) by pseudotime
pt <- meta_cell[meta_cell$sample %in% samp, ]
all.equal(rownames(pt), colnames(expr.agg))

pt.df <- data.frame()
for (i in 1:nrow(expr.agg)) {
  pt.tmp <- pt
  pt.tmp$comp <- rownames(expr.agg)[i]
  pt.tmp$expr <- expr.agg[i, ]
  if (nrow(pt.df) == 0) {
    pt.df <- pt.tmp
  } else {
    pt.df <- rbind(pt.df, pt.tmp)
  }
}

suppressPackageStartupMessages(library(ggplot2))
ggplot(data = pt.df, aes(x = pseudotime, y = expr, group = symptom, color = symptom)) +
  theme_classic() +
  facet_wrap(.~comp, scales = "free", nrow = 1) +
  # geom_point(alpha = 0.01, size = 0.1) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3)) +
  labs(title = NULL, x = "Pseudotime", y = "Expression", color = "Symptom") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right") +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = c(0, 5e4))
ggsave("Covid_Su/figure/metagene/clu/metagene_clu_fitted_group.pdf", width = 10, height = 2)

ggplot(data = pt.df, aes(x = pseudotime, y = expr, group = sample, color = symptom)) +
  theme_classic() +
  facet_wrap(.~comp, scales = "free", nrow = 1) +
  # geom_point(alpha = 0.01, size = 0.1) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = FALSE, linewidth = 0.3) +
  labs(title = NULL, x = "Pseudotime", y = "Expression", color = "Symptom") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right") +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = c(0, 5e4))
ggsave("Covid_Su/figure/metagene/clu/metagene_clu_fitted_sample.pdf", width = 10, height = 2)
