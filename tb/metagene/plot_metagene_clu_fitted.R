res_ftsvd <- readRDS("tb/tb_res.rds")

gene.embedding <- sweep(res_ftsvd$B.hat, 2, res_ftsvd$Lambda, "*")

hl <- hclust(as.dist(1 - cor(Matrix::t(gene.embedding))), method = "average")
gene.clu <- cutree(hl, k = 5)
table(gene.clu)

gene.clu.weight <- data.frame(row.names = names(gene.clu))
for (i in sort(unique(gene.clu))) {
  gene.clu.weight[, paste("Module", i)] <- (gene.clu == i) * 1 / sum(gene.clu == i)
}

meta = readRDS("/hpc/group/jilab/hz/tensordata/tb/design.rds")
meta <- data.frame(Sample_name = rownames(meta), gender = ifelse(meta[, "sex"] == 1, "Female", "Male"))
load("tb/tb_median_scaled.Rdata")
all.equal(colnames(expr_sc), names(pseudotime))

cellanno = sub(':.*', '', names(pseudotime))
meta_cell <- data.frame(row.names = names(pseudotime), pseudotime = pseudotime, sample = cellanno,
                        gender = meta[cellanno, "gender"])

# original expr is too large
expr.agg <- t(gene.clu.weight) %*% expr_sc # rank (+contrast) by pseudotime
pt <- meta_cell
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
ggplot(data = pt.df, aes(x = pseudotime, y = expr, group = gender, color = gender)) +
  theme_classic() +
  facet_wrap(.~comp, scales = "free", nrow = 1) +
  # geom_point(alpha = 0.01, size = 0.1) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3)) +
  labs(title = NULL, x = "Pseudotime", y = "Expression", color = "Gender") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = c(0, 3e5), labels = expression(0, 3 %*% 10^5))
ggsave("tb/figure/metagene/clu/metagene_clu_fitted_group.pdf", width = 9.5, height = 2.5)

ggplot(data = pt.df, aes(x = pseudotime, y = expr, group = sample, color = gender)) +
  theme_classic() +
  facet_wrap(.~comp, scales = "free", nrow = 1) +
  # geom_point(alpha = 0.01, size = 0.1) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = FALSE, linewidth = 0.3) +
  labs(title = NULL, x = "Pseudotime", y = "Expression", color = "Gender") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = c(0, 3e5), labels = expression(0, 3 %*% 10^5))
ggsave("tb/figure/metagene/clu/metagene_clu_fitted_sample.pdf", width = 9.5, height = 2.5)
