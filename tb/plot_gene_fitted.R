meta = readRDS("/hpc/group/jilab/hz/tensordata/tb/design.rds")
meta <- data.frame(Sample_name = rownames(meta), gender = ifelse(meta[, "sex"] == 1, "Female", "Male"))
load("tb/tb_median.Rdata")
all.equal(colnames(expr), names(pseudotime))

cellanno = sub(':.*', '', names(pseudotime))
meta_cell <- data.frame(row.names = names(pseudotime), pseudotime = pseudotime, sample = cellanno,
                        gender = meta[cellanno, "gender"])

# original expr is too large
gene_marked <- c("MYO1F", "ZEB2", "TTTY15", "ZFY",
                 "LTB", "JUNB", "XIST", "EIF1AX")

expr.marked <- expr[gene_marked, ]
pt <- meta_cell
all.equal(rownames(pt), colnames(expr.marked))

pt.df <- data.frame()
for (i in 1:nrow(expr.marked)) {
  pt.tmp <- pt
  pt.tmp$gene <- rownames(expr.marked)[i]
  pt.tmp$expr <- expr.marked[i, ]
  if (nrow(pt.df) == 0) {
    pt.df <- pt.tmp
  } else {
    pt.df <- rbind(pt.df, pt.tmp)
  }
}

pt.df$gene <- factor(pt.df$gene, levels = rownames(expr.marked))

suppressPackageStartupMessages(library(ggplot2))
ggplot(data = pt.df, aes(x = pseudotime, y = expr, group = gender, color = gender)) +
  theme_classic() +
  facet_wrap(.~gene, scales = "free", nrow = 2) +
  # geom_point(alpha = 0.01, size = 0.1) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3)) +
  labs(title = NULL, x = "Pseudotime", y = "Expression", color = "Gender") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom", strip.text = element_text(face = "italic")) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = c(0, 3e5), labels = expression(0, 3 %*% 10^5))
ggsave("tb/figure/topgene/topgene_fitted_group.pdf", width = 6, height = 3.5)

ggplot(data = pt.df, aes(x = pseudotime, y = expr, group = sample, color = gender)) +
  theme_classic() +
  facet_wrap(.~gene, scales = "free", nrow = 2) +
  # geom_point(alpha = 0.01, size = 0.1) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = FALSE, linewidth = 0.3) +
  labs(title = NULL, x = "Pseudotime", y = "Expression", color = "Gender") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom", strip.text = element_text(face = "italic")) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = c(0, 3e5), labels = expression(0, 3 %*% 10^5))
ggsave("tb/figure/topgene/topgene_fitted_sample.pdf", width = 6, height = 3.5)
