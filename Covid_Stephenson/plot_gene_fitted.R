expr = readRDS("Covid_Stephenson/data/saver_log2norm_sub.rds")
pseudotime = readRDS("Covid_Stephenson/data/tActivate_pseudotime.rds")
all.equal(colnames(expr), names(pseudotime))
meta <- readRDS("Covid_Stephenson/data/meta_sample.rds")
cellanno = readRDS("Covid_Stephenson/data/cellanno.rds")
meta_cell <- data.frame(row.names = names(pseudotime), pseudotime = pseudotime, sample = cellanno,
                        symptom = meta[cellanno, "Symptom_binary"], site = meta[cellanno, "Site"])

res_ftsvd <- readRDS("Covid_Stephenson/covid_stephenson_res.rds")
gene_marked <- c(rownames(res_ftsvd$B.hat)[order(-abs(res_ftsvd$B.hat[, 2]))[1:2]],
                 rownames(res_ftsvd$B.hat)[order(-abs(res_ftsvd$B.hat[, 3]))[1:2]])

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
ggplot(data = pt.df, aes(x = pseudotime, y = expr, group = site, color = site)) +
  theme_classic() +
  facet_wrap(.~gene, scales = "free", nrow = 1) +
  # geom_point(alpha = 0.01, size = 0.1) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3)) +
  labs(title = NULL, x = "Pseudotime", y = "Expression", color = "Site") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom", strip.text = element_text(face = "italic")) +
  scale_color_manual(values = c("Cambridge" = "#00AFBB", "Ncl" = "#E7B800", "Sanger" = "#FC4E07")) +
  scale_x_continuous(breaks = c(0, 6e4))
ggsave("Covid_Stephenson/figure/topgene/topgene_fitted_group.pdf", width = 6, height = 3)

ggplot(data = pt.df, aes(x = pseudotime, y = expr, group = sample, color = site)) +
  theme_classic() +
  facet_wrap(.~gene, scales = "free", nrow = 1) +
  # geom_point(alpha = 0.01, size = 0.1) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = FALSE, linewidth = 0.3) +
  labs(title = NULL, x = "Pseudotime", y = "Expression", color = "Site") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom", strip.text = element_text(face = "italic")) +
  scale_color_manual(values = c("Cambridge" = "#00AFBB", "Ncl" = "#E7B800", "Sanger" = "#FC4E07")) +
  scale_x_continuous(breaks = c(0, 6e4))
ggsave("Covid_Stephenson/figure/topgene/topgene_fitted_sample.pdf", width = 6, height = 3)
