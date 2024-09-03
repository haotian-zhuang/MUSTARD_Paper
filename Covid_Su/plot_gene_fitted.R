expr = readRDS("Covid_Su/data/saver_log2norm_sub.rds")
pseudotime = readRDS("Covid_Su/data/tActivate_pseudotime.rds") # 55,953
all.equal(colnames(expr), names(pseudotime))
meta <- readRDS('Covid_Su/data/meta_sample.rds')
cellanno = sub(':.*', '', names(pseudotime))
meta_cell <- data.frame(row.names = names(pseudotime), pseudotime = pseudotime, sample = cellanno,
                        symptom = meta[cellanno, "Symptom_binary"], time = meta[cellanno, "Time"])
samp = rownames(meta)[meta$Time == "First"]

gene_marked <- c("RPS12", "PRF1", "IFI44", "OAS2", "GZMK", "DUSP2",
                 "NKG7", "TCF7", "IRF7", "MX2", "CCL5", "JUN")

expr.marked <- expr[gene_marked, cellanno %in% samp]
pt <- meta_cell[meta_cell$sample %in% samp, ]
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
ggplot(data = pt.df, aes(x = pseudotime, y = expr, group = symptom, color = symptom)) +
  theme_classic() +
  facet_wrap(.~gene, scales = "free", nrow = 2) +
  # geom_point(alpha = 0.01, size = 0.1) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3)) +
  labs(title = NULL, x = "Pseudotime", y = "Expression", color = "Symptom") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right", strip.text = element_text(face = "italic")) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = c(0, 5e4))
ggsave("Covid_Su/figure/topgene/topgene_fitted_group.pdf", width = 10, height = 3)

ggplot(data = pt.df, aes(x = pseudotime, y = expr, group = sample, color = symptom)) +
  theme_classic() +
  facet_wrap(.~gene, scales = "free", nrow = 2) +
  # geom_point(alpha = 0.01, size = 0.1) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = FALSE, linewidth = 0.3) +
  labs(title = NULL, x = "Pseudotime", y = "Expression", color = "Symptom") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right", strip.text = element_text(face = "italic")) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = c(0, 5e4))
ggsave("Covid_Su/figure/topgene/topgene_fitted_sample.pdf", width = 10, height = 3)
