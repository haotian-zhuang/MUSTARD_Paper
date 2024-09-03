res_ftsvd <- readRDS("tb/tb_res.rds")
suppressMessages(source("function/plot.R"))

plot_time_loading(res = res_ftsvd, r = c(1, 9), linewidth = 0.8) +
  scale_x_continuous(breaks = c(0, 1e5, 2e5, 3e5), labels = expression(0, 10^5, 2 %*% 10^5, 3 %*% 10^5)) +
  scale_color_manual(values = c("#F8961E", "#9C9EDE"))
ggsave("tb/figure/loading/temporal/loading_pseudotime.pdf", width = 3.5, height = 2)

meta = readRDS("/hpc/group/jilab/hz/tensordata/tb/design.rds")
meta <- data.frame(Sample_name = rownames(meta), gender = ifelse(meta[, "sex"] == 1, "Female", "Male"))

plot_sample_loading(res = res_ftsvd, meta = meta, group = "gender", i = 1, j = 9, size = 0.8) +
  labs(color = "Gender", shape = "Gender") +
  scale_color_brewer(palette = "Set1")
ggsave("tb/figure/loading/sample/loading_sample_pc1&9.pdf", width = 3.5, height = 2)

plot_gene_loading_bar(res = res_ftsvd, i = 1, numtop = 20) + scale_fill_manual(values = c("#A55194", "#43AA8B")) + labs(title = NULL)
ggsave("tb/figure/loading/gene/loading_gene_pc1_top20.pdf", width = 3, height = 3)
plot_gene_loading_bar(res = res_ftsvd, i = 9, numtop = 20) + scale_fill_brewer(palette = "Set1") + labs(title = NULL)
ggsave("tb/figure/loading/gene/loading_gene_pc9_top20.pdf", width = 3, height = 3)

load("tb/tb_median_scaled.Rdata")
plot_aggregate_gene(res = res_ftsvd, expr_sc = expr_sc, pseudotime = pseudotime, cellanno = NULL, r = c(1, 9),
                    pct = 0.2, group = "gender", bws = 300, nrow = 1) +
  scale_x_continuous(breaks = c(0, 1e5, 2e5, 3e5), labels = expression(0, 10^5, 2 %*% 10^5, 3 %*% 10^5)) +
  scale_color_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1")
ggsave("tb/figure/aggregated_gene.pdf", width = 5.5, height = 2)
