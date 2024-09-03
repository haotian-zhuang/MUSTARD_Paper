mat <- readRDS("simulation/v3/multi_fitted.rds")
noise = 0.1
set.seed(2024)
mat <- mat + matrix(rnorm(length(mat), mean = 0, sd = noise), nrow = nrow(mat), ncol = ncol(mat))

pseudotime <- readRDS("Covid_Su/data/tActivate_pseudotime.rds")
cellanno = sub(':.*', '', names(pseudotime))
expr <- mat[, names(pseudotime)]

expr.agg <- matrix(0, nrow(expr), length(unique(cellanno)))
rownames(expr.agg) <- rownames(expr)
colnames(expr.agg) <- unique(cellanno)

for (i in unique(cellanno)) {
  expr.samp <- expr[, cellanno == i]
  expr.agg[, i] <- rowMeans(expr.samp)
  # expr.agg[, i] <- apply(expr.samp, 1, median)
}
expr.agg.sc <- (expr.agg - Matrix::rowMeans(expr.agg)) / matrixStats::rowSds(expr.agg)
range(expr.agg.sc)

# SVD----
svd_A <- Matrix::t(expr.agg.sc)
res_svd <- irlba::irlba(A = svd_A, nv = 3)
sample_load <- res_svd$u
gene_load <- res_svd$v

rownames(sample_load) <- rownames(svd_A)
rownames(gene_load) <- colnames(svd_A)
colnames(sample_load) <- colnames(gene_load) <- paste("Component", 1:3)

samp_group <- readRDS("simulation/v3/metasamp.rds")

# sample loading----
sample.loading <- reshape2::melt(sample_load)
colnames(sample.loading) <- c("Sample", "Comp", "Loading")
sample.loading$Group <- samp_group[sample.loading$Sample]

sample.loading <- sample.loading[order(sample.loading$Group, sample.loading$Sample, sample.loading$Comp), ]
sample.loading$Sample <- factor(sample.loading$Sample, levels = unique(sample.loading$Sample))

ggplot(sample.loading, aes(x = Sample, y = Loading, fill = Group)) +
  geom_bar(stat = "identity") + labs(x = "Sample", y = paste("Sample loading"), fill = "Sample") + theme_classic() +
  scale_fill_manual(values = c(brewer.pal(9, "Set1")[3], brewer.pal(11, "Spectral")[c(11, 10)])) +
  facet_grid(Comp ~ ., scales = "free") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom")
ggsave("simulation/v3/noise_0.1/sample_loading_pseudobulk.pdf", width = 4, height = 4)

# gene loading----
gene.loading <- reshape2::melt(gene_load)
colnames(gene.loading) <- c("Gene", "Comp", "Loading")
gene.loading$Set <- rep(rep(c("Set 1", "Set 2", "Set 3"), times = c(10, 10, 5)), 3)

ggplot(gene.loading, aes(x = Gene, y = Loading, fill = Set)) +
  geom_bar(stat = "identity") + labs(x = "Gene", y = paste("Gene loading"), fill = "Gene") + theme_classic() +
  scale_fill_manual(values = brewer.pal(8, "Set2")[1:3]) +
  facet_grid(Comp ~ ., scales = "free") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom")
ggsave("simulation/v3/noise_0.1/gene_loading_pseudobulk.pdf", width = 4, height = 4)
