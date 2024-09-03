AUC_null_pseudobulk <- readRDS("simulation/v3/predict/AUC/AUC_null_pseudobulk.rds")
AUC_real_pseudobulk <- readRDS("simulation/v3/predict/AUC/AUC_real_pseudobulk.rds")
AUC_null_tensor <- readRDS("simulation/v3/predict/AUC/AUC_null_tensor.rds")
AUC_real_tensor <- readRDS("simulation/v3/predict/AUC/AUC_real_tensor.rds")
AUC_nullWithinSamp_tensor <- readRDS("simulation/v3/predict/AUC/AUC_nullWithinSamp_tensor.rds")
AUC_nullTime <- readRDS("simulation/v3/predict/AUC/AUC_nullTime_tensor.rds")

AUC_real_pseudobulk$data <- AUC_real_tensor$data <- "Original"
AUC_null_pseudobulk$data <- AUC_null_tensor$data <- "Expression-based permuted"
AUC_nullTime$data <- "Pseudotime-based permuted"
AUC_nullWithinSamp_tensor$data <- "Sample-specific permuted"
AUC.pt <- rbind(AUC_real_pseudobulk, AUC_real_tensor, AUC_null_pseudobulk, AUC_null_tensor, AUC_nullWithinSamp_tensor, AUC_nullTime)
AUC.pt$data <- factor(AUC.pt$data, levels = c("Original", "Expression-based permuted", "Pseudotime-based permuted", "Sample-specific permuted"))

suppressPackageStartupMessages(library(ggplot2))
set.seed(2024)
ggplot(data = AUC.pt, aes(x = data, y = AUC, fill = data)) +
  # geom_boxplot(outlier.alpha = 0) +
  stat_summary(fun = median, geom = "bar") +
  geom_hline(yintercept = 0.5, alpha = 0.5, lty = "dashed") +
  geom_jitter(size = 0.4, position = position_jitter(0.1)) +
  theme_classic() + scale_fill_brewer(palette = "Dark2") +
  facet_grid(pred ~ method, scales = "free_x", space = "free_x") + theme(legend.position = "right") +
  labs(x = NULL, title = NULL, fill = "Simulated dataset") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
ggsave("simulation/v3/predict/AUC_method.pdf", width = 6.6, height = 3)
