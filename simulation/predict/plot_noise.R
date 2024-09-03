noise_seq <- c(2, 3, 4, 5)
AUC.list <- vector(mode = 'list', length = length(noise_seq))
names(AUC.list) <- paste0("noise_", noise_seq)

for (noise in noise_seq) {
  AUC <- readRDS(paste0("simulation/v3/predict/AUC/AUC_real_tensor_noise_", noise, ".rds"))
  AUC$noise <- noise
  AUC.list[[paste0("noise_", noise)]] <- AUC
}
AUC <- do.call(rbind, AUC.list)

library(ggplot2)
ggplot(data = AUC, aes(x = noise, y = AUC, group = noise), fill = "white") +
  geom_boxplot(outlier.alpha = 0) +
  scale_x_continuous(breaks = noise_seq) +
  geom_hline(yintercept = 0.5, alpha = 0.5, lty = "dashed") +
  geom_jitter(size = 0.4, position = position_jitter(0.1)) +
  theme_classic() +
  facet_wrap(~pred) +
  theme(legend.position = "none") + labs(x = "Noise standard deviation", title = NULL)
ggsave("simulation/v3/predict/AUC_noise.pdf", width = 4, height = 2)
