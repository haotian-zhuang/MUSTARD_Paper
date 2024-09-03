expr_sc <- readRDS("simulation/v3/noise_0.1/expr_median_scaled.rds")
samp_group <- readRDS("simulation/v3/metasamp.rds")
pseudotime <- as.numeric(sub('.*:', '', colnames(expr_sc)))
names(pseudotime) <- colnames(expr_sc)
str(expr_sc)
str(pseudotime)

cellanno = sub(':.*', '', names(pseudotime))
names(cellanno) <- names(pseudotime)
length(unique(cellanno))
samp <- names(samp_group)

source("function/scFTSVD.R")

set.seed(2023) ## not the same to 2024!!!
folds <- sample(cut(1:length(samp), breaks = 10, labels = F))
AUC_LR <- AUC_RF <- vector(length = length(unique(folds)))

start_time <- Sys.time()

for (i in 1:length(unique(folds))) {

print(sprintf("Perform the %dth validation", i))

sample_test = samp[folds == i]
sample_train = setdiff(samp, sample_test)

pseudotime_train = pseudotime[cellanno %in% sample_train]
pseudotime_test = pseudotime[cellanno %in% sample_test]

system.time( res_ftsvd <- ftsvd(expr = expr_sc, pseudotime = pseudotime_train,
                                r = 3, smooth = 1e-3) )[3]/60

sample_load_test = est_test_sample(expr_test = expr_sc, pseudotime_test = pseudotime_test, cellanno_test = NULL, res_tempted = res_ftsvd)

Atrain <- data.frame(res_ftsvd$A.hat)
Atrain$Group <- samp_group[rownames(Atrain)]
Atrain$regfact <- factor(Atrain$Group)
# LR_train <- glm(regfact ~ Component1 + Component2 + Component3,
#                  data = Atrain, family = binomial(link = "logit"))
LR_train <- nnet::multinom(regfact ~ Component1 + Component2 + Component3,
                           data = Atrain)

Atest <- data.frame(sample_load_test)
Atest$Group <- samp_group[rownames(Atest)]
# Atest$probLR <- predict(LR_train, newdata = Atest, type = "response")
prob_LR <- predict(LR_train, newdata = Atest, type = "probs")

RF_train <- randomForest::randomForest(regfact ~ Component1 + Component2 + Component3,
                                       data = Atrain)
# Atest$probRF <- predict(RF_train, newdata = Atest, type = "prob")[, 2]
prob_RF <- predict(RF_train, newdata = Atest, type = "prob")

Atest$regfact <- factor(Atest$Group)
# AUC_LR[i] = pROC::auc(pROC::roc(Atest$regfact, Atest$probLR, direction = "<"))
# AUC_RF[i] = pROC::auc(pROC::roc(Atest$regfact, Atest$probRF, direction = "<"))
AUC_LR[i] = pROC::multiclass.roc(Atest$regfact, prob_LR)$auc
AUC_RF[i] = pROC::multiclass.roc(Atest$regfact, prob_RF)$auc
}
Sys.time() - start_time

AUC.tensor <- data.frame(AUC = c(AUC_LR, AUC_RF), pred = rep(c("Logistic regression", "Random forest"), each = length(unique(folds))), method = "MUSTARD")
saveRDS(AUC.tensor, file = "simulation/v3/predict/AUC/AUC_real_tensor.rds")
