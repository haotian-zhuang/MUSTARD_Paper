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

samp_group <- readRDS("simulation/v3/metasamp.rds")
samp <- names(samp_group)

set.seed(2023)
folds <- sample(cut(1:length(samp), breaks = 10, labels = F))
AUC_LR <- AUC_RF <- vector(length = length(unique(folds)))

start_time <- Sys.time()

for (i in 1:length(unique(folds))) {
  
  print(sprintf("Perform the %dth validation", i))
  sample_test = samp[folds == i]
  sample_train = setdiff(samp, sample_test)
  
  svd_A <- Matrix::t(expr.agg.sc) ## 02/18/24
  svd_A_train <- svd_A[sample_train, ]
  svd_train <- irlba::irlba(A = svd_A_train, nv = 3)
  sample_load_train <- svd_train$u
  gene_load_train <- svd_train$v
  rownames(sample_load_train) <- rownames(svd_A_train)
  rownames(gene_load_train) <- colnames(svd_A_train)
  colnames(sample_load_train) <- colnames(gene_load_train) <- paste0('Component', 1:3)
  d_train <- svd_train$d
  D_inv_train <- diag(1 / d_train)
  rownames(D_inv_train) <- colnames(D_inv_train) <- paste0('Component', 1:3)
  svd_A_test <- svd_A[sample_test, ]
  sample_load_test = svd_A_test %*% gene_load_train %*% D_inv_train
  
  Atrain <- data.frame(sample_load_train)
  Atrain$Group <- samp_group[rownames(Atrain)]
  Atrain$regfact <- factor(Atrain$Group)
  
  # LR_train <- glm(regfact ~ Component1 + Component2 + Component3,
  #                 data = Atrain, family = binomial(link = "logit"))
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

AUC.pseudobulk <- data.frame(AUC = c(AUC_LR, AUC_RF), pred = rep(c("Logistic regression", "Random forest"), each = length(unique(folds))), method = "Pseudobulk-PCA")
saveRDS(AUC.pseudobulk, file = "simulation/v3/predict/AUC/AUC_real_pseudobulk.rds")
