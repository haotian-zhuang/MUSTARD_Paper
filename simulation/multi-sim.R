expr <- readRDS('Covid_Su/data/saver_log2norm_sub.rds')
str(expr)
pseudotime <- readRDS('Covid_Su/data/tActivate_pseudotime.rds')
str(pseudotime)
cellanno = sub(':.*', '', names(pseudotime))
names(cellanno) <- names(pseudotime)
length(unique(cellanno))
all.equal(colnames(expr), names(pseudotime))

res_ftsvd <- readRDS("Covid_Su/covid_su_res.rds")
upgene <- rownames(res_ftsvd$B.hat)[order(-(res_ftsvd$B.hat[, 1]))[1:10]]
downgene <- rownames(res_ftsvd$B.hat)[order((res_ftsvd$B.hat[, 1]))[1:10]]
peakgene <- rownames(res_ftsvd$B.hat)[order((res_ftsvd$B.hat[, 5]))[1:5]]

expr.sim <- matrix(NA, nrow = length(c(upgene, downgene, peakgene)), ncol = ncol(expr))
rownames(expr.sim) <- c(upgene, downgene, peakgene)
for (i in c(upgene, downgene, peakgene)) {
  fit.model <- mgcv::gam(expr[i, ] ~ s(pseudotime, k = 3))
  expr.sim[i, ] <- predict(fit.model)
}
colnames(expr.sim) <- colnames(expr)
saveRDS(expr.sim, file = "simulation/v3/fitted.rds")

# divide the 161 samples into three groups
set.seed(2024)
samp_group <- sample(rep(c("Group 1", "Group 2", "Group 3"), times = c(54, 54, 53)))
names(samp_group) <- unique(cellanno)
head(samp_group)
saveRDS(samp_group, file = "simulation/v3/metasamp.rds")

cell_group <- samp_group[cellanno]
names(cell_group) <- names(cellanno)
head(cell_group)

# produce patterns
expr.g1 <- expr.sim[, cell_group == "Group 1"]
expr.g2 <- expr.sim[, cell_group == "Group 2"]
expr.g3 <- expr.sim[, cell_group == "Group 3"]

permute_within_sample <- function(expr) {
  cellanno = sub(':.*', '', colnames(expr))
  datlist <- vector(mode = 'list', length = length(unique(cellanno)))
  names(datlist) <- unique(cellanno)
  for (i in unique(cellanno)) {
    expr.sample <- expr[, (cellanno == i)]
    
    # expr.sample.permuted <- expr.sample[, sample(colnames(expr.sample))]
    # colnames(expr.sample.permuted) <- colnames(expr.sample)
    
    expr.sample.permuted <- Matrix::t(apply(expr.sample, 1, sample))
    colnames(expr.sample.permuted) <- colnames(expr.sample)
    
    datlist[[i]] <- expr.sample.permuted
  }
  expr.permuted <- do.call(cbind, datlist)
  expr.permuted <- expr.permuted[, colnames(expr)]
}

set.seed(202401)
expr.g2[11:20, ] <- permute_within_sample(expr = expr.g2[11:20, ])

set.seed(202402)
expr.g2[21:25, ] <- permute_within_sample(expr = expr.g2[21:25, ])

set.seed(202403)
expr.g3[21:25, ] <- permute_within_sample(expr = expr.g3[21:25, ])
mat <- cbind(expr.g1, expr.g2, expr.g3)

saveRDS(mat, file = "simulation/v3/multi_fitted.rds")
