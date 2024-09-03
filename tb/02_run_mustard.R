load("tb/tb_median.Rdata")
highgene <- rownames(expr)[rowSums(expr > 0.1) > (ncol(expr) * 0.05)] # 6,127
expr <- expr[highgene, ]
range(rowSums(expr > 0.1)) # 454 9,038
sl <- matrixStats::rowSds(expr)
names(sl) <- rownames(expr)
ml <- Matrix::rowMeans(expr)

ms <- data.frame(sl = sl, ml = ml, row.names = rownames(expr))
loessfit <- loess(sl ~ ml, data = ms)
highgene <- rownames(expr)[loessfit$residuals > 0] # 1,790
library(ggplot2)
ggplot(data = cbind(ms, data.frame(fit = loessfit$fitted))) + geom_point(aes(x = ml, y = sl), size = 0.5, alpha = 0.8) +
  geom_line(aes(x = ml, y = fit), color = "blue")

expr_sc <- (expr - ml)/sl
expr_sc <- expr_sc[highgene, ]
range(expr_sc)
expr_sc[expr_sc > 10] = 10
expr_sc[expr_sc < -10] = -10
all.equal(colnames(expr_sc), names(pseudotime))

save(expr_sc, pseudotime, file = "tb/tb_median_scaled.Rdata")

load("tb/tb_median_scaled.Rdata")
cellanno <- sub(':.*', '', names(pseudotime))
length(unique(cellanno))

source('function/scFTSVD.R')
start_time <- Sys.time()
system.time( res_ftsvd <- ftsvd(expr = expr_sc, pseudotime = pseudotime, cellanno = NULL,
                                r = 10, smooth = 1e-3) )[3]/60
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken)

saveRDS(res_ftsvd, file = "tb/tb_res.rds")
