mat <- readRDS("simulation/v3/multi_fitted.rds")
noise = 0.1
set.seed(2024)
mat <- mat + matrix(rnorm(length(mat), mean = 0, sd = noise), nrow = nrow(mat), ncol = ncol(mat))

pseudotime <- readRDS("Covid_Su/data/tActivate_pseudotime.rds")
cellanno = sub(':.*', '', names(pseudotime))
expr <- mat[, names(pseudotime)]

interval_len = 50
interval_start = seq(from = min(pseudotime), to = max(pseudotime), length.out = interval_len+1)
interval_end = interval_start[-1]
interval_start = interval_start[-(length(interval_start))]
interval = cbind(interval_start, interval_end)
colnames(interval) = c('start', 'end')

agg_median <- function(expr, pseudotime){
  expr_agg <- apply(interval, 1, function(i){
    cols <- (pseudotime >= i['start'] & pseudotime < i['end'])
    if(sum(cols) > 0){
      expr_fall <- expr[, cols, drop = F]
      expr_median <- apply(expr_fall, 1, median)
    } else {
      expr_median <- rep(NA, nrow(expr))
    }
  })
  rownames(expr_agg) <- rownames(expr)
  colnames(expr_agg) <- apply(interval, 1, mean)
  expr_agg <- expr_agg[, colSums(is.na(expr_agg)) == 0]
  expr_agg
}

datlist <- vector(mode = 'list', length = length(unique(cellanno)))
names(datlist) <- unique(cellanno)
for (i in unique(cellanno)) {
  expr.sample <- expr[, (cellanno == i)]
  pseudotime.sample <- pseudotime[cellanno == i]
  expr.sample.agg <- agg_median(expr = expr.sample, pseudotime = pseudotime.sample)
  colnames(expr.sample.agg) <- paste(i, colnames(expr.sample.agg), sep = ':')
  datlist[[i]] <- expr.sample.agg
  print(paste0(which(unique(cellanno) == i), '/', length(unique(cellanno))))
}

expr.agg <- do.call(cbind, datlist)
str(expr.agg) # 25 * 6,824
pseudotime.agg <- as.numeric(sub('.*:', '', colnames(expr.agg)))
names(pseudotime.agg) <- colnames(expr.agg)
str(pseudotime.agg)
length(unique(sub(':.*', '', names(pseudotime.agg))))

expr_sc <- (expr.agg - Matrix::rowMeans(expr.agg)) / matrixStats::rowSds(expr.agg)
range(expr_sc)
saveRDS(expr_sc, file = "simulation/v3/noise_0.1/expr_median_scaled.rds")

source("function/scFTSVD.R")
system.time( res_ftsvd <- ftsvd(expr = expr_sc, pseudotime = pseudotime.agg, r = 3, smooth = 1e-3) )[3]/60
saveRDS(res_ftsvd, file = "simulation/v3/noise_0.1/sim_res.rds")
