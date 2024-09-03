expr <- readRDS("Covid_Stephenson/data/saver_log2norm_sub.rds")
str(expr)
pseudotime <- readRDS("Covid_Stephenson/data/tActivate_pseudotime.rds")
str(pseudotime)
cellanno = readRDS("Covid_Stephenson/data/cellanno.rds")
length(unique(cellanno))

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
  print(paste("Sample:", i))
  expr.sample <- expr[, (cellanno == i)]
  pseudotime.sample <- pseudotime[cellanno == i]
  expr.sample.agg <- agg_median(expr = expr.sample, pseudotime = pseudotime.sample)
  colnames(expr.sample.agg) <- paste(i, colnames(expr.sample.agg), sep = ':')
  datlist[[i]] <- expr.sample.agg
  print(paste0(which(unique(cellanno) == i), '/', length(unique(cellanno))))
}

expr.agg <- do.call(cbind, datlist)
str(expr.agg) # 23,339 4,614
pseudotime.agg <- as.numeric(sub('.*:', '', colnames(expr.agg)))
names(pseudotime.agg) <- colnames(expr.agg)
str(pseudotime.agg)
length(unique(sub(':.*', '', names(pseudotime.agg))))

expr <- expr.agg
pseudotime <- pseudotime.agg
save(expr, pseudotime, file = "Covid_Stephenson/covid_stephenson_median.Rdata")

# scaled
load("Covid_Stephenson/covid_stephenson_median.Rdata")
highgene <- rownames(expr)[rowSums(expr > 0.1) > (ncol(expr) * 0.05)] # 5,981
expr <- expr[highgene, ]
range(rowSums(expr > 0.1)) # 231 4,614
sl <- matrixStats::rowSds(expr)
names(sl) <- rownames(expr)
ml <- Matrix::rowMeans(expr)

ms <- data.frame(sl = sl, ml = ml, row.names = rownames(expr))
loessfit <- loess(sl ~ ml, data = ms)
highgene <- rownames(expr)[loessfit$residuals > 0] # 1,833
library(ggplot2)
ggplot(data = cbind(ms, data.frame(fit = loessfit$fitted))) + geom_point(aes(x = ml, y = sl), size = 0.5, alpha = 0.8) +
  geom_line(aes(x = ml, y = fit), color = "blue")

expr_sc <- (expr - ml)/sl
expr_sc <- expr_sc[highgene, ]
range(expr_sc)
expr_sc[expr_sc > 10] = 10
expr_sc[expr_sc < -10] = -10
all.equal(colnames(expr_sc), names(pseudotime))

save(expr_sc, pseudotime, file = "Covid_Stephenson/covid_stephenson_median_scaled.Rdata")
