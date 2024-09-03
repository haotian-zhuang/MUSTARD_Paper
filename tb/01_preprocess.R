library(rhdf5)
expr.path <- "/hpc/group/jilab/hz/tensordata/tb/exprpc2.h5"
h5ls(expr.path)
samp <- h5ls(expr.path, recursive = F)$name # 184
# batch = sub('.*_', '', samp)
# length(unique(batch)) # 38 too much to consider
genename <- h5read(expr.path, paste0(samp[1], "/feature"))
genename <- as.vector(genename)

pseudotime = readRDS("/hpc/group/jilab/hz/tensordata/tb/ptpc2.rds")
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

datlist <- vector(mode = 'list', length = length(samp))
names(datlist) <- samp
for (i in samp) {
  expr.sample <- h5read(expr.path, paste0(i, "/expr"))
  rownames(expr.sample) <- genename
  cell.sample <- h5read(expr.path, paste0(i, "/barcode"))
  pseudotime.sample <- pseudotime[cell.sample]
  expr.sample.agg <- agg_median(expr = expr.sample, pseudotime = pseudotime.sample)
  colnames(expr.sample.agg) <- paste(i, colnames(expr.sample.agg), sep = ':')
  datlist[[i]] <- expr.sample.agg
  print(i)
  gc()
}

expr.agg <- do.call(cbind, datlist)
str(expr.agg)
pseudotime.agg <- as.numeric(sub('.*:', '', colnames(expr.agg)))
names(pseudotime.agg) <- colnames(expr.agg)
str(pseudotime.agg)
length(unique(sub(':.*', '', names(pseudotime.agg))))

expr <- expr.agg
pseudotime <- pseudotime.agg
save(expr, pseudotime, file = "tb/tb_median.Rdata")
