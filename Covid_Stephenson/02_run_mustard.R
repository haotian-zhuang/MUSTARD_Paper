load("Covid_Stephenson/covid_stephenson_median_scaled.Rdata")
cellanno = sub(':.*', '', names(pseudotime))
length(unique(cellanno))
source("/hpc/group/jilab/hz/scFTSVD/function/scFTSVD.R")
start_time <- Sys.time()
system.time( res_ftsvd <- ftsvd(expr = expr_sc, pseudotime = pseudotime, cellanno = NULL,
                                r = 10, smooth = 1e-3) )[3]/60
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken)

saveRDS(res_ftsvd, file = "Covid_Stephenson/covid_stephenson_res.rds")
