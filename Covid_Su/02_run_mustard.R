load('Covid_Su/covid_su_median_scaled.Rdata')
meta = readRDS('Covid_Su/data/meta_sample.rds')
source('function/scFTSVD.R')
str(expr_sc)
str(pseudotime)
cellanno = sub(':.*', '', names(pseudotime))
length(unique(cellanno))

samp = rownames(meta)[meta$Time == 'First']
# samp = samp[stringr::str_detect(cellanno, '[.]2', negate = T)]
sum(cellanno %in% samp) # 3,719

start_time <- Sys.time()
system.time( res_ftsvd <- ftsvd(expr = expr_sc, pseudotime = pseudotime[cellanno %in% samp],
                                r = 5, smooth = 1e-3) )[3]/60
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken)

saveRDS(res_ftsvd, file = "Covid_Su/covid_su_res.rds")
