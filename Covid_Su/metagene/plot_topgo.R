res_ftsvd <- readRDS("Covid_Su/covid_su_res.rds")
# revise the signs of sample loadings and gene loadings Comp 3
res_ftsvd$A.hat[, 3] <- res_ftsvd$A.hat[, 3] * (-1)
res_ftsvd$B.hat[, 3] <- res_ftsvd$B.hat[, 3] * (-1)

gene.embedding <- sweep(res_ftsvd$B.hat, 2, res_ftsvd$Lambda, "*")

hl <- hclust(as.dist(1 - cor(Matrix::t(gene.embedding))), method = "average")
gene.clu <- cutree(hl, k = 5)
table(gene.clu)

suppressPackageStartupMessages(library(topGO))

diffgeneList <- sapply(sort(unique(gene.clu)), function(i) { names(gene.clu)[gene.clu == i] })
expr = readRDS("Covid_Su/data/saver_log2norm_sub.rds")

sep = ':.*'
resList <- lapply(diffgeneList, function(diffgene) {
  allgene <- rownames(expr)
  gl <- sub(sep, '', diffgene)
  back <- sub(sep, '', allgene)
  geneList <- factor(as.integer(back %in% gl))
  names(geneList) <- back
  suppressMessages({
  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
                geneSel = function(a) {a}, annot = annFUN.org,
                mapping = "org.Hs.eg.db", ID = "Symbol")
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  sigres <- GenTable(GOdata, classicFisher = resultFisher,
                     topNodes = length(resultFisher@score),
                     orderBy = "classicFisher", numChar = 1000)
  })
  sigres$classicFisher[sigres$classicFisher == "< 1e-30"] <- 0
  sigres <- sigres[sigres$Annotated >= 10, ]
  sigres$FDR <- p.adjust(sigres$classicFisher, method = "fdr")
  
  ptcount <- 0
  fc <- ((sigres[, "Significant"] + ptcount) / (sum(GOdata@allScores[GOdata@feasible] == 1) + ptcount)) /
    ((sigres[, "Annotated"] + ptcount) / (sum(GOdata@feasible) + ptcount))
  sigres <- data.frame(sigres, FC = fc)
  sigres <- sigres[order(sigres$FDR, -sigres$FC), ]
})
names(resList) <- sort(unique(gene.clu))

suppressPackageStartupMessages(library(ggplot2))
plotGOEnrich <- function(goRes, n = 5, sortByFDR = TRUE, fdr.cutoff = 0.05, fc.cutoff = 2){
  d <- lapply(names(goRes), function(i) {
    cbind(Cluster = i, goRes[[i]])
  })
  d <- do.call(rbind, d)    
  if (sortByFDR){
    d <- d[order(d$FDR, -d$FC), ]
  } else {
    d <- d[order(-d$FC, d$FDR), ]
  }
  cd <- do.call(rbind, sapply(sort(unique(d$Cluster)), function(i) {
    tmp <- d[d$Cluster == i, ]
    tmp <- tmp[tmp$FDR < fdr.cutoff & tmp$FC > fc.cutoff, , drop = F]
    if (nrow(tmp) > 0) {
      tmp[seq_len(min(n, nrow(tmp))), ]
    } else {
      NULL
    }
  }, simplify = F))
  ut <- unique(cd$Term)
  
  d <- d[d$Term %in% ut, c("Cluster", "Term", "FDR", "FC")]
  d$Cluster <- as.numeric(as.character(d$Cluster))
  
  d <- d[d$FDR < fdr.cutoff & d$FC > fc.cutoff, , drop = F]
  
  d$enc <- T
  fulld <- expand.grid(1:length(goRes), unique(d$Term))
  fulld <- fulld[!paste0(fulld[, 1], fulld[, 2]) %in% paste0(d[, 1], d[, 2]), ]
  fulld <- data.frame(Cluster = fulld[, 1], Term = as.character(fulld[, 2]), FDR = 1, FC = 0, enc = F)
  pd <- rbind(d, fulld)
  
  v1 <- sapply(unique(pd$Term), function(i) mean(pd$Cluster[pd$Term == i & pd$enc]))
  v2 <- sapply(unique(pd$Term), function(i) min(pd$Cluster[pd$Term == i & pd$enc]))
  
  pd$Term <- factor(as.character(pd$Term), levels = names(v2)[rev(order(v2, v1))])
  pd$enc <- ifelse(pd$enc, "Significant", "Non-significant")
  pd$Cluster <- as.character(pd$Cluster)
  
  FDR.label <- c("< 0.001", "< 0.01", "< 0.05", "> 0.05")
  # FDR.clv <- c(colorRampPalette(c("red", "mistyrose"))(3)[1:3], "grey")
  FDR.clv <- c("red3", "red", "mistyrose", "grey")
  names(FDR.clv) <- FDR.label
  pd$FDR.range <- cut(pd$FDR, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), right = F, labels = FDR.label)
  
  p <- ggplot(pd, aes(x = Cluster, y = Term, color = FDR.range, size = FC)) +
    geom_point() + theme_classic() + labs(x = "Module", y = NULL, size = "Fold change", color = "FDR") +
    scale_color_manual(values = FDR.clv) +
    theme(legend.position = "right") +
    scale_y_discrete(position = "left", labels = function(x) stringr::str_wrap(x, width = 50))
  return(p)
}

plotGOEnrich(goRes = resList)
ggsave("Covid_Su/figure/metagene/topgo_module.pdf", width = 5, height = 4)
