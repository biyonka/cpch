require(MASS)
library(robustbase)
library(ggplot2)
library(GEOquery)

setwd("~/cpch/DMD") #change directory to where the .rda and .csv data files are located

## Much of the following code is adapted heavily from pre-processing code written by
#Jingshu Wang, which was used to analyze these datasets for the AdaFilter paper
#https://arxiv.org/pdf/1610.03330.pdf
load("gds214.rda")
load("gds3027.rda")
load("gds1956.rda")
load("gds563.rda")

#download gene annotations that tell you about gene function
gpl246 <- read.csv("GPL246_annot.csv", sep = "\t")
gpl8300 <- read.csv("GPL8300_annot.csv", sep = "\t")
gpl96 <- read.csv("GPL96_annot.csv", sep = "\t")
gpl246 <- gpl246[order(gpl246$ID), ]
gpl96 <- gpl96[order(gpl96$ID), ]


#design matrices
X.3027 = metadata.3027$sample.X
X.563 <- metadata.563$sample.X
X.214 <-  metadata.214$sample.X
X.1956 <- metadata.1956$sample.X


## get preprocessed data from GDS data matrices
getYData <- function(gds.name, gplmat) {
  gds <- getGEO(gds.name)
  eset <- GDS2eSet(gds, do.log2= T)
  Ymat <- exprs(eset) 
  Ymat <- Ymat[order(rownames(Ymat)), ]
  Ymat.ori <- Ymat
  
  ## missing value imputation
  if (sum(is.na(Ymat)) > 0) {
    keep.idx <- apply(Ymat, 1, function(v) sum(is.na(v))) < ncol(Ymat) * 0.5
    Ymat <- Ymat[keep.idx, ]
    gplmat <- gplmat[keep.idx, ]
    require(impute)
    Ymat <- impute.knn(Ymat)$data
  }
  
  require(preprocessCore)
  ## background correction and log-transformation
  temp <- log(rma.background.correct(exp(Ymat)))
  ## quantile normalization
  temp <- normalize.quantiles(temp)
  colnames(temp) <- colnames(Ymat)
  rownames(temp) <- rownames(Ymat)
  Ymat <- temp
  
  
  
  return(list(Ymat = Ymat, gplmat= gplmat, Ymat.raw = Ymat.ori))
  
}


getP <- function(X.data, Y.data, do.cate = T,
                 r = -1,
                 contrasts = NULL) {
  
  if (do.cate) {
    require(cate)
    r2 <- est.confounder.num(~group|. - group, X.data, Y.data, method = "ed")
    print(r2)
    if (r < 0)
      r <- r2
    print(r)
    result <- cate(~group|. - group, X.data, Y.data, r = r)
    p1 <- result$beta.p.value
  } else
    p1 <- NULL
  
  
  library(limma)
  model.x <- model.matrix(~., X.data)
  fit1 <- lmFit(Y.data, model.x)
  if (is.null(contrasts))
    fit2 <- contrasts.fit(fit1, c(0, 1))
  else
    fit2 <- contrasts.fit(fit1, contrasts)
  
  limma.result <- eBayes(fit2)
  t2 <- limma.result$t
  t2 <- (t2 - median(t2)) / mad(t2 - median(t2))
  p2 <- 2 * pnorm(-abs(t2))
  return(list(z.stat = t2, p.limma = p2))
}



## reshape the annotation matrix
gplReshape <- function(gplmat) {
  gplmat <- data.frame(idx = 1:nrow(gplmat), ID = gplmat$ID, 
                       Gene.symbol = gplmat$Gene.symbol)
  gplmat <- gplmat[gplmat$Gene.symbol != "", ]
  dual.idx <- grep("///", gplmat$Gene.symbol)
  gplmat1 <- gplmat[-dual.idx, ]
  
  newdata <- gplmat[0, ]
  for (i in dual.idx) {
    symbol <- as.character(gplmat[i, 3])
    symbols <- strsplit(symbol, "///")[[1]]
    n <- length(symbols)
    temp <- gplmat[i, ]
    temp1 <- temp
    for (k in 1:(n-1))
      temp <- rbind(temp, temp1)
    temp$Gene.symbol <- symbols
    newdata <- rbind(newdata, temp)
  }
  
  gplmat1 <- rbind(newdata, gplmat1)
  gplmat1 <- gplmat1[order(gplmat1$Gene.symbol), ]
  return(gplmat1)
  
}


CompPVec_new <- function(X.data, Y.data, gplmat,
                         contrasts = NULL, r= -1,
                         use.limma = F) {
  if (use.limma)
    do.cate <- F
  else
    do.cate <- T
  
  result <- getP(X.data, Y.data, contrasts = contrasts, r = r,
                 do.cate = do.cate)
  if (use.limma)
    pvalue <- as.vector(result$p.limma)
  else 
    pvalue <- as.vector(result$p.cate)
  #  print(head(pvalue))
  pmat <- data.frame(ID1 = rownames(result$p.limma), z.stat = result$z.stat, 
                     pvalue = pvalue)
  gplmat['ID_check'] = gplmat$ID
  temp <- merge(pmat, gplmat, by.x = 'ID1', by.y = 'ID', all.y = TRUE)
  
  print(temp)
  pvec <-temp[, c('pvalue', 'z.stat', 'Gene.symbol')]
  pvec = pvec[order(pvec$Gene.symbol),]
  symbols <- pvec$Gene.symbol
  uni.symbols <- unique(symbols)
  p.group <- c()
  z.group <- c()
  index <- 1
  temp.group <- c()
  #for each unique gene, create a temp group that keeps track of the p-values that are repeated
  for (i in uni.symbols) {
    #  print(i)
    while (symbols[index] == i) {
      temp.group <- c(temp.group, pvec[index, 2]) #1
      index <- index + 1
      if (index > nrow(pvec))
        break
    }
    #average corresponding z-scores and rescale appropriately
    if (sum(is.na(temp.group)) == length(temp.group)){
      p.group = c(p.group, NA)
      z.group = c(z.group, NA)
    }
    else {
      k = length(temp.group[!is.na(temp.group)])
      z.agg = (k^0.5) * mean(temp.group, na.rm = TRUE)
      p <- 2 * pnorm(-abs(z.agg)) 
      p.group <- c(p.group, p)
      z.group = c(z.group, z.agg)
    }  
    temp.group <- c()
  }
  
  pvec <- data.frame(Gene.symbol = uni.symbols,
                     pvalue = p.group, zstat = z.group)
  
  
}
### compute p-values and combine by gene symbol
gpl246 <- gplReshape(gpl246)
gpl8300 <- gplReshape(gpl8300)
gpl96 <- gplReshape(gpl96)


#############
# w Avging #
############

r <- -1
use.limma <- T


p1 <- CompPVec_new(X.214, Y.214, gpl246, r= r, use.limma = use.limma)
p2 <- CompPVec_new(X.563, Y.563, gpl8300, r= r, use.limma = use.limma)
p3 <- CompPVec_new(X.1956, Y.1956, gpl96, r=r, use.limma = use.limma)
p4 <- CompPVec_new(X.3027, Y.3027, gpl96, contrasts = c(0, 1, 0, 0), r=r,
                   use.limma = use.limma)
p34 <- merge(p3, p4, by = "Gene.symbol")
p.list <- list(p1, p2, p3, p4)

pmat_test <- merge(merge(p1, p2, by = "Gene.symbol", all = T), 
                   p34, by = "Gene.symbol", all = T)

#############
# Run Tests #
#############
r = 4
q = 0.05

library(adaFilter)
num_shared_test = pmat_test[complete.cases(pmat_test),]
rownames(num_shared_test) = num_shared_test$Gene.symbol
num_shared_test = num_shared_test[, c(2, 4, 6, 8)]


result <-adaFilter(num_shared_test, r, alpha = q)
#View(result)
discoveries = result[result$decision == 1,]
print(nrow(discoveries))

#run classical methods
result_f = ClassicalMTP(num_shared_test, r = r, alpha = 0.05, method = "Fisher")
sum(result_f$decision == 1)
result_s = ClassicalMTP(num_shared_test, r = r, alpha = 0.05, method = "Simes")
sum(result_s$decision == 1)
result_b = ClassicalMTP(num_shared_test, r = r, alpha = 0.05, method = "Bonferroni")
sum(result_b$decision == 1)

#look at p-value histogram
hist(result_s$PC.pvalues)
hist(result_f$PC.pvalues)
hist(result_b$PC.pvalues)

colnames(pmat_test) = c('Gene.symbol', 'pvalue.1', 'ts.1', 'pvalue.2', 'ts.2', 'pvalue.3', 'ts.3','pvalue.4', 'ts.4' )
pmat_to_save = pmat_test[complete.cases(pmat_test),]
pmat_to_save['ada_decision_r_4'] = result$decision
hist(pmat_to_save$ts.4)

#change directory to where you want the csv to be saved
write.csv(pmat_to_save, '~/cpch/DMD/dmd.csv', row.names = FALSE) 


