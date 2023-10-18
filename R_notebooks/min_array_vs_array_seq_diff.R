
# what to do about the min array values that max seq?


library(tidyverse)
library(hexbin)

resdf_m <- read_csv('../../data/feature_pairs_males_v3.csv')
dfa_m <- read_csv('../../data/Males_CV_Data_v3.csv')
dfb_m <- read_csv('../../data/Males_RNASeq_Data_v3.csv')

dfa <- dfa_m[dfa_m$sample %in% dfb_m$sample, ]
dfb <- dfb_m
rownames(dfa) <- dfa$sample
rownames(dfb) <- dfb$sample
dfa <- dfa[dfb$sample,]


### PLOT ###
gene <-  'RNF13' #'SOX11' #  'SEMA5A'  # PTCD1  LPIN1
qplot(x=as.vector(dfa[,gene])[[1]], y=as.vector(dfb[,gene])[[1]], 
      color=dfa$cluster.group, ylab='RNA-seq', xlab='Array', main=gene)



mean_bin_diff <- function(dfa, dfb, gene, numbins) {

  #first get values
  x <- as.vector(dfa[,gene])[[1]] 
  y <- as.vector(dfb[,gene])[[1]]
  
  # then scale to 0-1
  xs <- x/max(x) # scale(x) # 
  ys <- y/max(y) # scale(y) # 
  
  # then bin the values and get the index to binlevel
  cx <- cut_number(xs,n=numbins)
  
  
  res0 <- sapply (1:numbins, function(bi){
    l1 <- levels(cx)[bi]
    mean( abs( xs[cx==l1] - ys[cx==l1] ) )
  })
  
  res1 <- sapply (1:numbins, function(bi){
    l1 <- levels(cx)[bi]
    IQR( x[cx==l1] ) 
  })
  
  # then get mean diff between the two.  
  list(MinVal=min(dfa[,gene]), BinDiff=res0, IQRs=res1)
}



bigres0 <- lapply (2:10343, function(i) {
  gene <- colnames(dfa)[i]
  mean_bin_diff(dfa,dfb,gene,5)
})

names(bigres0) <- colnames(dfa)[2:10343]


# plot ... min val vs bin1 diff
x <- unlist( lapply(bigres0, function(a) a[[1]]) )
y <- unlist( lapply(bigres0, function(a) a[[2]][1]))
z <- unlist( lapply(bigres0, function(a) a[[3]][1]))

idx <- which(z > 0.5)

qplot(x=x[idx],y=y[idx],color=z[idx],
      xlab='Minimum array expression level', ylab='scaled difference between rna-seq and array')



bigres0[ which(y < 0.2)[4] ] 


