
# what to do about the min array values that max seq?


library(tidyverse)

# read in array and rna-seq data
resdf_m <- read_csv('../results/feature_selection_table_male.csv')
dfa_m <- read_csv('../data/Males_Array_Data.csv.gz')
dfb_m <- read_csv('../data/Males_RNASeq_Data.csv.gz')
dfa_m$sample <- str_sub(dfa_m$sample, start=1, end=12)

# get them in order
dfa <- dfa_m[dfa_m$sample %in% dfb_m$sample, ]
dfb <- dfb_m
rownames(dfa) <- dfa$sample
rownames(dfb) <- dfb$sample
dfa <- dfa[dfb$sample,]


### PLOT ###
gene <-  'SOX11' #'SOX11' #  'SEMA5A'  # PTCD1  LPIN1
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
x <- unlist( lapply(bigres0, function(a) a[[1]]) )  # min value
y <- unlist( lapply(bigres0, function(a) a[[2]][1]))
z <- unlist( lapply(bigres0, function(a) a[[3]][1]))

idx <- which(z > 0.1)

qplot(x=x[idx],y=y[idx],
      xlab='Minimum array expression level', ylab='scaled difference between rna-seq and array') +
  geom_vline(xintercept = 4.6)


bigres0[ which(y < 0.2)[4] ] 


selected_genes <- names(x)[x > 4.6]


getgenecols <- function(gene_i, gene_j, dfa_m, x)  {
  gene_names <- colnames(dfa_m)[c(gene_i, gene_j)]
  return(gene_names)
}


makenewcols <- function(gene_i, gene_j, dfa_m, x)  {
  gene_names <- colnames(dfa_m)[c(gene_i, gene_j)]
  minvals <- x[gene_names]
  return(minvals)
}


minvallist <- lapply(1:nrow(resdf_m), function(i) {
  makenewcols(unlist(resdf_m[i, 'gene_i']), unlist(resdf_m[i, 'gene_j']), dfa_m, x)
})


genelist <- lapply(1:nrow(resdf_m), function(i) {
  getgenecols(unlist(resdf_m[i, 'gene_i']), unlist(resdf_m[i, 'gene_j']), dfa_m, x)
})


minvalmat <- do.call('rbind', minvallist)
genemat <- do.call('rbind', genelist)


newmat <- resdf_m
newmat['gene_i_symbol'] <- genemat[,1]
newmat['gene_j_symbol'] <- genemat[,2]
newmat['gene_i_min'] <- minvalmat[,1]
newmat['gene_j_min'] <- minvalmat[,2]




getpairs <- function(pairmat, cluster_c, thresh) {

    filtmat <- newmat[ (newmat$gene_i_min > thresh) & (newmat$gene_j_min > thresh) & (newmat$cluster == cluster_c), ]
    
    idx <- order(filtmat$dist_ij_cluster, decreasing = T)[1:25]
    jdx <- order(filtmat$dist_ij_cluster, decreasing = F)[1:25]

    upgenes <- c()
    for (i in idx) {
      upgenes <- c(upgenes, as.character(filtmat[i, 'gene_i_symbol']))
      upgenes <- c(upgenes, as.character(filtmat[i, 'gene_j_symbol']))
    }
    dngenes <- c()
    for (j in jdx) {
      dngenes <- c(dngenes, as.character(filtmat[j, 'gene_i_symbol']))
      dngenes <- c(dngenes, as.character(filtmat[j, 'gene_j_symbol']))
    }

    return( c(upgenes, dngenes) )    
}

genepairs <- list(
  'cluster1'=getpairs(newmat, 'cluster1', 4.6),
  'cluster2'=getpairs(newmat, 'cluster2', 4.6),
  'cluster3'=getpairs(newmat, 'cluster3', 4.6),
  'cluster4'=getpairs(newmat, 'cluster4', 4.6),
  'cluster5'=getpairs(newmat, 'cluster5', 4.6)
)


thresh <- 4.6
cluster_c <- 'cluster2'
filtmat <- newmat[ (newmat$gene_i_min > thresh) & (newmat$gene_j_min > thresh) & (newmat$cluster == cluster_c), ]









