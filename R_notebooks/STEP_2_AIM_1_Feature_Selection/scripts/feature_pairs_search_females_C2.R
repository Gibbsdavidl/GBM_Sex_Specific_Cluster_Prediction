
#library(tidyverse)

dat <- read.csv('version3/Females_CV_Data_v3.csv')

## We will work on one cluster at a time.
labels <- c('cluster1','cluster2','cluster3','cluster4','cluster5')

# number of genes in the table
ngenes <- ncol(dat)

# the pair count/index
k <- 1

# list of results
resl <- list()

#for (this_label in labels) {
this_label <- labels[2]

    #  index into this cluster
    idx <- dat$cluster.group == this_label
    jdx <- dat$cluster.group != this_label
    
    labdat <- dat[idx,]
    notdat <- dat[jdx,]

    # number of samples in each category
    m <- sum(idx)
    n <- sum(jdx)

    # for each pair of genes
    for (i in 2:(ngenes-1)){
      for (j in (i+1):ngenes) {

        # proportion of samples showing pattern gene1 > gene2    
        prop1 <- (sum(labdat[,i] > labdat[,j])/m) # proportion i>j in this cluster
        prop2 <- (sum(notdat[,i] > notdat[,j])/n) # proportion i>j in others

        # if the difference in proportions, between cluster and not-cluster
        # is greater than 75%, then record the gene pair.
        if (abs(prop1-prop2) > 0.66) {

              # next line: distnance between the two genes, within sample grouping.
              dist1 <- sum(labdat[,i] - labdat[,j])
              dist2 <- sum(notdat[,i] - notdat[,j])

              #capture this result        
              resl[[k]] <- c(this_label,k,i,j,
                             prop1,prop2,
                             dist1, dist2,
                             (prop1-prop2)) # want diff to be 1 or -1
              print(k)
              k <- k+1
        } 
      }
    } # end gene pair loop
      # going to new cluster

# transform the list into a data frame.
resdf <- as.data.frame(do.call('rbind', resl))

print(dim(resdf))

# write out the results.
colnames(resdf) <- c('cluster','pair','gene_i','gene_j','prop_this_cluster','prop_not_cluster',
                     'dist_ij_cluster','dist_ij_not_cluster','prop_diff')
save(resdf, file=paste0('female_resdf_',this_label,'.rda'))
write.csv(resdf, file=paste0('feature_pairs_female_',this_label,'.csv'))
