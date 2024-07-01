

#!/usr/bin/env Rscript
args = c('1', '0.005') #commandArgs(trailingOnly=TRUE)

## MALES OR FEMALES
dat <- read.csv('../data/jive_training_array_data_F.csv.gz')

# drop the rownames
dat <- dat[,-1]

# We will work on one cluster at a time.
labels <- c('cluster1','cluster2','cluster3','cluster4','cluster5')

# number of genes in the table
ngenes <- ncol(dat)

# the pair count/index
k <- 1

# list of results
resl <- list()

#for (this_label in labels) {
this_label <- labels[as.numeric(args[1])]

# difference in props
this_thresh <- as.numeric(args[2])

#  index into this cluster
idx <- dat$ClusterLabel == this_label
jdx <- dat$ClusterLabel != this_label

# data labeled as THIS cluster
labdat <- dat[idx,]
# data labels as NOT THIS cluster
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
    if (abs(prop1-prop2) > this_thresh) {

          # next line: distnance between the two genes, within sample grouping.
          dist1 <- (labdat[,i] - labdat[,j])
          dist2 <- (notdat[,i] - notdat[,j])

          t_dist1   <- t.test(dist1, na.rm = T)$statistic
          med_dist1 <- median(dist1, na.rm = T)
          t_dist2   <- t.test(dist2, na.rm = T)$statistic
          med_dist2 <- median(dist2, na.rm=T)
          n1        <- length(dist1)
          n2        <- length(dist2)

          #capture this result        
          resl[[k]] <- c(this_label,
                         k,i,j,
                         prop1,prop2,
                         (prop1-prop2),
                         t_dist1, t_dist2,
                         med_dist1, med_dist2,
                         n1, n2) # want diff to be 1 or -1
          
          
          k <- k+1
    }
  }
} # end gene pair loop


# transform the list into a data frame.
resdf <- as.data.frame(do.call('rbind', resl))

#print(dim(resdf))

# write out the results.
colnames(resdf) <- c('cluster','pair','gene_i','gene_j',
                     'prop_this_cluster','prop_not_cluster','prop_diff',
                     't_test_cluster','t_test_not_cluster',
                     'med_dist_cluster','med_dist_not_cluster', 'n1', 'n2')

save(resdf, file=paste0('/users/dgibbs/proj/gbm_sex/results/female_resdf_',this_label,'.rda'))

write.csv(resdf, file=paste0('/users/dgibbs/proj/gbm_sex/results/feature_pairs_females_',this_label,'.csv'), quote=F)

