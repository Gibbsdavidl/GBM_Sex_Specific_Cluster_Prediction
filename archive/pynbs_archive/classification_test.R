
#devtools::install_github('gibbsdavidl/robencla')

library(tidyverse)
library(robencla)

###################################################################################

### PULLING OUT SOME GENES
#setwd("E:/Work/GBM_clusters/feature_selection")
dfa <- read.csv('data/Males_CV_Data_v3.csv')
dfb <- read.csv('data/Males_RNASeq_Data_v3.csv')
resdf <- read.csv('data/feature_pairs_males_v3.csv')

#qplot(x=prop_diff, y=prop_this_cluster,  data=resdf, colour=cluster, alpha=0.5) #  %>% dplyr::filter(cluster=='cluster4')

#CLUSTER1
resdf %>% dplyr::filter(cluster=='cluster1') %>% arrange( desc(abs(dist_ij_cluster))) %>% head(n=20)
resdf %>% dplyr::filter(cluster=='cluster1') %>% arrange( desc(abs(prop_diff))) %>% head(n=20)
##CLUSTER2
resdf %>% dplyr::filter(cluster=='cluster2') %>% arrange( desc(abs(dist_ij_cluster))) %>% head()
resdf %>% dplyr::filter(cluster=='cluster2') %>% arrange( desc(abs(prop_diff))) %>% head()
#CLUSTER3
resdf %>% dplyr::filter(cluster=='cluster3') %>% arrange( desc(abs(dist_ij_cluster))) %>% head(n=10)
resdf %>% dplyr::filter(cluster=='cluster3') %>% arrange( desc(abs(prop_diff))) %>% head(n=10)
#CLUSTER4
resdf %>% dplyr::filter(cluster=='cluster4') %>% arrange( desc(abs(dist_ij_cluster))) %>% head(n=15)
resdf %>% dplyr::filter(cluster=='cluster4') %>% arrange( desc(abs(prop_diff))) %>% head(n=15)
#CLUSTER5
resdf %>% dplyr::filter(cluster=='cluster5') %>% arrange( desc(abs(dist_ij_cluster))) %>% head(n=20)
resdf %>% dplyr::filter(cluster=='cluster5') %>% arrange( desc(abs(prop_diff))) %>% head(n=20)

cluster1 <- c(3575,9646,2437,3575,5683,9725,7907,9145,9717,9749)
cluster2 <- c(6684,9747,6403,8092,2825,6717,2825,6717,6403,8092)
cluster3 <- c(2825,6532,556,4376,6532,10253,556,1013,4787,5661)
cluster4 <- c(760,4638,1760,3489,3489,6061,8710,10057,128,232)
cluster5 <- c(424,664,664,3169,4351,7775,5568,10102,3538,5568,4351,9320,9536,10089)


genelist <- c()

genelist <-  (c(
  colnames(dfa)[cluster1],
  colnames(dfa)[cluster2],
  colnames(dfa)[cluster3],
  colnames(dfa)[cluster4],
  colnames(dfa)[cluster5]
))

#########################################################################################


##allgenes <- unique(c(genelist, unlist(sigs)))

#dfa2 <- dfa[,c('sample','cluster.group',allgenes)]  # (dfa$sample %in% dfb$sample)
#dfb2 <- dfb[,c('sample','cluster.group',allgenes)]

# our classifier object named Anne.
anne <- Robencla$new("Anne")

# xgboost parameters to pass to each sub-classifier in the ensembles
params <- list(max_depth=12,    # "height" of the tree, 6 is actually default. I think about 12 seems better.  (xgboost parameter)
               eta=0.2,        # this is the learning rate. smaller values slow it down, more conservative   (xgboost parameter)
               nrounds=24,     # number of rounds of training, lower numbers less overfitting (potentially)  (xgboost parameter)
               nthreads=5,     # parallel threads
               gamma=0.2,        # Minimum loss reduction required to again partition a leaf node. higher number ~ more conservative (xgboost parameter)
               lambda=1.2,     # L2 regularization term on weights, higher number ~ more conservative (xgboost parameter)
               alpha=0.2,     # L1 regularization term on weights. higher number ~ more conservative (xgboost parameter)
               verbose=0,
               train_perc=0.8,
               combine_function='median',
               size=11
)

# First we use the training data
anne$autotrain(data_frame = dfa,
               label_name='cluster_group',
               sample_id = 'sample',
               data_mode=c('pairs'),  # pairs,allpairs, sigpairs,quartiles,tertiles,binary,ranks,original #
               signatures=NULL,
               pair_list=genelist,  # subset to these genes.
               params=params)

# now we apply the classifier to a test set.
anne$autotest(data_frame = dfb,#[,colnames(dfa2)], #'data/Males_RNASeq_Data_v3.csv',
              label_name='cluster_group',
              sample_id = 'sample')

table(Pred=anne$results()$BestCall, True=anne$test_label)

anne$classification_metrics(use_cv_results=FALSE)

# with small set 1
# True
# Pred       cluster1 cluster2 cluster3 cluster4 cluster5
# cluster1       13        1        0        2        1
# cluster2        1       12        0        1        2
# cluster3        1        3       16        2        3
# cluster4        1        0        0        2        0
# cluster5        3        0        1        0        7


resl <- list()

for (i in 1:20) {

  print(i)
  idx <- sample(1:72,size=5,replace = F)
  sampleidx <- dfb$sample[idx]

  dfa3 <- dfa[!(dfa$sample %in% sampleidx), ]
  dfb3 <- dfb[c(idx), ]

  # our classifier object named Anne.
  anne <- Robencla$new("Anne")

  # xgboost parameters to pass to each sub-classifier in the ensembles
  params <- list(max_depth=12,    # "height" of the tree, 6 is actually default. I think about 12 seems better.  (xgboost parameter)
                 eta=0.2,        # this is the learning rate. smaller values slow it down, more conservative   (xgboost parameter)
                 nrounds=24,     # number of rounds of training, lower numbers less overfitting (potentially)  (xgboost parameter)
                 nthreads=5,     # parallel threads
                 gamma=0.2,        # Minimum loss reduction required to again partition a leaf node. higher number ~ more conservative (xgboost parameter)
                 lambda=1.2,     # L2 regularization term on weights, higher number ~ more conservative (xgboost parameter)
                 alpha=0.2,     # L1 regularization term on weights. higher number ~ more conservative (xgboost parameter)
                 verbose=0,
                 train_perc=0.8,
                 combine_function='median',
                 size=11
  )

  # First we use the training data
  anne$autotrain(data_frame = dfa3,
                 label_name='cluster_group',
                 sample_id = 'sample',
                 data_mode=c('pairs'),  # pairs,sigpairs,quartiles,tertiles,binary,ranks,original #
                 signatures=NULL,
                 pair_list=genelist,  # subset to these genes.
                 params=params)

  # now we apply the classifier to a test set.
  anne$autotest(data_frame = dfb3,#[,colnames(dfa2)], #'data/Males_RNASeq_Data_v3.csv',
                label_name='cluster_group',
                sample_id = 'sample')

  #table(Pred=anne$results()$BestCall, True=anne$test_label)
  resl[[i]] <- data.frame(Pred=anne$results()$BestCall, True=anne$test_label)
}

resdf <- do.call("rbind", resl)
table(resdf$Pred, resdf$True)


# pairs only
#          cluster1 cluster2 cluster3 cluster4 cluster5
# cluster1      149       21        0       18       13
# cluster2       82      208      106       76       39
# cluster3       12        1      123        0       12
# cluster5       19        0       12        0      109
#
#
# tertile only
#          cluster1 cluster2 cluster3 cluster4 cluster5
# cluster1      133       39       61       22       35
# cluster2       86      198       31       66       34
# cluster3        0        0      138        0       19
# cluster4        1        0        0        0        0
# cluster5       19        3       18        0       97

