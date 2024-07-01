

# try random forest

library(tidyverse)
library(randomForest)

load("E:/Work/GBM_clusters/GBM_JIVE_ANALYSIS/data/F_processed_data.rda")

# removing signature genes that are not available in our data sets
male.jive.gene.list <- readRDS("data/Males_jive_cluster_genes.rds") %>% dplyr::filter(value %in% colnames(dat_f))

genes <- male.jive.gene.list$value

# in case we want to test it #
idx <- sample(1:139, size=125, replace = F)
jdx <- setdiff(1:139, idx)
train <- dat_f[idx,genes]
test <- dat_f[jdx,genes]
y_train <- dat_f$ClusterLabel[idx]
y_test <- dat_f$ClusterLabel[jdx]

randomForest(x=train, y=as.factor(y_train),
             xtest=test, ytest=as.factor(y_test))

# Call:
#   randomForest(x = train, y = as.factor(y_train), xtest = test,      ytest = as.factor(y_test)) 
# Type of random forest: classification
# Number of trees: 500
# No. of variables tried at each split: 17
# 
# OOB estimate of  error rate: 16%
# Confusion matrix:
#   cluster1 cluster2 cluster3 cluster4 cluster5 class.error
# cluster1       34        0        0        0        3  0.08108108
# cluster2        0       17        0        1        5  0.26086957
# cluster3        1        0       10        0        2  0.23076923
# cluster4        2        2        0       13        1  0.27777778
# cluster5        0        3        0        0       31  0.08823529
# Test set error rate: 0%
# Confusion matrix:
#   cluster1 cluster2 cluster3 cluster4 cluster5 class.error
# cluster1        4        0        0        0        0           0
# cluster2        0        1        0        0        0           0
# cluster3        0        0        1        0        0           0
# cluster4        0        0        0        2        0           0
# cluster5        0        0        0        0        6           0


