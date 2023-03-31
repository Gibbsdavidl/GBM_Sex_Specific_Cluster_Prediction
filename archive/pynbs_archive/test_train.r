library(tidyverse)
library(robencla)

arr_f <- read_csv('data/Females_Array_Data_CW.csv')
rna_f <- read_csv('data/Females_RNASeq_Data_CW.csv')
load("data/Females_feature_sel_genelist.rda")

idx <- sample(1:nrow(rna_f),size=5,replace = F)
sampleidx <- rna_f$sample[idx]

dfa3 <- arr_f[!(arr_f$sample %in% sampleidx), ]
dfb3 <- rna_f[c(idx), ]

dim(dfa3)
dim(dfb3)

# our classifier object named Anne.
anne <- Robencla$new("Anne")

# xgboost parameters to pass to each sub-classifier in the ensembles
params <- list(max_depth=12,    # "height" of the tree, 6 is actually default. I think about 12 seems better.  (xgboost parameter)
               eta=0.2,        # this is the learning rate. smaller values slow it down, more conservative   (xgboost parameter)
               nrounds=24,     # number of rounds of training, lower numbers less overfitting (potentially)  (xgboost parameter)
               nthreads=4,     # parallel threads
               gamma=0.2,        # Minimum loss reduction required to again partition a leaf node. higher number ~ more conservative (xgboost parameter)
               lambda=1.2,     # L2 regularization term on weights, higher number ~ more conservative (xgboost parameter)
               alpha=0.2,     # L1 regularization term on weights. higher number ~ more conservative (xgboost parameter)
               verbose=0,
               train_perc=0.8,
               combine_function='median',
               size=11)

print('training')

# First we use the training data
anne$autotrain(data_frame = dfa3,
             label_name='cluster_group',
             sample_id = 'sample',
             data_mode=c('pairs'),  # pairs,sigpairs,quartiles,tertiles,binary,ranks,original #
             pair_list=genelist,  # subset to these genes.
             params=params)

print('testing')

anne$autotest(data_frame = dfb3, #[,colnames(dfa2)], #'data/Males_RNASeq_Data_v3.csv',
            label_name='cluster_group',
            sample_id = 'sample')


