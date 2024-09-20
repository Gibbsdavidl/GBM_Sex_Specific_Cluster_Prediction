devtools::install_github("gibbsdavidl/robencla", ref = 'v0.3.4')

library(robencla)
library(readr)

dat.f <- read_csv('jive_training_array_data_F_v2.csv')

load('females_genepairs.rda')

# removing signature genes that are not available in our data sets
#female.jive.gene.list <- readRDS("../data/Female_jive_cluster_genes.rds") %>% dplyr::filter(value %in% colnames(dat.f))
#female.jive.gene.list <-  female.jive.gene.list[-c(18,289),]  # "RP11-35N6_1"

metrics_list <- list()
i <- 1

for (md in c(3,6,8,10,12,16,18)) {
  for (ei in c(0.05, 0.1, 0.2, 0.3, 0.5) ) {
    for (ni in c(3, 6, 12, 18, 24, 32)) {
      for (li in c(1.0, 1.2, 1.5, 1.8, 2)) {
        for (ai in c(0.0, 0.2, 0.5, 0.8, 1.0)) {
          
          # our classifier object named Roberta
          buffy <- Robencla$new("buffy")
          
          # xgboost parameters to pass to each sub-classifier in the ensembles
          params <- list(max_depth=md,   # "height" of the tree, 6 is actually default. I think about 12 seems better.  (xgboost parameter)
                         eta=ei,        # this is the learning rate. smaller values slow it down, more  conservative   (xgboost parameter)
                         nrounds=ni,     # number of rounds of training, lower numbers less overfitting (potentially)  (xgboost parameter)
                         early_stopping_rounds=2,
                         nthreads=8,     # parallel threads
                         gamma=0.2,      # Minimum loss req'd to again partition a leaf node. higher number ~ more conservative (xgboost)
                         lambda=1.2,     # L2 regularization term on weights, higher number ~ more conservative (xgboost parameter)
                         alpha=0.2,      # L1 regularization term on weights. higher number ~ more conservative (xgboost parameter)
                         verbose=0,
                         train_perc=0.8,
                         combine_function='median',
                         size=11
          )
          
          
          # First we use the training data
          buffy$autocv(data_frame=dat.f,
                       label_name='ClusterLabel',
                       sample_id = 'sample',
                       drop_list = 'Sex',
                       data_mode=c('pairs'), # pairs, allpairs, sigpairs, quartiles, tertiles, binarize, ranks, original #
                       signatures=NULL,
                       pair_list=genepairs,  # subset to these genes.
                       params=params,
                       cv_rounds=10
          )
          
          metrics <- buffy$classification_metrics()
          
          metrics_list[[i]] <- list(params, metrics)
          i <- i+1
          print(i)
        }
        
      }
    }
  }
}
