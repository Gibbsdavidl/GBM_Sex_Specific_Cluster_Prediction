library(tidyverse)
library(robencla)
library(survival)
library(survminer)



# Take the data set, and filter out pair genes that are not present.

clean_genepairs_list <- function(genepairs, dat_x) {
  genepairs_clean <- list()
  i <- 1
  for (pi in names(genepairs)) {
    genes <- genepairs[[pi]]
    missing_genes <- unique( genes [! genes %in% colnames(dat_x)] )
    print(paste0(pi, '  ', missing_genes))
    newpairlist <- c()
    for (j in seq.int(from=1,to=length(genes), by=2)) {
      if ( (genes[j] %in% missing_genes) | (genes[j+1] %in% missing_genes) ) {
        # then don't add them!
        print(paste0("removing from pair list:  ", genes[j], ' ', genes[j+1]))
      } else {
        newpairlist <- c(newpairlist, genes[j], genes[j+1])
      }
    }
    genepairs_clean[[pi]] <- newpairlist
  }
  return(genepairs_clean)
}

# pairlist is a list of vectors (paired terms)
# n is the number of pairs for each entry desired

shorter_plist <- function(pairlist, n) {
  newlist <- list()
  for (ni in names(pairlist)) {
    newlist[[ni]] <- pairlist[[ni]][1:(2*n)]
  }
  return(newlist)  
}


# read in the data
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2)
#arr_f <- read_csv('../data/jive_training_array_data_F.csv.gz')
#cpt_f <- read_csv('../data/cptac3_gbm_rnaseq_table_unstranded_Female.csv.gz')
#load('../data/tcga_array_data_F.rda')
#arr_f <- dat_f

#load('../results/results_no_array_min/females_genepairs.rda')
#gset1 <- genepairs

#load('../results/females_genepairs_val.rda')
#gset2 <- genepairs

#load('../results/females_genepairs_rna.rda')
#gset3 <- genepairs

# make sure the gene pairs are all present
gpairs <- shorter_plist(gset3, 4)  # gset2 is NOT predictive!!

# check they're present in the data
lapply(gpairs, function(x) x[! (x %in% colnames(arr_f ))] )
lapply(gpairs, function(x) x[! (x %in% colnames(cpt_f))] ) # missing "FABP5"  "CD24"   "LGALS3"

# clean the pairs
genepairs_cl <- clean_genepairs_list(gpairs, cpt_f) # only one with missing genes

gs1 <- c(unlist(unique(genepairs_cl)), "Barcode", "ClusterLabel")
gs2 <- c(unlist(unique(genepairs_cl)), "CaseID", "Survival", "Censored")


dat_f <- arr_f[,gs1]
val_f <- cpt_f[,gs2]


# "fix" the overall survival in days, paitents still alive are listed as "--"
val_f_surv <- as.numeric(val_f$Survival)
val_f <- val_f[val_f_surv > 10,]
val_f_surv <- as.numeric(val_f$Survival)
val_f$Survival <- val_f_surv

# what if ...
dat_f$ClusterLabel2 <- ifelse(dat_f$ClusterLabel == 'cluster3', yes = 'C3', no='notC3')
c3_pairs <- list()
c3_pairs[['C3']] <- genepairs_cl$cluster3
c3_pairs[['notC3']] <- c(genepairs_cl$cluster1, genepairs_cl$cluster2, genepairs_cl$cluster4, genepairs_cl$cluster5)


# xgboost parameters to pass to each sub-classifier in the ensembles
# xgboost parameters to pass to each sub-classifier in the ensembles
params1 <- list(max_depth=6,   # "height" of the tree, 6 is actually default. I think about 12 seems better.  (xgboost parameter)
                eta=0.3,        # this is the learning rate. smaller values slow it down, more  conservative   (xgboost parameter)
                nrounds=24,     # number of rounds of training, lower numbers less overfitting (potentially)  (xgboost parameter)
                early_stopping_rounds=2,
                nthreads=8,     # parallel threads
                gamma=1.2,      # Minimum loss req'd to again partition a leaf node. higher number ~ more conservative (xgboost)
                lambda=2.5,     # L2 regularization term on weights, higher number ~ more conservative (xgboost parameter)
                alpha=1.2,      # L1 regularization term on weights. higher number ~ more conservative (xgboost parameter)
                verbose=0,
                train_perc=0.8,
                combine_function='median',
                size=11
)



mod_f <- robencla::Robencla$new('f_model')

# mod_f$train(data_frame=dat_f,
#             label_name='ClusterLabel', #'ClusterLabel2',
#             sample_id = 'Barcode',
#             drop_list = c(),
#             data_mode=c('pairs'),    # pairs, allpairs, sigpairs, quartiles, tertiles, binarize, ranks, original #
#             signatures=NULL,         #
#             pair_list=genepairs_cl,  # c3_pairs   # subset to these genes.
#             params=params1
# )


mod_f$train(data_frame=dat_f,
            label_name='ClusterLabel',
            sample_id = 'Barcode',
            drop_list = c('ClusterLabel2'),
            data_mode=c('pairs'),    # pairs, allpairs, sigpairs, quartiles, tertiles, binarize, ranks, original #
            signatures=NULL,         #
            pair_list=genepairs_cl,  #    # subset to these genes.
            params=params1
)


# run the model
mod_f$predict(data_frame = val_f, 
              sample_id = 'CaseID',
              drop_list = c('Survival', 'Censored'))



# check the predictions
table(mod_f$results()$BestCall)

head(mod_f$results())

# reorder the tables
rownames(val_f) <- val_f$CaseID
pheno_reorg <- val_f[mod_f$results()$SampleIDs, c("CaseID","Survival","Censored")]

# now bind in the results to the phenotype data
resdf <- cbind(pheno_reorg, mod_f$results())
#resdf <- na.omit(resdf)
resdf$Survival <- as.numeric(resdf$Survival)

#resdf$C3 <- unlist(resdf$C3)
#resdf$notC3 <- unlist(resdf$notC3)
resdf$cluster1 <- unlist(resdf$cluster1)
resdf$cluster2 <- unlist(resdf$cluster2)
resdf$cluster3 <- unlist(resdf$cluster3)
resdf$cluster4 <- unlist(resdf$cluster4)
resdf$cluster5 <- unlist(resdf$cluster5)

resdf$CensoredCode <- ifelse(resdf$Censored == 'Alive', yes=0, no=1)

modfit <- survfit(Surv(Survival, CensoredCode)~BestCalls,data=resdf)
ggsurvplot(modfit, pval = T)

surp <- ggsurvplot(modfit, pval = T, xlim=c(0,1700) )
surp
ggsave("../figures/survplot_cptac3_4pairs_refined_features_F.pdf", surp$plot, height = 5, width = 10)

# scores for the samples called cluster3
#boxplot(as.numeric(resdf$cluster3)~as.factor(resdf$BestCalls))


#boxplot(as.numeric(resdf$C3)~as.factor(resdf$BestCalls))

