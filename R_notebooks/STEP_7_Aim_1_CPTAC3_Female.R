
#devtools::install_github("gibbsdavidl/robencla", ref = 'v0.3.4')

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
shorter_and_clean <- function(pairlist, n) {
  newlist <- list()
  for (ni in names(pairlist)) {
    newlist[[ni]] <- str_replace_all(pairlist[[ni]][1:(2*n)], '\\.', '_')
  }
  return(newlist)  
}


load('data/F_processed_data.rda')
load('results/females_genepairs.rda')

gpairs <- shorter_and_clean(genepairs, 3)  
# clean the pairs
genepairs_cl <- clean_genepairs_list(gpairs, dat_f)  # only one with missing genes
genepairs_cl <- clean_genepairs_list(gpairs, cpt_f)  # only one with missing genes


# xgboost parameters to pass to each sub-classifier in the ensembles
# xgboost parameters to pass to each sub-classifier in the ensembles
params1 <- list(max_depth=12,   # "height" of the tree, 6 is actually default. I think about 12 seems better.  (xgboost parameter)
                eta=0.45,        # this is the learning rate. smaller values slow it down, more  conservative   (xgboost parameter)
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


mod_f$train(data_frame=dat_f,
            label_name='ClusterLabel',
            sample_id = 'Barcode',
            drop_list = c(), #c('ClusterLabel2'),
            data_mode=c('pairs'),    # pairs, allpairs, sigpairs, quartiles, tertiles, binarize, ranks, original #
            signatures=NULL,         #
            pair_list=genepairs_cl,  #    # subset to these genes.
            params=params1
)


# run the model
mod_f$predict(data_frame = cpt_f, 
              sample_id = 'Barcode')

#"sample", "diag__days_to_last_known_disease_status", "demo__vital_status", "demo__gender"

# check the predictions
table(mod_f$results()$BestCall)
#cluster1 cluster2 cluster3 cluster5 
#14       15        2       13 

# reorder the tables
cpt_f2 <- as.data.frame(cpt_f)
rownames(cpt_f2) <- cpt_f2$Barcode
pheno_reorg <- cpt_f2[mod_f$results()$SampleIDs, c("Barcode", "diag__days_to_last_known_disease_status", "demo__vital_status", "demo__gender")]

# now bind in the results to the phenotype data
resdf <- cbind(pheno_reorg, mod_f$results())
#resdf <- na.omit(resdf)
resdf$Survival <- as.numeric(resdf$diag__days_to_last_known_disease_status)

#resdf$C3 <- unlist(resdf$C3)
#resdf$notC3 <- unlist(resdf$notC3)
resdf$cluster1 <- unlist(resdf$cluster1)
resdf$cluster2 <- unlist(resdf$cluster2)
resdf$cluster3 <- unlist(resdf$cluster3)
resdf$cluster4 <- unlist(resdf$cluster4)
resdf$cluster5 <- unlist(resdf$cluster5)

resdf$CensoredCode <- ifelse(resdf$demo__vital_status == 'Alive', yes=1, no=2)

modfit <- survfit(Surv(Survival, CensoredCode)~BestCalls,data=resdf)
ggsurvplot(modfit, pval = T)

surp <- ggsurvplot(modfit, pval = T, xlim=c(0,1700) )
surp

ggsave("figures/survplot_cptac3_F.pdf", surp$plot, height = 8, width = 8)

mod_f$importance()

# $cluster3
# # A tibble: 3 × 4
# Feature       MedGain MedCover MedFreq
# <chr>           <dbl>    <dbl>   <dbl>
# 1 RBP1_X_BMP2     0.671    0.551   0.444
# 2 C1QL1_X_FABP5   0.190    0.271   0.333
# 3 MT1M_X_FERMT1   0.126    0.187   0.25 

