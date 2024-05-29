# feature search

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


run_model <- function(train, test, genepairs_list, params, npairs) {
  
  gpairs <- shorter_plist(genepairs_list, npairs)  
  genepairs_cl <- clean_genepairs_list(gpairs, train)
  genepairs_cl <- clean_genepairs_list(genepairs_cl, test)  # only one with missing genes
  
  
  mod_f <- robencla::Robencla$new('f_model')
  
  
  mod_f$train(data_frame=train,
              label_name='ClusterLabel', #'ClusterLabel2',
              sample_id = 'Barcode',
              drop_list = c(), #c('Sex','ClusterLabel2'),
              data_mode=c('pairs'),    # pairs, allpairs, sigpairs, quartiles, tertiles, binarize, ranks, original #
              signatures=NULL,         # 
              pair_list=genepairs_cl,  # c3_pairs   # subset to these genes.
              params=params
  )
  
  # run the model
  mod_f$predict(data_frame = test, 
                sample_id = 'Barcode')
                #drop_list = c('Sex','Survival', 'Censored'))
  
  return(mod_f)
}



# name the data used
train <- dat_f
test <- val_f

# feature processing
load('../results/females_genepairs_rna.rda')
gset3 <- genepairs
gpairs <- shorter_plist(gset3, 12)  
lapply(gpairs, function(x) x[! (x %in% colnames(train ))] )
lapply(gpairs, function(x) x[! (x %in% colnames(test))] ) # missing "FABP5"  "CD24"   "LGALS3"
genepairs_cl <- clean_genepairs_list(gpairs, test)  # only one with missing genes

# run the model
mod_f <- run_model(train,test,genepairs_cl)

# check the preds
table(mod_f$results()$BestCall)

# build the predicted output table
resdf <- data.frame(Survival=mod_f$test_data$Survival, 
                    Censored=mod_f$test_data$Censored, 
                    Barcode=mod_f$test_data$Barcode,
                    Pred=mod_f$results()$BestCall)
head(resdf)

#resdf$CensoredCode <- ifelse(resdf$Censored == 'alive', yes=1, no=2)
resdf$CensoredCode <- resdf$Censored
modfit <- survfit(Surv(Survival, CensoredCode)~Pred,data=resdf)
ggsurvplot(modfit, pval = T, xlim=c(0,1700))

surp <- ggsurvplot(modfit, pval = T, xlim=c(0,1700) )
ggsave("../figures/survplot_validation_refined_features_F_8by8.pdf", surp$plot, height = 8, width = 8)


impdf <- mod_f$importance()

df1 <- cbind(data.frame(Data=c('TCGA Array'), Cluster=c('cluster1')), impdf$cluster1)
df2 <- cbind(data.frame(Data=c('TCGA Array'), Cluster=c('cluster2')), impdf$cluster2)
df3 <- cbind(data.frame(Data=c('TCGA Array'), Cluster=c('cluster3')), impdf$cluster3)
df4 <- cbind(data.frame(Data=c('TCGA Array'), Cluster=c('cluster4')), impdf$cluster4)
df5 <- cbind(data.frame(Data=c('TCGA Array'), Cluster=c('cluster5')), impdf$cluster5)

df <- rbind(df1,df2,df3,df4,df5)
write.csv(df, file = '../results/feature_importance/tcga_array.csv')



##########################
##########################

# train with Val_f

all(val_test_dat$Barcode == resdf$Barcode)  # TRUE
val_test_dat <- mod_f$test_data
val_test_dat$ClusterLabel <- resdf$Pred


# run the model
mod_f2 <- run_model(val_test_dat,train,genepairs_cl)



# build the predicted output table
resdf <- data.frame(Barcode=mod_f2$test_data$Barcode,
                    Pred=mod_f2$results()$BestCall,
                    True=mod_f2$test_data$ClusterLabel)
head(resdf)

table(resdf$True, resdf$Pred)

#         cluster1 cluster2 cluster3 cluster4 cluster5
#cluster1       34        3        1        2        1
#cluster2        1       20        0        1        2
#cluster3        0        0       13        0        1
#cluster4        2        0        0       18        0
#cluster5        1        1        0        1       37

impdf <- mod_f2$importance()

df1 <- cbind(data.frame(Data=c('Valid Array'), Cluster=c('cluster1')), impdf$cluster1)
df2 <- cbind(data.frame(Data=c('Valid Array'), Cluster=c('cluster2')), impdf$cluster2)
df3 <- cbind(data.frame(Data=c('Valid Array'), Cluster=c('cluster3')), impdf$cluster3)
df4 <- cbind(data.frame(Data=c('Valid Array'), Cluster=c('cluster4')), impdf$cluster4)
df5 <- cbind(data.frame(Data=c('Valid Array'), Cluster=c('cluster5')), impdf$cluster5)
df <- rbind(df1,df2,df3,df4,df5)

write.csv(df, file = '../results/feature_importance/validation_array.csv')


##########################
##########################


# train with tcga-rna-seq



# run the model
mod_f3 <- run_model(tcga_f,train,genepairs_cl)



# build the predicted output table
resdf <- data.frame(Barcode=mod_f3$test_data$Barcode,
                    Pred=mod_f3$results()$BestCall,
                    True=mod_f3$test_data$ClusterLabel)
head(resdf)

table(True=resdf$True, Pred=resdf$Pred)

#         Pred
#True     cluster1 cluster2 cluster3 cluster4 cluster5
#cluster1       37        2        0        2        0
#cluster2        6        9        0        1        8
#cluster3        1        0        8        0        5
#cluster4        0        0        0       20        0
#cluster5        3        3        0        3       31

impdf <- mod_f3$importance()
df1 <- cbind(data.frame(Data=c('TCGA RNAseq'), Cluster=c('cluster1')), impdf$cluster1)
df2 <- cbind(data.frame(Data=c('TCGA RNAseq'), Cluster=c('cluster2')), impdf$cluster2)
df3 <- cbind(data.frame(Data=c('TCGA RNAseq'), Cluster=c('cluster3')), impdf$cluster3)
df4 <- cbind(data.frame(Data=c('TCGA RNAseq'), Cluster=c('cluster4')), impdf$cluster4)
df5 <- cbind(data.frame(Data=c('TCGA RNAseq'), Cluster=c('cluster5')), impdf$cluster5)
df <- rbind(df1,df2,df3,df4,df5)

write.csv(df, file = '../results/feature_importance/tcga_rnaseq.csv')


##############################
##############################

# train with tempus



# name the data used
train <- dat_f
test <- tmp_f

# feature processing
load('../results/females_genepairs_rna.rda')
gset3 <- genepairs

# set model parameters
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


# run the model
mod_f <- run_model(train,test,gset3,params1,12)

# check the preds
table(mod_f$results()$BestCall)

#cluster1 cluster2 cluster3 cluster4 cluster5 
#38        3       10        9       26 

# build the predicted output table
resdf <- data.frame(Survival=mod_f$test_data$Survival, 
                    Censored=mod_f$test_data$Censored, 
                    Barcode=mod_f$test_data$Barcode,
                    Pred=mod_f$results()$BestCall)
head(resdf)

resdf$CensoredCode <- ifelse(resdf$Censored == 'alive', yes=1, no=2)
#resdf$CensoredCode <- resdf$Censored
modfit <- survfit(Surv(Survival, CensoredCode)~Pred,data=resdf)
ggsurvplot(modfit, pval = T, xlim=c(0,1700))


surp <- ggsurvplot(modfit, pval = T, xlim=c(0,1700) )
ggsave("../figures/survplot_validation_refined_features_F_8by8.pdf", surp$plot, height = 8, width = 8)


### now bind in the preds to the data

resdf$ClusterLabel <- resdf$Pred
tmp_f_pred <- dplyr::inner_join(tmp_f, resdf)
# Joining with `by = join_by(Barcode, Censored, Survival)`

# run the model
# name the data used
train <- tmp_f_pred
test <- dat_f
mod_f2 <- run_model(train,test,gset3,params1,12)


# build the predicted output table
resdf2 <- data.frame(Barcode=mod_f2$test_data$Barcode,
                    Pred=mod_f2$results()$BestCall,
                    True=mod_f2$test_data$ClusterLabel)
head(resdf2)

table(resdf2$True, resdf2$Pred)

#         cluster1 cluster3 cluster4 cluster5
#cluster1       33        1        3        4
#cluster2        8        0        0       16
#cluster3        0       13        0        1
#cluster4        4        0       11        5
#cluster5        0        0        1       39


impdf <- mod_f2$importance()

df1 <- cbind(data.frame(Data=c('Tempus'), Cluster=c('cluster1')), impdf$cluster1)
df2 <- cbind(data.frame(Data=c('Tempus'), Cluster=c('cluster2')), impdf$cluster2)
df3 <- cbind(data.frame(Data=c('Tempus'), Cluster=c('cluster3')), impdf$cluster3)
df4 <- cbind(data.frame(Data=c('Tempus'), Cluster=c('cluster4')), impdf$cluster4)
df5 <- cbind(data.frame(Data=c('Tempus'), Cluster=c('cluster5')), impdf$cluster5)

df <- rbind(df1,df2,df3,df4,df5)
write.csv(df, file = '../results/feature_importance/tempus.csv')

