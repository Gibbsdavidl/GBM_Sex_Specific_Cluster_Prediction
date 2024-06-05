

x <- matrix(data=rnorm(25*50), nrow = 50)


### gene names
cnames <- c()
for (i in 1:25) {
  g <- paste0('g',as.character(sample(1:1231241, 1)))
  cnames <- c(cnames,g)
}
x <- as.data.frame(x)
colnames(x) <- cnames

### sample names
rnames <- c()
for (i in 1:50) {
  s <- paste0('s',as.character(sample(1:1231241, 1)))
  rnames <- c(rnames,s)
}
x$sample <- rnames


### labels
l <- sample(1:2, size = 50, replace = T)
x$label <- l

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


mod1 <- robencla::Robencla$new('f_model')

mod1$autocv(data_frame = x,
            label_name='label', 
            sample_id = 'sample',
            data_mode=c('allpairs'),    # pairs, allpairs, sigpairs, quartiles, tertiles, binarize, ranks, original #
            cv_rounds = 5,
            pair_list=cnames,  # c3_pairs   # subset to these genes.
            params=params1)


mod1$cv_results
table(mod1$cv_results$BestCalls, mod1$cv_results$Label)

