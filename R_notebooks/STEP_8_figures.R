library(ggplot2)
library(pheatmap)
library(forcats)
library(dplyr)



load("E:/Work/GBM_clusters/GBM_JIVE_ANALYSIS/data/F_processed_data.rda")
load('results/genepair_features/females_genepairs_curated_across_datasets.rda')


load("E:/Work/GBM_clusters/GBM_JIVE_ANALYSIS/data/M_processed_data.rda")
load('data/m_genepairs.rda')

f_val_pred <- read_csv("results/F_validation_predictions_using_rna_pairs.csv", 
                       col_types = cols(...1 = col_skip()))


f_tmp_pred <- read_csv("results/F_Tempus_results.csv", 
                       col_types = cols(...1 = col_skip()))


##################################################

# figure 1
annot_f <- read.csv('results/annotdf_f.csv.gz')

df <- data.frame()
for (cl in unique(annot_f$cluster)) {
  df_new <- data.frame(Cluster=cl, Score=sort(annot_f$score[annot_f$cluster == cl ]))  # & abs(annot_f$prop_diff) > 0.5 same plot
  df_new$Feature_Pair <- 1:nrow(df_new)
  df <- rbind(df, df_new)
}

# Reorder following the value of another column:
df %>%
  ggplot( aes(x=Feature_Pair, y=Score, colour = Cluster)) +
  geom_point(size = 2) +
  xlab("Feature Pair ID") +
  ylab("Score") +
  ggtitle("Female Feature-Pair Scores") +
  ylim(c(0,2)) +
  theme_bw()


ggsave("figures/Fig1_F_Scores.pdf", width = 5, height = 5)



annot_m <- read.csv('results/annotdf_m.csv.gz')

df <- data.frame()
for (cl in unique(annot_m$cluster)) {
  df_new <- data.frame(Cluster=cl, Score=sort(annot_m$score[annot_m$cluster == cl]))
  df_new$Feature_Pair <- 1:nrow(df_new)
  df <- rbind(df, df_new)
}

# Reorder following the value of another column:
df %>%
  ggplot( aes(x=Feature_Pair, y=Score, colour = Cluster)) +
  geom_point(size = 2) +
  xlab("Feature Pair ID") +
  ylab("Score") +
  ggtitle("Male Feature-Pair Scores") +
  ylim(c(0,2)) +
  theme_bw()

ggsave("figures/Fig1_M_Scores.pdf", width = 5, height = 5)

#################################################################################

# using validation data in 4_5 to predict back to TCGA 
fmat <- matrix(data=c(41,4,0,2,0,0,19,0,0,0,0,0,11,0,0,0,1,0,18,0,0,0,3,0,40), byrow = T, ncol=5)
colnames(fmat) <- c("cluster1","cluster2","cluster3","cluster4","cluster5")
rownames(fmat) <- c("cluster1","cluster2","cluster3","cluster4","cluster5")
pheatmap(fmat, cluster_rows = F, cluster_cols = F, scale='row')

mmat <- matrix(data=c(29,0,0,1,0,
                      2,80,1,1,2,
                      1,0,24,1,0,
                      6,5,8,23,3,
                      0,0,2,0,31), byrow=T, ncol=5)
colnames(fmat) <- c("cluster1","cluster2","cluster3","cluster4","cluster5")
rownames(fmat) <- c("cluster1","cluster2","cluster3","cluster4","cluster5")
pheatmap(mmat, cluster_rows = F, cluster_cols = F, scale = 'row'  )

#################################################################################

# Plotting expression changes for important feature-pairs

plot_genepair <- function(a, b, data, res0, fileout) {
  genea <- unlist(data[,a])
  geneb <- unlist(data[,b])
  df <- data.frame(Barcode=data$Barcode, GeneA=genea, GeneB=geneb)
  df <- inner_join(df, res0, join_by('Barcode'=='SampleIDs'))
  g <- ggplot(data = df, mapping = aes(x=BestCalls, y=log(GeneA/GeneB), colour = BestCalls)) +
    geom_boxplot() + 
    geom_jitter() +
    ggtitle(paste0(a, " > ", b))
  if (!is.na(fileout)) {
    ggsave(file = fileout, plot = g, width = 5, height = 5)
  } 
  return(g)
}

load( file='results/val_array_data_with_labels.rda')
f_val_pred <- read.csv('results/females_validation_results_table.csv')

# RBP1XBMP2
genea <- "RBP1"
geneb <- "BMP2"
fileout <- 'figures/f_genepairs_RBP1_BMP2.pdf'
fileout <- NA
plot_genepair(genea,geneb,val_f,f_val_pred,fileout)

#SCN3A_X_DYNLT3
genea <- "SCN3A"
geneb <- "DYNLT3"
fileout <- 'figures/f_genepairs_SCN3A_DYNLT3.pdf'
fileout <- NA
plot_genepair(genea,geneb,val_f,f_val_pred,fileout)

#BNIP3_X_CD74
genea <- "BNIP3"
geneb <- "CD74"
fileout <- 'figures/f_genepairs_BNIP3_CD74.pdf'
fileout <- NA
plot_genepair(genea,geneb,val_f,f_val_pred,fileout)

# EPHB1_X_IGFBP2
genea <- "EPHB1"
geneb <- "IGFBP2"
fileout <- 'figures/f_genepairs_EPHB1_IGFBP2.pdf'
fileout <- NA
plot_genepair(genea,geneb,val_f,f_val_pred,fileout)

# NNMT_X_SH3GL2
genea <- "NNMT"
geneb <- "SH3GL2"
fileout <- 'figures/f_genepairs_NNMT_SH3GL2.pdf'
fileout <- NA
plot_genepair(genea,geneb,val_f,f_val_pred,fileout)

#EMP3_X_FERMT1
genea <- "EMP3"
geneb <- "FERMT1"
fileout <- 'figures/f_genepairs_val_EMP3_FERMT1.pdf'
fileout <- NA
plot_genepair(genea,geneb,val_f,res0,fileout)

#PDPN_X_DLL3
genea <- "PDPN"
geneb <- "DLL3"
fileout <- 'figures/f_genepairs_val_PDPN_DLL3.pdf'
fileout <- NA
plot_genepair(genea,geneb,val_f,f_val_pred,fileout)




### same in RNA ###

genea <- "RBP1"
geneb <- "BMP2"
fileout <- 'figures/f_genepairs_tempus_RBP1_BMP2.pdf'
fileout <- NA
plot_genepair(genea,geneb,tmp_f,f_tmp_pred,fileout)

genea <- "EPHB1"
geneb <- "IGFBP2"
fileout <- 'figures/f_genepairs_tempus_EPHB1_IGFBP2.pdf'
fileout <- NA
plot_genepair(genea,geneb,tmp_f,f_tmp_pred,fileout)

genea <- "NNMT"
geneb <- "SH3GL2"
fileout <- 'figures/f_genepairs_tempus_NNMT_SH3GL2.pdf'
fileout <- NA
plot_genepair(genea,geneb,tmp_f,f_tmp_pred,fileout)

#EMP3_X_FERMT1
genea <- "EMP3"
geneb <- "FERMT1"
fileout <- 'figures/f_genepairs_tempus_EMP3_FERMT1.pdf'
fileout <- NA
plot_genepair(genea,geneb,tmp_f,f_tmp_pred,fileout)

#PDPN_X_DLL3
genea <- "PDPN"
geneb <- "DLL3"
fileout <- 'figures/f_genepairs_tempus_PDPN_DLL3.pdf'
fileout <- NA
plot_genepair(genea,geneb,tmp_f,f_tmp_pred,fileout)


#### plot in males ####

plot_genepair_m <- function(a, b, data, fileout) {
  genea <- unlist(data[,a])
  geneb <- unlist(data[,b])
  df <- data.frame(Barcode=data$Barcode, BestCalls=data$ClusterLabel, GeneA=genea, GeneB=geneb)
  g <- ggplot(data = df, mapping = aes(x=BestCalls, y=log(GeneA/GeneB), colour = BestCalls)) +
    geom_boxplot() + 
    geom_jitter() +
    ggtitle(paste0(a, " > ", b))
  if (!is.na(fileout)) {
    ggsave(file = fileout, plot = g, width = 5, height = 5)
  } 
  return(g)
}


genea <- "RBP1"
geneb <- "BMP2"
fileout <- 'figures/m_genepairs_tcga_RBP1_BMP2.pdf'
fileout <- NA
plot_genepair(genea,geneb,dat_m,fileout)


genea <- "EMP3"
geneb <- "FERMT1"
fileout <- 'figures/m_genepairs_tcga_EMP3_FERMT1.pdf'
fileout <- NA
plot_genepair_m(genea,geneb,dat_m,fileout)




##################################################################################
# RNA vs array

bs <- intersect(dat_f$Barcode, tcga_f$SampleBarcode)
rownames(dat_f) <- dat_f$Barcode
rownames(tcga_f) <- tcga_f$SampleBarcode

gene <- "FERMT1"

df <- data.frame(RNAseq=as.numeric(unlist(tcga_f[bs, gene])),
                 Array=as.numeric(unlist(dat_f[bs, gene])),
                 ClusterLabel=dat_f[bs,"ClusterLabel"])

ggplot(data=df, aes(x=Array,y=Array/RNAseq, color=ClusterLabel)) +
  geom_point(size = 2) +
  xlab("Expression Array") +
  ylab("Array / RNA-seq") +
  ggtitle("Female gene expression ratio") +
  theme_bw()

ggsave("female_array_vs_rnaseq_Ratio.pdf", width = 5, height = 5)



gene1 <- "FERMT1"
gene2 <- "EMP3"

df <- data.frame(RNAseq1=as.numeric(unlist(tcga_f[bs, gene1])),
                 Array1=as.numeric(unlist(dat_f[bs, gene1])),
                 RNAseq2=as.numeric(unlist(tcga_f[bs, gene2])),
                 Array2=as.numeric(unlist(dat_f[bs, gene2])),
                 ClusterLabel=dat_f[bs,"ClusterLabel"])

ggplot(data=df, aes(x=Array1<Array2,y=RNAseq1<RNAseq2)) + #, #color=ClusterLabel)) +
  geom_bin2d() +
  stat_bin2d(geom = "text", aes(label = ..count..), binwidth = 1) +
  scale_fill_gradient(low = "white", high = "red") +
  xlab("Array, EMP3 > FERMT1") +
  ylab("RNA-seq, EMP3 > FERMT1") +
  ggtitle("Female gene expression ratio") +
  theme_bw()

ggsave("female_array_vs_rnaseq_RatioCount.pdf", width = 5, height = 5)







