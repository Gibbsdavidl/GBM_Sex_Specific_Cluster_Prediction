library(ggplot2)
library(pheatmap)
library(forcats)

load('results/results_median_min/females_genepairs.rda')
gp1 <- genepairs
load('results/females_genepairs_val.rda')
gp2 <- genepairs
load('results/females_genepairs_rna.rda')
gp3 <- genepairs
load('results/females_genepairs_refined.rda')
gp4 <- genepairs


load("E:/Work/GBM_clusters/GBM_JIVE_ANALYSIS/data/F_processed_data.rda")

##################################################

# figure 1
annot_f <- read.csv('results/annotdf_f_valdiff.csv.gz')

df <- data.frame()
for (cl in unique(annot_f$cluster)) {
  df_new <- data.frame(Cluster=cl, Score=sort(annot_f$score[annot_f$cluster == cl]))
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
  theme_bw()


ggsave("figures/Fig1_F_Scores.pdf", width = 5, height = 5)



annot_m <- read.csv('results/annotdf_m_valdiff.csv.gz')

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
  theme_bw()

ggsave("figures/Fig1_M_Scores.pdf", width = 5, height = 5)

#################################################################################

# using validation data in 4_5 to predict back to TCGA 
fmat <- matrix(data=c(41,4,0,2,0,0,19,0,0,0,0,0,11,0,0,0,1,0,18,0,0,0,3,0,40), byrow = T, ncol=5)
colnames(fmat) <- c("cluster1","cluster2","cluster3","cluster4","cluster5")
rownames(fmat) <- c("cluster1","cluster2","cluster3","cluster4","cluster5")
pheatmap(fmat, cluster_rows = F, cluster_cols = F, scale='row')

mmat <- matrix(data=c(30, 13, 6, 171), byrow=T, ncol=2)
colnames(mmat) <- c("cluster3","not_cluster3")
rownames(mmat) <- c("cluster3","not_cluster3")
pheatmap(mmat, cluster_rows = F, cluster_cols = F, scale = 'row'  )
