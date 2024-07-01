
library(tidyverse)
library(robencla)
library(survival)
library(survminer)

load('../data/F_processed_data.rda')
load('../data/f_genepairs.rda')

fgenepairs[[1]]

# first pair in cluster 3 is 
#"EMP3"    "FERMT1"

dat_f$Barcode12 <- str_sub(dat_f$Barcode, 1, 12)




df <- inner_join(dat_f[,c("Barcode12","EMP3","FERMT1")], 
                 tcga_f[,c("Barcode","EMP3","FERMT1")],
                 by=join_by("Barcode12"=="Barcode"))

head(df)

qplot(x=df$EMP3.x, y=df$EMP3.y)


dat_df <- as.data.frame(dat_f)

res0 <- sapply(1:ncol(dat_df), FUN = function(i)  
    if (is.numeric(dat_df[,i])) {
      cor(dat_df$EMP3, dat_df[,i], use = "pairwise.complete.obs")
    } else {NA}
  )

plot(sort(res0))

sort(res0,decreasing = T)[1:10]
# 0.628

colnames(dat_df)[which(res0 > 0.60)]

#[1] "PDPN"   "VAMP5"  "CHI3L1" "EMP3"   "CAVIN1" "EFEMP2" "ANXA1"  "ANXA5"  "HRH1"  
#[10] "ILK"    "LGALS1" "PTX3"   "TIMP1"  "TRIP6" 
#[1]  665  779 1065 2083 3901 4263 4276 4330 4492 4739 5067 8221 9598 9706

cor(dat_f$EMP3, dat_f$PDPN)
#>   res0[665]
#[1] 0.6545142

dat_df <- as.data.frame(dat_f)

gs <- c("PDPN", "VAMP5", "CHI3L1", "EMP3",  "CAVIN1", "EFEMP2", "ANXA1",  
        "ANXA5",  "HRH1",  "ILK", "LGALS1", "PTX3", "TIMP1",  "TRIP6" )

gs <- colnames(dat_df)[which(res0 > 0.60)]


medsum <- function(x) {
  med <- median(x)
  bot <- median(abs(x - med))
  return(med/bot)
}


xdat1 <- apply(dat_df[,gs], 1, sum) / 16
ydat1 <- apply(tcga_f[,gs], 1, sum) / 16

xdat2 <- apply(dat_df[,gs], 1, medsum) 
ydat2 <- apply(tcga_f[,gs], 1, medsum) 


dat_df$gs_avg <- xdat2
tcga_f$gs_avg <- ydat2

df <- inner_join(dat_df[,c("Barcode12","EMP3","FERMT1", "gs_avg")], 
                 tcga_f[,c("Barcode","EMP3","FERMT1", "gs_avg")],
                 by=join_by("Barcode12"=="Barcode"))


qplot(x=df$gs_avg.x, y=df$gs_avg.y)


qplot(x=df$EMP3.x, y=df$EMP3.y)



res1 <- sapply(1:ncol(dat_df), FUN = function(i)  
  if (is.numeric(dat_df[,i])) {
    cor(dat_df$FERMT1, dat_df[,i], use = "pairwise.complete.obs")
  } else {NA}
)

plot(sort(res1))

sort(res1,decreasing = T)[1:10]
# 0.628

gs <- colnames(dat_df)[which(res1 > 0.56)]
gs

xdat <- apply(dat_df[,gs], 1, sum) / length(gs)
ydat <- apply(tcga_f[,gs], 1, sum) / length(gs)

dat_df$gs_avg <- xdat
tcga_f$gs_avg <- ydat

df <- inner_join(dat_df[,c("Barcode12","EMP3","FERMT1", "gs_avg")], 
                 tcga_f[,c("Barcode","EMP3","FERMT1", "gs_avg")],
                 by=join_by("Barcode12"=="Barcode"))


qplot(x=df$FERMT1.x, y=df$FERMT1.y)

qplot(x=df$gs_avg.x, y=df$gs_avg.y)





