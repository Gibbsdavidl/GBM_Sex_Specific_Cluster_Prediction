### R code from vignette source 'SCAN.vignette.Rnw'
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("SCAN.UPC")

###################################################
### code chunk number 1: setup
###################################################
## do work in temporary directory
pwd <- getwd()  # after setting to source file ### setwd(tempdir())


###################################################
### code chunk number 2: download-geo-direct (eval = FALSE)
###################################################
## library(GEOquery)
## tmpDir = tempdir()
## library(GEOquery)
## getGEOSuppFiles("GSM555237", makeDirectory=FALSE, baseDir=tmpDir)
#  celFilePath = file.path(pwd, "GSM326661.CEL.gz")


###################################################
### code chunk number 3: download-normalize (eval = FALSE)
###################################################
library(SCAN.UPC)

library(doParallel) 
registerDoParallel(cores=6)


## First we run SCAN on the raw array files

norm13041_GBM = SCAN("E:/Work/GBM_clusters/Array Processing/GSE13041_RAW/GBM/*.CEL.gz")
norm13041_MDA = SCAN("E:/Work/GBM_clusters/Array Processing/GSE13041_RAW/MDA_UCSF/*.CEL.gz")
norm13041_TB = SCAN("E:/Work/GBM_clusters/Array Processing/GSE13041_RAW/TB/*.CEL.gz")
save(norm13041_GBM, file = 'norm13041_GBM.rda')
save(norm13041_MDA, file = 'norm13041_MDA.rda')
save(norm13041_TB, file = 'norm13041_TB.rda')

norm16011 = SCAN("E:/Work/GBM_clusters/Array Processing/GSE16011_RAW/*.CEL.gz")
save(norm16011, file = 'norm16011.rda')

dirpath = "E:/Work/GBM_clusters/Array Processing/GSE13041_RAW/GBM/"


# Then this collection is actually three different array platforms #
# so we sort them out into a list, by length of the number of features

eset1 = SCAN("E:/Work/GBM_clusters/Array Processing/GSE13041_RAW/GBM/GSM326661.CEL.gz")
fnames1 <- featureNames(eset1)
lenset <- unique(length(fnames1))
esetList <- list(eset1)

for (file in list.files(path = "E:/Work/GBM_clusters/Array Processing/GSE13041_RAW/GBM/", full.names = T)) {
  
  print(file)
  eset2 = SCAN(file)
  fnames2 <- featureNames(eset2)
  
  if (length(fnames2) %in% lenset) {
    # which part of the list is it?
    idx <- which(lenset == length(fnames2))  
    esetList[[idx]] <- combine(esetList[[idx]], eset2)
    print(esetList[[idx]])
    
  } else {
    # it's a new size of feature names
    lenset <- c(lenset, length(fnames2))
    idx <- which(lenset == length(fnames2))  
    esetList[[idx]] <- eset2
    print("NEW LENGTH!!")
    print(esetList[[idx]])
  }
}

save(esetList, file = 'norm13041_esetList.rda')




## TCGA Processing

datalist <- list()
dirlist <- c()
for (dir in list.dirs(path = "E:/Work/GBM_clusters/GBM_JIVE_ANALYSIS/data/GBM.Gene_Expression.Level_1/", full.names = T)) {
  if (dir != "E:/Work/GBM_clusters/GBM_JIVE_ANALYSIS/data/GBM.Gene_Expression.Level_1/")
  {
    print(dir)
    dirlist <- c(dirlist, dir)
    datalist[[dir]] <- SCAN(paste(dir,"*.CEL",sep='/'))
    
    }
  
}
save(datalist, dirlist, file='training_array_scan_proc.rda')


# then we combine the list of expression sets
tcga_gbm <- datalist[[1]]
for (i in 2:length(datalist)){
  tcga_gbm <- combine(tcga_gbm, datalist[[i]])
}
save(tcga_gbm, file='data/TCGA_GBM_array_scan_reprocessed.rda')


# getting the barcodes out from the data.freeze files

m <- exprs(tcga_gbm)
m_cols <- colnames(m)

library(readr)
library(stringr)
datafreeze <- read_delim('data/GBM.Gene_Expression.Level_1/data.freeze.txt',delim = '\t')

sampleids <- c()
unkn <- 1
for (mi in m_cols) {
  idx  <- which(str_detect(string = datafreeze$FILE_LOCATION_URL, pattern = mi))
  if ( length(idx) > 0 ) {
    sampleids <- c(sampleids, datafreeze$SAMPLE_BARCODE[idx])
  } else {
    print(idx)
    print(mi)
    sampleids <- c(sampleids, paste("unknown_",unkn,sep='_'))
    unkn <- unkn+1
  }
}
colnames(m) <- sampleids


#https://bioconductor.org/packages/release/data/annotation/html/pd.ht.hg.u133a.html
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL3921
#load('data/validation_array_data/symbols13041.rda')

probemap_gpl3921 <- read_delim('data/GBM.Gene_Expression.Level_1/GPL3921-25447.txt',delim = '\t', comment = '#')

sum(rownames(m) %in% probeMap3$ID)
# 22277 matches as stated in the GPL3921 page

# trying probmap 3
symbolslist3 <- c()
unkn <- 1
for (pi in rownames(m)) {
  idx <- which(probeMap3$ID == pi)
  si <- probeMap3[idx,2]
  if (si == "") {
    symbolslist3 <- c(symbolslist3, paste("unknown",as.character(unkn),sep = "_"))
    unkn <- unkn+1
  } else {
    symbolslist3 <- c(symbolslist3, si)
  }
}

# trying gpl3921 # should be the same?
symbolslist <- c()
unkn <- 1
for (pi in rownames(m)) {
  idx <- which(probemap_gpl3921$ID == pi)
  si <- probemap_gpl3921[idx,"Gene Symbol"]
  if (is.na(si)) {
    symbolslist <- c(symbolslist, paste("unknown",as.character(unkn),sep = "_"))
    unkn <- unkn+1
  } else {
    symbolslist <- c(symbolslist, si)
  }
}


### symbolslist3 was better, had newer symbol names.
### the other GPL had old definitions.

rownames(m) <- symbolslist3
tcga_gbm <- m
write.csv(m, file = 'data/tcga_gbm_array_scan_proc.csv')
save(tcga_gbm, file = 'data/tcga_gbm_array_scan_proc.rda')

###################################################
### code chunk number 4: scan-geo (eval = FALSE)
###################################################
## normalized = SCAN("GSM555237")


###################################################
### code chunk number 5: download-normalize2 (eval = FALSE)
###################################################
## normalized = SCAN(celFilePath, outFilePath="output_file.txt")


###################################################
### code chunk number 6: download-brainarray (eval = FALSE)
###################################################
## install.packages("hgu95ahsentrezgprobe_15.0.0.tar.gz", repos=NULL, type="source")


###################################################
### code chunk number 7: install-brainarray (eval = FALSE)
###################################################
## pkgName = InstallBrainArrayPackage(celFilePath, "15.0.0", "hs", "entrezg")


###################################################
### code chunk number 8: scan-brainarray (eval = FALSE)
###################################################
## normalized = SCAN(celFilePath, probeSummaryPackage=pkgName)


###################################################
### code chunk number 9: scan-twocolor-main (eval = FALSE)
###################################################
## SCAN_TwoColor("GSM1072833", "output_file.txt")


###################################################
### code chunk number 10: upc-microarray (eval = FALSE)
###################################################
## upc1 = UPC("GSM555237")
## 
## upc2 = UPC_TwoColor("GSM1072833")


###################################################
### code chunk number 11: upc-rnaseq (eval = FALSE)
###################################################
## upc3 = UPC_RNASeq("ReadCounts.txt", "Annotation.txt")


###################################################
### code chunk number 12: scan-parallel (eval = FALSE)
###################################################
## library(doParallel)
## registerDoParallel(cores=2)
## result = SCAN("GSE22309")
