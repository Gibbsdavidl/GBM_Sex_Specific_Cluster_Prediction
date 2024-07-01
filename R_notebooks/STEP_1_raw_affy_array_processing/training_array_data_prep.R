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

## Data was downloaded from the GDC legacy site under publications.


###################################################
### code chunk number 3: download-normalize (eval = FALSE)
###################################################
library(SCAN.UPC)
library(doParallel) 
library(readr)
library(stringr)

if(!require("hgu133a.db", quietly = TRUE))
  BiocManager::install("hgu133a.db")

registerDoParallel(cores=6)


# where the raw affy cel files are located.
cel_file_path <- "E:/Work/GBM_clusters/GBM_JIVE_ANALYSIS/data/GBM.Gene_Expression.Level_1/"
# the raw dir is a collection of dirs.
dirlist <- c()
# each sample will be an element in the list
datalist <- list()



# for each dir in the raw affy cel dir...
for (dir in list.dirs(path = cel_file_path, full.names = T)) {
  if (dir != cel_file_path)
  {
    print(dir)
    dirlist <- c(dirlist, dir)
    datalist[[dir]] <- SCAN(paste(dir,"*.CEL",sep='/'))
  }
}
# then we'll save the list of processed data.
save(datalist, dirlist, file='training_array_scan_proc.rda')



# then we'll combine the list of data sets
tcga_gbm <- datalist[[1]]
num_samples <- dim(datalist[[1]])[2]
for (i in 2:length(datalist)){
  print(dim(datalist[[i]]))
  num_samples <- num_samples + dim(datalist[[i]])[2]
  tcga_gbm <- combine(tcga_gbm, datalist[[i]])
}
print(num_samples)   # 537 + 23 = 560 samples all with 22277 features.
print(dim(tcga_gbm)) # 560 samples
save(tcga_gbm, file='data/TCGA_GBM_array_scan_reprocessed.rda')
#load('data/Array Data/TCGA_GBM_array_scan_reprocessed.rda')



load(file='data/Array Data/TCGA_GBM_array_scan_reprocessed.rda')
# then extract the expression matrix
m <- exprs(tcga_gbm) # matrix of intensities
m_cols <- colnames(m)


# links the file names to the GSM192382 IDs that are used in phenotype data
datafreeze <- read_delim('data/Array Data/validation_array_resources//data.freeze.txt',delim = '\t')
sampleids <- c()
# number of non-matched files
unkn <- 1
for (mi in m_cols) {
  idx  <- which(str_detect(string = datafreeze$FILE_LOCATION_URL, pattern = mi))
  # if there's a single match
  if ( length(idx) == 1 ) {
    sampleids <- c(sampleids, datafreeze$ALIQUOT_BARCODE[idx])
  # but sometimes there's two rows that match
  } else if ( length(idx) == 2 ){
    print("DOUBLE MATCH")
    id_list <- datafreeze$ALIQUOT_BARCODE[idx]
    if (id_list[1] == id_list[2]) {
      sampleids <- c(sampleids, id_list[1])
    } else {
      print("DOUBLE IDS DONT MATCH")
    }
  # otherwise 
  } else {
    print(idx)
    print(mi)
    sampleids <- c(sampleids, paste("unknown_",unkn,sep='_'))
    unkn <- unkn+1
  }
}


# spot check:
{
  i <- 410
  mi <- m_cols[i]
  idx  <- which(str_detect(string = datafreeze$FILE_LOCATION_URL, pattern = mi))
  print(datafreeze[idx,])
  print(sampleids[i])
}

#https://bioconductor.org/packages/release/data/annotation/html/pd.ht.hg.u133a.html
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL3921
#load('data/validation_array_data/symbols13041.rda')

dim(m)

# the stat used to filter out probes
Xvar <- apply(m, 1, var, na.rm=T)

# https://rdrr.io/github/Bioconductor/genefilter/man/findLargest.html
# https://support.bioconductor.org/p/23397/
fidx <- genefilter::findLargest(gN = rownames(m), testStat = Xvar, data = "hgu133a")


# filter the probes
Xfilt <- m[fidx,]

xs <- hgu133aSYMBOL
geneids <- sapply(rownames(Xfilt), function(a) xs[[a]])
sum(duplicated(geneids))

# bring in the genes
exprmat <- as.data.frame(t(Xfilt))
colnames(exprmat) <- geneids


### symbolslist3 was better than gpl3921, had newer symbol names.
### the other GPL had old definitions.

barcodeclusters <- read_csv('data/barcode_cluster_labels.csv')

# taking the aliquots that are not listed as "non-canonical"
x <- datafreeze[ (datafreeze$SAMPLE_BARCODE == "TCGA-06-0137-01A") & (datafreeze$DATATYPE == 'Expression-Genes'),] ## portioning study
# data error with this aliquot
x <- datafreeze[ (datafreeze$SAMPLE_BARCODE == "TCGA-06-0138-01A") & (datafreeze$DATATYPE == 'Expression-Genes'),] ## portioning study
# TCGA-06-0138-01A-02R-0233-01
x <- datafreeze[ (datafreeze$SAMPLE_BARCODE == "TCGA-06-0145-01A") & (datafreeze$DATATYPE == 'Expression-Genes'),] ## portioning study
# TCGA-06-0145-01A-01R-0219-01 
x <- datafreeze[ (datafreeze$SAMPLE_BARCODE ==  "TCGA-06-0176-01A") & (datafreeze$DATATYPE == 'Expression-Genes'),] ## portioning study
# TCGA-06-0176-01A-02R-0233-01
x <- datafreeze[ (datafreeze$SAMPLE_BARCODE ==  "TCGA-06-0154-01A") & (datafreeze$DATATYPE == 'Expression-Genes'),] ## portioning study
# TCGA-06-0154-01A-02R-0233-01
x <- datafreeze[ (datafreeze$SAMPLE_BARCODE ==  "TCGA-06-0168-01A") & (datafreeze$DATATYPE == 'Expression-Genes'),] ## portioning study
# TCGA-06-0168-01A-01R-0233-01
x <- datafreeze[ (datafreeze$SAMPLE_BARCODE ==  "TCGA-06-0156-01A") & (datafreeze$DATATYPE == 'Expression-Genes'),] ## portioning study
# TCGA-06-0156-01A-01R-0233-01
x <- datafreeze[ (datafreeze$SAMPLE_BARCODE ==  "TCGA-06-0148-01A") & (datafreeze$DATATYPE == 'Expression-Genes'),] ## portioning study
# TCGA-06-0148-01A-01R-0219-01

aliquot_list <- c('TCGA-06-0137-01A-01R-0219-01', 'TCGA-06-0137-01A-02R-0219-01','TCGA-06-0137-01A-03R-0219-01',
                  'TCGA-06-0138-01A-01R-0233-01',
                  'TCGA-06-0145-01A-02R-0219-01','TCGA-06-0145-01A-03R-0219-01','TCGA-06-0145-01A-04R-0219-01', 
                  'TCGA-06-0145-01A-05R-0219-01','TCGA-06-0145-01A-06R-0219-01',
                  'TCGA-06-0176-01A-03R-0233-01',
                  'TCGA-06-0154-01A-03R-0233-01',
                  'TCGA-06-0168-01A-02R-0233-01', 
                  'TCGA-06-0156-01A-02R-0233-01','TCGA-06-0156-01A-03R-0233-01',
                  'TCGA-06-0148-01A-03R-0233-01','TCGA-06-0148-01A-04R-0233-01','TCGA-06-0148-01A-05R-0233-01')
# which ones to remove
idx <- which(!sampleids %in% aliquot_list)

# subset to matching barcodes to JIVE
shortbarcodes <- (str_sub(sampleids[idx],start = 1,end=15))

which(duplicated(shortbarcodes))  # NONE!!!!!!!!!!!!

# subset the rows
exprmat <- exprmat[idx,]
rownames(exprmat) <- shortbarcodes

# get the barcodes that were in the JIVE study
jdx <- (rownames(exprmat) %in% barcodeclusters$Barcode)
sum(jdx)

barcodeclusters[!(barcodeclusters$Barcode %in% rownames(tcga_gbm)),]
# missing
# Barcode         ClusterLabel Sex  
#<chr>           <chr>        <chr>
#  1 TCGA-06-0216-01 cluster4     F    
x <- datafreeze[ (datafreeze$PATIENT_BARCODE == 'TCGA-06-0216') & (!is.na(datafreeze$DATA_LEVEL)) ,]
## gene expression not in the data freeze

# merge in the cluster labels
library(dplyr)
exprmat['Barcode'] <- rownames(exprmat)
jive_train_array <- inner_join(x = exprmat, y= barcodeclusters, by=c('Barcode' = 'Barcode' ))

# split by sex
jive_train_array_M <- jive_train_array[jive_train_array$Sex == 'M',]
jive_train_array_F <- jive_train_array[jive_train_array$Sex == 'F',]

write.csv(jive_train_array_M, file = 'data/jive_training_array_data_M.csv')
write.csv(jive_train_array_F, file = 'data/jive_training_array_data_F.csv')


