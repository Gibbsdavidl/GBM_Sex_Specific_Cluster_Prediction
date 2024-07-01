

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("GEOquery", quietly = TRUE))
  BiocManager::install("GEOquery")

if(!require("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")

if(!require("SCAN.UPC", quietly = T)) 
  BiocManager::install("SCAN.UPC")

if(!require('pd.hg.u133.plus.2', quietly = T)) 
  BiocManager::install('pd.hg.u133.plus.2')

if(!require('hgu133plus2.db', quietly = T)) 
  BiocManager::install("hgu133plus2.db")

library(stringr)
library(readr)
library(doParallel) 
registerDoParallel(cores=5)

GSE108474_1 = SCAN("C:/Users/dgibbs/ISB_Work/GBM_Sex/REMBRANDT/array_cels/set1/*.CEL.gz")
GSE108474_2 = SCAN("C:/Users/dgibbs/ISB_Work/GBM_Sex/REMBRANDT/array_cels/set2/*.CEL.gz")
GSE108474_3 = SCAN("C:/Users/dgibbs/ISB_Work/GBM_Sex/REMBRANDT/array_cels/set3/*.CEL.gz")

GSE108474_4 = SCAN("C:/Users/dgibbs/ISB_Work/GBM_Sex/REMBRANDT/array_cels/set4/*.CEL.gz")
GSE108474_5 = SCAN("C:/Users/dgibbs/ISB_Work/GBM_Sex/REMBRANDT/array_cels/set5/*.CEL.gz")
GSE108474_6 = SCAN("C:/Users/dgibbs/ISB_Work/GBM_Sex/REMBRANDT/array_cels/set6/*.CEL.gz")

GSE108474_7 = SCAN("C:/Users/dgibbs/ISB_Work/GBM_Sex/REMBRANDT/array_cels/set7/*.CEL.gz")
GSE108474_8 = SCAN("C:/Users/dgibbs/ISB_Work/GBM_Sex/REMBRANDT/array_cels/set8/*.CEL.gz")
GSE108474_9 = SCAN("C:/Users/dgibbs/ISB_Work/GBM_Sex/REMBRANDT/array_cels/set9/*.CEL.gz")
GSE108474_10 = SCAN("C:/Users/dgibbs/ISB_Work/GBM_Sex/REMBRANDT/array_cels/set10/*.CEL.gz")


# save as a list
GSE108474_list <- list()
GSE108474_list[[1]] <- GSE108474_1
GSE108474_list[[2]] <- GSE108474_2
GSE108474_list[[3]] <- GSE108474_3
GSE108474_list[[4]] <- GSE108474_4
GSE108474_list[[5]] <- GSE108474_5
GSE108474_list[[6]] <- GSE108474_6
save(GSE108474_list, file = 'C:/Users/dgibbs/ISB_Work/GBM_Sex/REMBRANDT/GSE_list.rda')
rm(GSE108474_1,GSE108474_2,GSE108474_3,GSE108474_4,GSE108474_5,GSE108474_6)

load('C:/Users/dgibbs/ISB_Work/GBM_Sex/REMBRANDT/GSE_list.rda')
GSE108474_list1 <- GSE108474_list
load('C:/Users/dgibbs/ISB_Work/GBM_Sex/REMBRANDT/GSE_list2.rda')
GSE108474_list2 <- GSE108474_list

# then combine into a single.
GSE108474 <- combine(GSE108474_list1[[1]], GSE108474_list1[[2]])
GSE108474 <- combine(GSE108474, GSE108474_list1[[3]])
GSE108474 <- combine(GSE108474, GSE108474_list1[[4]])
GSE108474 <- combine(GSE108474, GSE108474_list1[[5]])
GSE108474 <- combine(GSE108474, GSE108474_list1[[6]])
GSE108474 <- combine(GSE108474, GSE108474_list2[[7]])
GSE108474 <- combine(GSE108474, GSE108474_list2[[8]])
GSE108474 <- combine(GSE108474, GSE108474_list2[[9]])
GSE108474 <- combine(GSE108474, GSE108474_list2[[10]])
rm(GSE108474_list1, GSE108474_list2, GSE108474_list)
gc()


# pull in resource tables
clin <- read_delim('C:/Users/dgibbs/ISB_Work/GBM_Sex/REMBRANDT/GSE108474_REMBRANDT_clinical.data.txt')
bios <- read_delim('C:/Users/dgibbs/ISB_Work/GBM_Sex/REMBRANDT/GSE108474_REMBRANDT_biospecimen_mapping_GEO.txt')
arry <- read_delim('C:/Users/dgibbs/ISB_Work/GBM_Sex/REMBRANDT/GPL570-55999.txt',comment='#')
acce <- read_delim('C:/Users/dgibbs/ISB_Work/GBM_Sex/REMBRANDT/GSE108474_sample.tsv')

# create the expression matrix
X <- exprs(GSE108474)
dim(X)

# the stat used to filter out probes
Xvar <- apply(X, 1, var, na.rm=T)

# https://rdrr.io/github/Bioconductor/genefilter/man/findLargest.html
# https://support.bioconductor.org/p/23397/
fidx <- findLargest(gN = rownames(X), testStat = Xvar, data = "hgu133plus2")

length(fidx)   # 20857

# filter the probes
Xfilt <- X[fidx,]

# check the ID mapping
x <- hgu133plus2ENTREZID
geneids <- sapply(rownames(Xfilt), function(a) x[[a]])
sum(duplicated(geneids))
# [1] 0

xs <- hgu133plus2SYMBOL
geneids <- sapply(rownames(Xfilt), function(a) xs[[a]])
sum(duplicated(geneids))
# [1] 0

exprmat <- as.data.frame(t(Xfilt))
colnames(exprmat) <- geneids

rm(GSE108474, X, Xfilt)
gc()

# getting the IDs in format to match to clinical table
split_sampids <- str_split_fixed(rownames(exprmat),'_',n = 3)
sampleids <- as.character(split_sampids[,2])
geoids <- as.character((split_sampids[,1]))

# get the index into the biospecimen table.
# and then get the subject ID
cidx <- match(geoids, acce$Accession)
all(geoids == acce$Accession[cidx])
#TRUE

# get the subject IDs for each accession
subjectids <- acce$Title[cidx]

# lay in the row names
rownames(exprmat) <- subjectids

clinidx <- match(rownames(exprmat), clin$SUBJECT_ID)

# check the ordering is correct
data.frame(rownames(exprmat), clin$SUBJECT_ID[clinidx])
all(rownames(exprmat) == clin$SUBJECT_ID[clinidx], na.rm = T)
#[1] TRUE

# bring in the clin factors
exprmat$SubjectID <- clin$SUBJECT_ID[clinidx]
exprmat$Sex       <- clin$GENDER[clinidx]
exprmat$Censored  <- clin$EVENT_OS[clinidx]
exprmat$Survival  <- clin$OVERALL_SURVIVAL_MONTHS[clinidx]

exprmat[,20855:20861]

# DUPLICATES
naidx <- which(is.na(exprmat$SubjectID))
#exprmat[naidx, 20851:20861]
exprmat <- exprmat[-naidx, ]

# remove normal samples and those
naidx <- which(is.na(exprmat$Censored))
#exprmat[naidx, 20851:20861]
exprmat <- exprmat[-naidx, ]

# remove samples without sex
naidx <- which(is.na(exprmat$Sex))
#exprmat[naidx, 20851:20861]
exprmat <- exprmat[-naidx, ]

exprmat_F <- exprmat[exprmat$Sex == 'FEMALE',]
exprmat_M <- exprmat[exprmat$Sex == 'MALE',]

write_csv(exprmat_F, file = 'data/Array Data/validation_array_gse108474_F.csv')
write_csv(exprmat_M, file = 'data/Array Data/validation_array_gse108474_M.csv')









