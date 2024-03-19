
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("GEOquery", quietly = TRUE))
  BiocManager::install("GEOquery")

if(!require("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")

library(readr)
library(stringr)
library(dplyr)

# the validation array data from validation array data prep.R
load('data/Array Data/jive_validation_array_data.rda')

# doesn't have gender!
meta16011 <- read_csv('data/Array Data/validation_array_data/GSE16011_metadata.csv')

# array metadata https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL8542
array16011 <- read_delim('data/Array Data/validation_array_data/GPL8542_probe_map_entrez_id.txt',delim = '\t',comment = '#')

lapply(gset16011, dim)
#Features  Samples 
#   17527      284 
fvarLabels(gset16011[[1]]) <- make.names(fvarLabels(gset16011[[1]]))

ph1 <- phenoData(gset16011[[1]])
phenotype_data1 <- pData(ph1)
#sym1 <- fData(gset16011[[1]])$
id1 <- fData(gset16011[[1]])$ID
#probeMap1 <- data.frame(ID=id1, Symbol=sym1)
pdf1 <- phenotype_data1[, c('geo_accession', 'characteristics_ch1.1', 'characteristics_ch1.2',
                            'gender:ch1','histology:ch1','tissue:ch1')]


table(meta16011$grade)
#G1  G2  G3  G4 
# 8  24  83 154 

table(pdf1$`histology:ch1`)
#A (grade II)  A (grade III)  control GBM (grade IV)  OA (grade II) OA (grade III)  OD (grade II) OD (grade III)   PA (grade I) 
#      13             16         8            159              3             25              8             44              8 

# first set
ex1 <- exprs(gset16011[[1]])
mat1 <- as.data.frame(ex1)
mat1[1:5,1:5]
sum(rownames(mat1) %in% array16011$ID)
which(rownames(mat1) == '1_at')
##
idx <- match(table=array16011$ID, x=rownames(mat1))
all(rownames(mat1) == array16011$ID[idx])
##
symbols <- mapIds(org.Hs.eg.db,
                 keys=as.character(array16011$ORF)[idx],
                 column="SYMBOL",
                 keytype="ENTREZID",
                 multiVals="first")
##
# some of the symbols are NAs
sum(duplicated(symbols))
#[1] 331
sum(is.na(symbols))
#[1] 332

# adding the symbols and transposing
ndx <- which(is.na(symbols))
for (ni in ndx) {
  symbols[ni] <- paste0('NA_',ni)
}
rownames(mat1) <- symbols
mat1 <- as.data.frame(t(mat1))

# add in the sample IDs
mat1$SampleID <- rownames(mat1)
mat1[1:5,17525:17528]

# extract out the sample id without the .CEL
meta16011_sampleid <- str_sub(meta16011$name, start=1, end=9)
meta16011$geo_accession <- meta16011_sampleid

# but not all samples are in the meta16011 table, because some are controls etc
midx <- match(table=mat1$SampleID, x=meta16011_sampleid)
data.frame(meta16011_sampleid, mat1$SampleID[midx])

# join the two phenotype tables together on geo_accession
pheno <- inner_join(meta16011, pdf1)

all(pheno$geo_accession == meta16011$geo_accession)
# [1] TRUE

# now we match the meta16011 table
mat2 <- mat1[midx,]

all(pheno$geo_accession == mat2$SampleID)
#[1] TRUE

mat2$Sex <- pheno$`gender:ch1`
mat2$Survival <- pheno$survival_months
mat2$Censored <- pheno$censored
mat2[1:5,17525:17531]



# GCOM1 alias of POLR2M (GRINL1A)
# GCOM1, MYZAP-POLR2M Combined Locus
i <- which(str_detect(colnames(mat2), "GCOM1"))
#[1] 2395
colnames(mat2)[i] <- 'POLR2M'

# and we can separte by sex
mat_f <- mat2[mat2$Sex == 'Female',]
mat_m <- mat2[mat2$Sex == 'Male',]

# can just try it for now....
write.csv(mat_f, 'data/Array Data/jive_validation_16011_F.csv')
write.csv(mat_m, 'data/Array Data/jive_validation_16011_M.csv')






