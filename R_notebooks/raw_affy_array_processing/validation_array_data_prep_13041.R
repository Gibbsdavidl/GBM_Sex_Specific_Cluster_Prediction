

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("GEOquery", quietly = TRUE))
  BiocManager::install("GEOquery")

if(!require("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")

if(!require("hgu133a.db", quietly = TRUE))
  BiocManager::install("hgu133a.db")


library(readr)

#From the JIVE paper
#GSE13041 consists of 174 GBM
#patients (105 males and 69 females) profiled on Affymetrix human genome U133A and U133
#plus 2.0 (samples profiled on Affymetrix human genome U95Av2 were excluded because of the
#          small number of overlapping genes).

load('data/Array Data/jive_validation_array_data.rda')

meta13041_3 <- read_delim('data/Array Data/validation_array_resources/GPL96-57554.txt',delim='\t',comment='#')

##                                                 ##   
# The list item [[2]] is not used in the JIVE paper #
##                                                 ##


# make proper column names to match toptable 
fvarLabels(gset13041[[1]]) <- make.names(fvarLabels(gset13041[[1]]))
#fvarLabels(gset13041[[2]]) <- make.names(fvarLabels(gset13041[[2]])) 
fvarLabels(gset13041[[3]]) <- make.names(fvarLabels(gset13041[[3]]))

ex1 <- exprs(gset13041[[1]]) # [1] 54675    27
#ex2 <- exprs(gset13041[[2]])
ex3 <- exprs(gset13041[[3]]) # [1] 22283   191

phenotype_data1 <- pData(  phenoData(gset13041[[1]])  )
#phenotype_data2 <- pData(ph2)
phenotype_data3 <- pData(  phenoData(gset13041[[3]])  )


pdf1 <- phenotype_data1[, c('geo_accession', 'characteristics_ch1.1', 'characteristics_ch1.2',
                            'characteristics_ch1.3', 'characteristics_ch1.4',
                            'characteristics_ch1.5', 'characteristics_ch1.6')]

#pdf2 <- phenotype_data2[, c('geo_accession', 'characteristics_ch1.1', 'characteristics_ch1.2',
#                            'characteristics_ch1.3', 'characteristics_ch1.4',
#                            'characteristics_ch1.5', 'characteristics_ch1.6')]

pdf3 <- phenotype_data3[, c('geo_accession', 'characteristics_ch1.1', 'characteristics_ch1.2',
                            'characteristics_ch1.3', 'characteristics_ch1.4',
                            'characteristics_ch1.5', 'characteristics_ch1.6')]


table(pdf3[pdf3$characteristics_ch1.6 == 'gender: F', 'characteristics_ch1.4'])


days1 <- as.numeric(
  gsub(pattern = 'tts\\(days\\): ', 
       replacement = '', 
       pdf1$characteristics_ch1.1)
)

#days2 <- as.numeric(
#  gsub(pattern = 'tts\\(days\\): ', 
#       replacement = '', 
#       pdf2$characteristics_ch1.1)
#)


days3 <- as.numeric(
  gsub(pattern = 'tts\\(days\\): ', 
       replacement = '', 
       pdf3$characteristics_ch1.1)
)

#status: censoring status 1=censored, 2=dead
vitalstatus1 <- ifelse(pdf1$characteristics_ch1.2 == 'vital status: ALIVE', yes=1, no=2)
#vitalstatus2 <- ifelse(pdf2$characteristics_ch1.2 == 'vital status: ALIVE', yes=1, no=2)
vitalstatus3 <- ifelse(pdf3$characteristics_ch1.2 == 'vital status: ALIVE', yes=1, no=2)

pheno1 <- data.frame('id'=pdf1$geo_accession, 'days'=days1, 'sex'=pdf1$characteristics_ch1.6, status=vitalstatus1)
#pheno2 <- data.frame('id'=pdf2$geo_accession, 'days'=days2, 'sex'=pdf2$characteristics_ch1.6, status=vitalstatus2)
pheno3 <- data.frame('id'=pdf3$geo_accession, 'days'=days3, 'sex'=pdf3$characteristics_ch1.6, status=vitalstatus3)


#save(pheno1,pheno3, file='data/Array Data/validation_array_data/pheno13041.rda')
#save(probeMap1,probeMap3, file='data/Array Data/validation_array_data/symbols13041.rda')

#Then we need to convert the esets to data.frames

#load('data/Array Data/validation_array_resources/')
#load('data/Array Data/validation_array_data/symbols13041.rda')
#load('data/Array Data/validation_array_data/esets/norm13041_esetList.rda')
#load('data/Array Data/validation_array_data/esets/norm13041_MDA.rda')
#load('data/Array Data/validation_array_data/esets/norm13041_TB.rda')

dim(ex3)
#[1] 22283   191

ex3[1:5,1:5]


# the stat used to filter out probes
Xvar <- apply(ex3, 1, var, na.rm=T)

# https://rdrr.io/github/Bioconductor/genefilter/man/findLargest.html
# https://support.bioconductor.org/p/23397/
fidx <- genefilter::findLargest(gN = rownames(ex3), testStat = Xvar, data = "hgu133a")


# filter the probes
Xfilt <- ex3[fidx,]

# check the ID mapping
x <- hgu133aENTREZID
geneids <- sapply(rownames(Xfilt), function(a) x[[a]])
sum(duplicated(geneids))
# [1] 0

xs <- hgu133aSYMBOL
geneids <- sapply(rownames(Xfilt), function(a) xs[[a]])
sum(duplicated(geneids))
# [1] 0

exprmat <- as.data.frame(t(Xfilt))
colnames(exprmat) <- geneids
#[1] 0

exprmat <- as.data.frame(t(Xfilt))
colnames(exprmat) <- geneids


# going to take the probe with the highest median measure,
# after looking at ACTN2 and CSH1, where just a couple probes 
# have expression and most don't

############ 
# Now to add in the phenotype

idx <- match(rownames(exprmat), pheno3$id)
all(rownames(exprmat) == pheno3$id[idx])

exprmat$SampleID <- pheno3$id[idx]
exprmat$Sex      <- pheno3$sex[idx]
exprmat$Survival <- pheno3$days[idx]
exprmat$Censored <- pheno3$status[idx]

exprmat_F <- exprmat[exprmat$Sex == 'gender: F',]
exprmat_M <- exprmat[exprmat$Sex == 'gender: M',]

write_csv(exprmat_F, file='data/Array Data/jive_validation_13041_ex3_F_v2.csv')
write_csv(exprmat_M, file='data/Array Data/jive_validation_13041_ex3_M_v2.csv')


