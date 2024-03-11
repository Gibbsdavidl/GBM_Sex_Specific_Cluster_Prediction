

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("GEOquery", quietly = TRUE))
  BiocManager::install("GEOquery")

if(!require("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")

# load series and platform data from GEO
#gset13041 <- getGEO("GSE13041", GSEMatrix =TRUE, AnnotGPL=TRUE)
#gset16011 <- getGEO("GSE16011", GSEMatrix =TRUE, AnnotGPL=TRUE)
#gsetPheno <- read_delim('../data/validation_metadata/GBM_meta_analysis.txt')

#save(gset13041, gset16011, file = '../data/jive_validation_array_data.rda')

load('data/jive_validation_array_data.rda')

# The list item [[2]] is not used in the JIVE paper 

# make proper column names to match toptable 
fvarLabels(gset13041[[1]]) <- make.names(fvarLabels(gset13041[[1]]))
fvarLabels(gset13041[[2]]) <- make.names(fvarLabels(gset13041[[2]])) 
fvarLabels(gset13041[[3]]) <- make.names(fvarLabels(gset13041[[3]]))

#ex1 <- exprs(gset13041[[1]]) # expression will come from scan
#ex2 <- exprs(gset13041[[2]])
#ex3 <- exprs(gset13041[[3]])

ph1 <- phenoData(gset13041[[1]])
ph2 <- phenoData(gset13041[[2]])
ph3 <- phenoData(gset13041[[3]])

phenotype_data1 <- pData(ph1)
phenotype_data2 <- pData(ph2)
phenotype_data3 <- pData(ph3)

sym1 <- fData(gset13041[[1]])$Gene.symbol
sym2 <- fData(gset13041[[2]])$Gene.symbol
sym3 <- fData(gset13041[[3]])$Gene.symbol

id1 <- fData(gset13041[[1]])$ID
id2 <- fData(gset13041[[2]])$ID
id3 <- fData(gset13041[[3]])$ID

probeMap1 <- data.frame(ID=id1, Symbol=sym1)
probeMap2 <- data.frame(ID=id2, Symbol=sym2)
probeMap3 <- data.frame(ID=id3, Symbol=sym3)

#df1 <- as.data.frame(t(ex1))
#colnames(df1) <- sym1
#df1[['sample']] <- rownames(df1)
#df2 <- as.data.frame(t(ex2))
#colnames(df2) <- sym2
#df3 <- as.data.frame(t(ex3))
#colnames(df3) <- sym3


pdf1 <- phenotype_data1[, c('geo_accession', 'characteristics_ch1.1', 'characteristics_ch1.2',
                            'characteristics_ch1.3', 'characteristics_ch1.4',
                            'characteristics_ch1.5', 'characteristics_ch1.6')]

pdf2 <- phenotype_data2[, c('geo_accession', 'characteristics_ch1.1', 'characteristics_ch1.2',
                            'characteristics_ch1.3', 'characteristics_ch1.4',
                            'characteristics_ch1.5', 'characteristics_ch1.6')]

pdf3 <- phenotype_data3[, c('geo_accession', 'characteristics_ch1.1', 'characteristics_ch1.2',
                            'characteristics_ch1.3', 'characteristics_ch1.4',
                            'characteristics_ch1.5', 'characteristics_ch1.6')]


table(pdf3[pdf3$characteristics_ch1.6 == 'gender: F', 'characteristics_ch1.4'])


days1 <- as.numeric(
  gsub(pattern = 'tts\\(days\\): ', 
       replacement = '', 
       pdf1$characteristics_ch1.1)
)

days2 <- as.numeric(
  gsub(pattern = 'tts\\(days\\): ', 
       replacement = '', 
       pdf2$characteristics_ch1.1)
)


days3 <- as.numeric(
  gsub(pattern = 'tts\\(days\\): ', 
       replacement = '', 
       pdf3$characteristics_ch1.1)
)

#status: censoring status 1=censored, 2=dead
vitalstatus1 <- ifelse(pdf1$characteristics_ch1.2 == 'vital status: ALIVE', yes=1, no=2)
vitalstatus2 <- ifelse(pdf2$characteristics_ch1.2 == 'vital status: ALIVE', yes=1, no=2)
vitalstatus3 <- ifelse(pdf3$characteristics_ch1.2 == 'vital status: ALIVE', yes=1, no=2)

pheno1 <- data.frame('id'=pdf1$geo_accession, 'days'=days1, 'sex'=pdf1$characteristics_ch1.6, status=vitalstatus1)
pheno2 <- data.frame('id'=pdf2$geo_accession, 'days'=days2, 'sex'=pdf2$characteristics_ch1.6, status=vitalstatus2)
pheno3 <- data.frame('id'=pdf3$geo_accession, 'days'=days3, 'sex'=pdf3$characteristics_ch1.6, status=vitalstatus3)


save(pheno1,pheno2,pheno3, file='data/validation_array_data/pheno13041.rda')
save(probeMap1,probeMap2,probeMap3, file='data/validation_array_data/symbols13041.rda')

#Then we need to convert the esets to data.frames

load('../data/validation_array_data/pheno13041.rda')
load('../data/validation_array_data/symbols13041.rda')
load('../data/validation_array_data/esets/norm13041_esetList.rda')
load('../data/validation_array_data/esets/norm13041_MDA.rda')
load('../data/validation_array_data/esets/norm13041_TB.rda')

# pheno3
esetList[[1]]  # 105 samples 22283 features
norm13041_MDA  # 55 samples  22283 features
norm13041_TB   # 31 samples  22283 features
# 105+31+55 # 191

# pheno2
esetList[[2]]  # 49 samples  12625 features

#pheno1
esetList[[3]]  # 27 samples  54675 features



# first set
mat1 <- as.data.frame(esetList[[3]])
#data.frame(colnames(mat1), probeMap3$ID)[10000:10050,]  # matches
colnames(mat1) <- probeMap1$Symbol
geoids <- rownames(mat1)
geoids <- unlist(sapply(strsplit(geoids,'\\.'), function(a) a[1]))
rownames(mat1) <- geoids
rm(geoids)
save(mat1, file='../data/validation_array_data/GSE13041_mat1.rda')


# second set
mat2 <- as.data.frame(esetList[[2]])
#data.frame(colnames(mat1), probeMap3$ID)[10000:10050,]  # matches
colnames(mat2) <- probeMap2$Symbol
geoids <- rownames(mat2)
geoids <- unlist(sapply(strsplit(geoids,'\\.'), function(a) a[1]))
rownames(mat2) <- geoids
rm(geoids)
save(mat2, file='../data/validation_array_data/GSE13041_mat2.rda')


# third set
mat3_1 <- as.data.frame(esetList[[1]])
mat3_2 <- as.data.frame(norm13041_MDA)
mat3_3 <- as.data.frame(norm13041_TB)
mat3 <- rbind(mat3_1, mat3_2, mat3_3)
rm(mat3_1, mat3_2, mat3_3)
#data.frame(colnames(mat1), probeMap3$ID)[10000:10050,]  # matches
colnames(mat3) <- probeMap3$Symbol
geoids <- rownames(mat3)
geoids <- unlist(sapply(strsplit(geoids,'\\.'), function(a) a[1]))
rownames(mat3) <- geoids
rm(geoids)
save(mat3, file='../data/validation_array_data/GSE13041_mat3.rda')



library(robencla)

load('../data/validation_array_data/GSE13041_mat3.rda')

# load up the model
load('../models/females_tcga_array_step4_5.rda')

# load up the gene pairs used
load('../results/females_genepairs.rda')

# check what genes might be missing
unlist(genepairs) [! unlist(genepairs) %in% colnames(mat1) ]

# check what genes might be missing
unlist(genepairs) [! unlist(genepairs) %in% colnames(mat2) ]

# check what genes might be missing
unlist(genepairs) [! unlist(genepairs) %in% colnames(mat3) ]


### mat2 has too many missing genes ###

# cluster23 cluster213 
#   "PTCD1"     "CCT3" 

# mat3 for females
# cluster23 cluster213 cluster415 
#   "PTCD1"     "CCT3"    "CRIM1" 

PTCD1_alias <- c("PTCD1") 
sapply(PTCD1_alias, function(a) {
  any(str_detect(string = colnames(mat1), pattern=a))
})
colnames(mat1) [
  which(str_detect(string = colnames(mat1), pattern="PTCD1"))
]
which(str_detect(string = colnames(mat1), pattern="PTCD1")) # 28241 32076
colnames(mat1)[28241] <- "PTCD1" # ""ATP5J2-PTCD1///PTCD1"


PTCD1_alias <- c("PTCD1") 
sapply(PTCD1_alias, function(a) {
  any(str_detect(string = colnames(mat3), pattern=a))
})
colnames(mat3) [
  which(str_detect(string = colnames(mat3), pattern="PTCD1"))
]
which(str_detect(string = colnames(mat3), pattern="PTCD1"))
colnames(mat3)[18320] <- "PTCD1" # ""ATP5J2-PTCD1///PTCD1"



CCT3_alias <- c("CCT3") 
sapply(CCT3_alias, function(a) {
  any(str_detect(string = colnames(mat1), pattern=a))
})
colnames(mat1) [
  which(str_detect(string = colnames(mat1), pattern="CCT3"))
]
which(str_detect(string = colnames(mat1), pattern="CCT3"))
colnames(mat1)[10359] <- "CCT3" # "LOC101927137///CCT3"



CCT3_alias <- c("CCT3") 
sapply(CCT3_alias, function(a) {
  any(str_detect(string = colnames(mat3), pattern=a))
})
colnames(mat3) [
  which(str_detect(string = colnames(mat3), pattern="CCT3"))
]
which(str_detect(string = colnames(mat3), pattern="CCT3"))
colnames(mat3)[438] <- "CCT3" # "LOC101927137///CCT3"



CRIM1_alias <- c("CRIM1") 
sapply(CRIM1_alias, function(a) {
  any(str_detect(string = colnames(mat3), pattern=a))
})
colnames(mat3) [
  which(str_detect(string = colnames(mat3), pattern="CRIM1"))
]
which(str_detect(string = colnames(mat3), pattern="CRIM1"))
cor(mat3[,2079], mat3[,2080])
colnames(mat3)[2079] <- "CRIM1" # ""ATP5J2-PTCD1///PTCD1"


### build the data set with the gene pairs

expr_f1 <- mat1[, unique(unlist(genepairs))]
expr_f1[["sample"]] <- rownames(expr_f1)

expr_f3 <- mat3[, unique(unlist(genepairs))]
expr_f3[["sample"]] <- rownames(expr_f3)


### split into male and females
rownames(pheno1) <- pheno1$id
pheno1 <- pheno1[expr_f1$sample, ]
all(pheno1$id == expr_f1$sample)

rownames(pheno3) <- pheno3$id
pheno3 <- pheno3[expr_f3$sample, ]
all(pheno3$id == expr_f3$sample)

f1_idx <- which(pheno1$sex == 'gender: F')
f1_ids <- pheno1$id[pheno1$sex == 'gender: F']
expr_f1 <- expr_f1[f1_ids,]

f3_idx <- which(pheno3$sex == 'gender: F')
f3_ids <- pheno3$id[pheno3$sex == 'gender: F']
expr_f3 <- expr_f3[f3_ids,]

# save data
save(expr_f1,expr_f3, file='../data/expr_f_13041.rda')


