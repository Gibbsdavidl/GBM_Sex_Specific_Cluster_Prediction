

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

#From the JIVE paper
#GSE13041 consists of 174 GBM
#patients (105 males and 69 females) profiled on Affymetrix human genome U133A and U133
#plus 2.0 (samples profiled on Affymetrix human genome U95Av2 were excluded because of the
#          small number of overlapping genes).

load('data/Array Data/jive_validation_array_data.rda')

meta13041_1 <- read_delim('data/Array Data/validation_array_data/GPL96-57554.txt',delim='\t',comment='#')

##
# The list item [[2]] is not used in the JIVE paper ##
##

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


save(pheno1,pheno3, file='data/Array Data/validation_array_data/pheno13041.rda')
save(probeMap1,probeMap3, file='data/Array Data/validation_array_data/symbols13041.rda')

#Then we need to convert the esets to data.frames

#load('data/Array Data/validation_array_data/pheno13041.rda')
#load('data/Array Data/validation_array_data/symbols13041.rda')
#load('data/Array Data/validation_array_data/esets/norm13041_esetList.rda')
#load('data/Array Data/validation_array_data/esets/norm13041_MDA.rda')
#load('data/Array Data/validation_array_data/esets/norm13041_TB.rda')

dim(ex3)
#[1] 22283   191

# map IDs to symbols
gidx <- match(rownames(ex3), meta13041_1$ID)
all(rownames(ex3) == meta13041_1$ID[gidx])
#[1] TRUE

ex3_entrez <-  meta13041_1$ENTREZ_GENE_ID[gidx]

sum(duplicated(meta13041_1$ENTREZ_GENE_ID))
#[1] 9038

ex3_entrez_take1 <- as.character(
  sapply(ex3_entrez, function(x) str_split(x, pattern=' /// ')[[1]][1] )
)

ex3_symbols <- mapIds(org.Hs.eg.db,
                  keys=as.character(ex3_entrez_take1),
                  column="SYMBOL",
                  keytype="ENTREZID",
                  multiVals="first")


##
# many genes have more than one probe mapped to it
sum(duplicated(ex3_symbols))
#[1] 9038
sum(is.na(ex3_symbols))
#[1] 75
sum(ex3_symbols == "NULL")
#[1] 1310

# start with removing the NAs
ndx <- which(is.na(ex3_symbols))
for (ni in ndx) {
  ex3_symbols[ni] <- paste0('NA_',ni)
}

# then removing the NULLs
ndx <- which(ex3_symbols == "NULL")
for (ni in ndx) {
  ex3_symbols[ni] <- paste0('NULL_',ni)
}


# going to take the probe with the highest median measure,
# after looking at ACTN2 and CSH1, where just a couple probes 
# have expression and most don't

ex3_mat <- as.data.frame(t(ex3))

# get the unique list of genes that have more than one probe
gene_dups <- unlist(unique(ex3_symbols[duplicated(ex3_symbols)]))

# for each gene
for (gi in gene_dups) {
  # get the index for this gene
  cidx <- which(ex3_symbols == gi)
  # compute the variance for each probe
  mexpr <- apply(ex3_mat[,cidx], 2, var, na.rm=T)
  # determine which one to take
  idx <- cidx[mexpr == max(mexpr)]
  if (length(idx) > 1) {
    print("MORE THAN ONE MAX!!")
    idx <- idx[1]
  }
  # collection probes that are less 
  jdx <- cidx[mexpr < max(mexpr,na.rm = T)]
  # for each that are less than the median, mark them
  for (ji in jdx) {
    ex3_symbols[ji] <- paste0(ex3_symbols[ji],'_',ji) 
  }
  
}


# lay in the symbols
colnames(ex3_mat) <- ex3_symbols


############ 
# Now to add in the phenotype

idx <- match(rownames(ex3_mat), pheno3$id)
all(rownames(ex3_mat) == pheno3$id[idx])

ex3_mat$SampleID <- pheno3$id[idx]
ex3_mat$Sex <- pheno3$sex[idx]
ex3_mat$Survival <- pheno3$days[idx]
ex3_mat$Censored <- pheno3$status

ex3_mat_F <- ex3_mat[ex3_mat$Sex == 'gender: F',]
ex3_mat_M <- ex3_mat[ex3_mat$Sex == 'gender: M',]

write_csv(ex3_mat_F, file='data/Array Data/jive_validation_13041_ex3_F.csv')
write_csv(ex3_mat_M, file='data/Array Data/jive_validation_13041_ex3_M.csv')


