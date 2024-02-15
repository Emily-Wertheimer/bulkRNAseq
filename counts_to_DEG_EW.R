# --------------------------------------------------------------------------- IMPORT LIBRARIES

# install packages & dependencies
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
BiocManager::install("edgeR")
library(limma)
library(readr)
library(dplyr)
install.packages("tibble")
library(tibble)
library(edgeR)
library(data.table)
BiocManager::install("variancePartition")


# --------------------------------------------------------------------------- IMPORT RAW COUNTS DATA
setwd('/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq')

# import data - opens csv and formats so that ensmbl ID is row names
gene_counts.dt <- fread("/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/2_counts/count_matrices/merged/Count_Matrix_Annotated_FINAL.txt")
gene_counts.df <- as.data.frame(gene_counts.dt)

annotations.df <- read.csv("/gpfs/gibbs/pi/huckins/ekw28/resources/gene_annotations.csv")
# pick subset of annotation data that we want
annotations.df <- subset(annotations.df, select=c(ensembl, Chromosome, Gene_name, description))

gene_counts.df <- merge(gene_counts.df, annotations.df, by.x='Gene_name', by.y='Gene_name', all.x = TRUE)

# This adds a suffix to make row names unique
gene_counts.df$ensembl <- make.unique(gene_counts.df$ensembl)
# make the ensmbl ID into the row name
rownames(gene_counts.df) <- gene_counts.df$ensembl

# --------------------------------------------------------------------------- CREATE METADATA CSV

# get sample names file
setwd('/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/sample_dir/')
sampleNames_txt <- read_csv(file="/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/sample_dir/unique_filenames.csv", col_names=FALSE)
sampleNames_vec <- as.character(sampleNames_txt$X1)

nSamples <- length(sampleNames_vec)

# Create a metadata data frame with columns initialized to NA of appropriate types
metadata.df <- data.frame(
  sampleName = rep(NA_character_, nSamples),  # Use NA_character_ for character columns
  cellType = rep(NA_character_, nSamples),
  Donor = rep(NA_character_, nSamples),
  Treatment = rep(NA_character_, nSamples),
  Well = rep(NA_character_, nSamples),
  Plate = rep(NA_character_, nSamples),
  Sex = rep(NA_character_, nSamples), RIN = rep(NA_complex_),
  stringsAsFactors = FALSE)  # Ensure that character data does not get converted to factors


# fill first col w/ sample ids
metadata.df$sampleName <- sampleNames_vec

# drop dumb appendages
metadata.df$sampleName <- sub("_.*$", "", metadata.df$sampleName)


## fill other cols using coded sampleName
# cell type
metadata.df$cellType[substr(metadata.df$sampleName, 1, 3) == "318"] <- 'iGABA'
metadata.df$cellType[substr(metadata.df$sampleName, 1, 3) == "NGN"] <- 'iGlut'
metadata.df$cellType[substr(metadata.df$sampleName, 1, 3) == "5As"] <- 'iAstro'
metadata.df$cellType[substr(metadata.df$sampleName, 1, 3) == "553"] <- 'iGABA'

# donor
metadata.df$Donor[substr(metadata.df$sampleName, 1, 3) == "318"] <- '3182'
metadata.df$Donor[substr(metadata.df$sampleName, 1, 5) == "NGN2A"] <- '2607'
metadata.df$Donor[substr(metadata.df$sampleName, 1, 5) == "NGN2B"] <- '553'
metadata.df$Donor[substr(metadata.df$sampleName, 1, 3) == "5As"] <- '553'
metadata.df$Donor[substr(metadata.df$sampleName, 1, 3) == "553"] <- '553'

# tx
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 5) == "31821"] <- 'veh'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 5) == "31822"] <- 'ghre'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 5) == "31823"] <- 'lep'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 5) == "31824"] <- 'sema'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 5) == "31825"] <- 'ghreSema'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 5) == "31826"] <- 'lepSema'

metadata.df$Treatment[substr(metadata.df$sampleName, 1, 6) == 'NGN2A1'] <- 'veh'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 6) == "NGN2A2"] <- 'ghre'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 6) == "NGN2A3"] <- 'lep'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 6) == "NGN2A4"] <- 'sema'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 6) == "NGN2A5"] <- 'ghreSema'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 6) == "NGN2A6"] <- 'lepSema'

metadata.df$Treatment[substr(metadata.df$sampleName, 1, 6) == "NGN2B1"] <- 'veh'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 6) == "NGN2B2"] <- 'ghre'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 6) == "NGN2B3"] <- 'lep'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 6) == "NGN2B4"] <- 'sema'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 6) == "NGN2B5"] <- 'ghreSema'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 6) == "NGN2B6"] <- 'lepSema'

metadata.df$Treatment[substr(metadata.df$sampleName, 1, 4) == "5531"] <- 'veh'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 4) == "5532"] <- 'ghre'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 4) == "5533"] <- 'lep'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 4) == "5534"] <- 'sema'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 4) == "5535"] <- 'ghreSema'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 4) == "5536"] <- 'lepSema'

metadata.df$Treatment[substr(metadata.df$sampleName, 1, 4) == "5As1"] <- 'veh'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 4) == "5As2"] <- 'ghre'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 4) == "5As3"] <- 'lep'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 4) == "5As4"] <- 'sema'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 4) == "5As5"] <- 'ghreSema'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 4) == "5As6"] <- 'lepSema'

# well 
metadata.df$Well[substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)) == "A1"] <- "A1"
metadata.df$Well[substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)) == "A2"] <- "A2"
metadata.df$Well[substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)) == "A3"] <- "A3"
metadata.df$Well[substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)) == "A4"] <- "A4"

metadata.df$Well[substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)) == "B1"] <- "B1"
metadata.df$Well[substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)) == "B2"] <- "B2"
metadata.df$Well[substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)) == "B3"] <- "B3"
metadata.df$Well[substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)) == "B4"] <- "B4"

metadata.df$Well[substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)) == "C1"] <- "C1"
metadata.df$Well[substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)) == "C2"] <- "C2"
metadata.df$Well[substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)) == "C3"] <- "C3"
metadata.df$Well[substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)) == "C4"] <- "C4"

metadata.df$Well[metadata.df$cellType == 'iAstro' & substr(metadata.df$sampleName, 5, 5) == "A" & substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)-1) == "1"] <- "A1"
metadata.df$Well[metadata.df$cellType == 'iAstro' & substr(metadata.df$sampleName, 5, 5) == "A" & substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)-1) == "2"] <- "A2"
metadata.df$Well[metadata.df$cellType == 'iAstro' & substr(metadata.df$sampleName, 5, 5) == "A" & substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)-1) == "3"] <- "A3"
metadata.df$Well[metadata.df$cellType == 'iAstro' & substr(metadata.df$sampleName, 5, 5) == "A" & substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)-1) == "4"] <- "A4"

metadata.df$Well[metadata.df$cellType == 'iAstro' & substr(metadata.df$sampleName, 5, 5) == "B" & substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)-1) == "1"] <- "B1"
metadata.df$Well[metadata.df$cellType == 'iAstro' & substr(metadata.df$sampleName, 5, 5) == "B" & substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)-1) == "2"] <- "B2"
metadata.df$Well[metadata.df$cellType == 'iAstro' & substr(metadata.df$sampleName, 5, 5) == "B" & substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)-1) == "3"] <- "B3"
metadata.df$Well[metadata.df$cellType == 'iAstro' & substr(metadata.df$sampleName, 5, 5) == "B" & substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)-1) == "4"] <- "B4"

metadata.df$Well[metadata.df$cellType == 'iAstro' & substr(metadata.df$sampleName, 5, 5) == "C" & substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)-1) == "1"] <- "C1"
metadata.df$Well[metadata.df$cellType == 'iAstro' & substr(metadata.df$sampleName, 5, 5) == "C" & substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)-1) == "2"] <- "C2"
metadata.df$Well[metadata.df$cellType == 'iAstro' & substr(metadata.df$sampleName, 5, 5) == "C" & substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)-1) == "3"] <- "C3"
metadata.df$Well[metadata.df$cellType == 'iAstro' & substr(metadata.df$sampleName, 5, 5) == "C" & substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)-1) == "4"] <- "C4"

metadata.df$Well[metadata.df$cellType == 'iGABA' & metadata.df$Donor == '553' & substr(metadata.df$sampleName, 5, 5) == "A" & substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)-1) == "1"] <- "A1"
metadata.df$Well[metadata.df$cellType == 'iGABA' & metadata.df$Donor == '553' & substr(metadata.df$sampleName, 5, 5) == "A" & substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)-1) == "2"] <- "A2"
metadata.df$Well[metadata.df$cellType == 'iGABA' & metadata.df$Donor == '553' & substr(metadata.df$sampleName, 5, 5) == "A" & substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)-1) == "3"] <- "A3"
metadata.df$Well[metadata.df$cellType == 'iGABA' & metadata.df$Donor == '553' & substr(metadata.df$sampleName, 5, 5) == "A" & substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)-1) == "4"] <- "A4"

metadata.df$Well[metadata.df$cellType == 'iGABA' & metadata.df$Donor == '553' & substr(metadata.df$sampleName, 5, 5) == "B" & substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)-1) == "1"] <- "B1"
metadata.df$Well[metadata.df$cellType == 'iGABA' & metadata.df$Donor == '553' & substr(metadata.df$sampleName, 5, 5) == "B" & substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)-1) == "2"] <- "B2"
metadata.df$Well[metadata.df$cellType == 'iGABA' & metadata.df$Donor == '553' & substr(metadata.df$sampleName, 5, 5) == "B" & substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)-1) == "3"] <- "B3"
metadata.df$Well[metadata.df$cellType == 'iGABA' & metadata.df$Donor == '553' & substr(metadata.df$sampleName, 5, 5) == "B" & substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)-1) == "4"] <- "B4"

metadata.df$Well[metadata.df$cellType == 'iGABA' & metadata.df$Donor == '553' & substr(metadata.df$sampleName, 5, 5) == "C" & substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)-1) == "1"] <- "C1"
metadata.df$Well[metadata.df$cellType == 'iGABA' & metadata.df$Donor == '553' & substr(metadata.df$sampleName, 5, 5) == "C" & substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)-1) == "2"] <- "C2"
metadata.df$Well[metadata.df$cellType == 'iGABA' & metadata.df$Donor == '553' & substr(metadata.df$sampleName, 5, 5) == "C" & substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)-1) == "3"] <- "C3"
metadata.df$Well[metadata.df$cellType == 'iGABA' & metadata.df$Donor == '553' & substr(metadata.df$sampleName, 5, 5) == "C" & substr(metadata.df$sampleName, nchar(metadata.df$sampleName)-1, nchar(metadata.df$sampleName)-1) == "4"] <- "C4"

metadata.df$Plate[substr(metadata.df$sampleName, nchar(metadata.df$sampleName), nchar(metadata.df$sampleName)) == "A"] <- "A"
metadata.df$Plate[substr(metadata.df$sampleName, nchar(metadata.df$sampleName), nchar(metadata.df$sampleName)) == "B"] <- "B"
metadata.df$Plate[substr(metadata.df$sampleName, nchar(metadata.df$sampleName), nchar(metadata.df$sampleName)) == "C"] <- "C"

# sex
metadata.df$Sex[metadata.df$Donor == '2607'| metadata.df$Donor == '553'] <- "XY"
metadata.df$Sex[metadata.df$Donor == '3182'] <- "XX"

## RIN
# get RIN from YGCA
library(readxl)
glut_gaba.ygca.qc <- read_excel("/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/0_YGCA_QC/RQ21708_Wertheimer_samples 1-59_07262023.xlsx")
astro.ygca.qc <- read_excel("/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/0_YGCA_QC/RQ22282_ekw28_1-18_09182023.xlsx")

# drop leading "_" in sample ID
glut_gaba.ygca.qc$`Sample ID` <- sub(".*?_", "", glut_gaba.ygca.qc$`Sample ID`)
astro.ygca.qc$`Sample Name` <- sub(".*?_", "", astro.ygca.qc$`Sample Name`)

# make big table with RIN and sample ID 
comb.rin.sampID <- data.frame(sampleName = rep(NA_character_, nSamples), RIN = rep(NA_complex_),
                              stringsAsFactors = FALSE)  # Ensure that character data does not get converted to factors

# extract RIN & combine w/ metadata.df
new.vec.ID <- c(glut_gaba.ygca.qc$`Sample ID`, astro.ygca.qc$`Sample Name`)
new.vec.RIN <- c(glut_gaba.ygca.qc$RQN, astro.ygca.qc$RIN)

comb.rin.sampID$sampleName <- new.vec.ID 
comb.rin.sampID$RIN <- new.vec.RIN

metadata.df <- merge(metadata.df, comb.rin.sampID, by='sampleName')
metadata.df <- subset(metadata.df, select = -RIN.y)

# --------------------------------------------------------------------------- SEPARATE BY CELL TYPE
## gene counts
# remove leading X if applicable 
colnames(gene_counts.df) <- gsub("^X", "", sub("_.*", "", colnames(gene_counts.df)))

# subset by sample name code 
gene_counts.df.iGABA <- gene_counts.df[, substr(colnames(gene_counts.df), 1, 3) %in% c('318', '553')]
gene_counts.df.iGlut <- gene_counts.df[, substr(colnames(gene_counts.df), 1, 3) == 'NGN']
gene_counts.df.iAstro <- gene_counts.df[, substr(colnames(gene_counts.df), 1, 3) == '5As']


## meta
metadata.df.iGABA <- metadata.df[metadata.df$cellType == 'iGABA',]
metadata.df.iGlut <- metadata.df[metadata.df$cellType == 'iGlut',]
metadata.df.iAstro <- metadata.df[metadata.df$cellType == 'iAstro',]

# --------------------------------------------------------------------------- REMOVE LOW CPM GENES

# this takes away bias
#  if a gene is expressed in 1 transcript but doubles, it's not important but
#  the analysis will think it is. So we should remove it.
# THIS ANALYSIS: we want genes with at least 20 counts in at least 2 of our conditions

# make sure all data in gene_counts.df is numeric (ie no gene name, description, or chromosome cols); also remove GeneID
#gene_counts.df.iGABA.numeric <- subset(gene_counts.df.iGABA, select = c(-Gene, -GeneID, -description, -Chromosome, -ensembl)) #-Gene_name,
#gene_counts.df.iGlut.numeric <- subset(gene_counts.df.iGlut, select = c(-Gene, -GeneID, -description, -Chromosome, -ensembl)) # -Gene_name,
#gene_counts.df.iAstro.numeric <- subset(gene_counts.df.iAstro, select = c(-Gene, -GeneID, -description, -Chromosome, -ensembl)) #-Gene_name,

# convert to matrix
gene_counts_matrix_iGABA <- as.matrix(gene_counts.df.iGABA)
gene_counts_matrix_iGlut <- as.matrix(gene_counts.df.iGlut)
gene_counts_matrix_iAstro <- as.matrix(gene_counts.df.iAstro)

# run cpm 
gene_cpm.mtx.iGABA <- cpm(gene_counts_matrix_iGABA)
gene_cpm.mtx.iGlut <- cpm(gene_counts_matrix_iGlut)
gene_cpm.mtx.iAstro <- cpm(gene_counts_matrix_iAstro)

# high pass filter
genes_over_cpm_threshold.iGABA.mtx <- gene_cpm.mtx.iGABA > 0.5  # 0.5 is our threshold; this returns a matrix of booleans 
genes_over_cpm_threshold.iGlut.mtx <- gene_cpm.mtx.iGlut > 0.5 
genes_over_cpm_threshold.iAstro.mtx <- gene_cpm.mtx.iAstro > 0.5 

# how many of the samples have cpm >0.5 for each gene?
samples_over_thresh.iGABA.int <- rowSums(genes_over_cpm_threshold.iGABA.mtx)
samples_over_thresh.iGlut.int <- rowSums(genes_over_cpm_threshold.iGlut.mtx)
samples_over_thresh.iAstro.int <- rowSums(genes_over_cpm_threshold.iAstro.mtx)


# subset the matrix only for the genes that are >0.5 cpm in more than 2 samples
keep.iGABA <- samples_over_thresh.iGABA.int > 2
keep.iGlut <- samples_over_thresh.iGlut.int > 2
keep.iAstro <- samples_over_thresh.iAstro.int > 2


gene_counts_aftercpm.df.iGABA <- gene_counts.df.iGABA[keep.iGABA,]
gene_counts_aftercpm.df.iGlut <- gene_counts.df.iGlut[keep.iGlut,]
gene_counts_aftercpm.df.iAstro <- gene_counts.df.iAstro[keep.iAstro,]



## check for validation that the 0.5 cpm corresponds to 20 counts
#plot1 <- plot(gene_cpm.mtx[,9], gene_counts.df[,9], ylim=c(0,50), xlim=c(0,3))
#plot1<-abline(v=0.5, h=20)

# --------------------------------------------------------------------------- DGE OBJECT CREATION

# The DGElist holds all of the data you want to analyze including:
#    counts, library size to normalize to, normalization factors, experimental
#    conditions, gene annotations

data.dge.iGABA <-  DGEList(gene_counts_aftercpm.df.iGABA) 
data.dge.iGlut <-  DGEList(gene_counts_aftercpm.df.iGlut) 
data.dge.iAstro <-  DGEList(gene_counts_aftercpm.df.iAstro) 


# add experimental condition metadata to the object
data.dge.iGABA$samples$group <- metadata.df.iGABA$Treatment 
data.dge.iGlut$samples$group <- metadata.df.iGlut$Treatment 
data.dge.iAstro$samples$group <- metadata.df.iAstro$Treatment 

## add annotation data to the dge object
# match the annotations to the dge list object rownames
annotations.df.iGABA <- annotations.df[match(rownames(data.dge.iGABA), annotations.df$ensembl),]
annotations.df.iGlut <- annotations.df[match(rownames(data.dge.iGlut), annotations.df$ensembl),]
annotations.df.iAstro <- annotations.df[match(rownames(data.dge.iAstro), annotations.df$ensembl),]


data.dge.iGABA$genes <- annotations.df.iGABA
data.dge.iGlut$genes <- annotations.df.iGlut
data.dge.iAstro$genes <- annotations.df.iAstro


# --------------------------------------------------------------------------- TMM NORMALIZATION
# TMM normalizes library sizes. If you do RNA seq on 2 samples and happen to get
# better reads on one, it will have a bigger library size.
# TMM accounts for this by defining the right library size, and then normalizing
# all other library sizes to that one. So it may multiply libraries that are
# too small by, say, 1.3, and those that are too large by, say, 0.8

data.dge.iGABA <- calcNormFactors(data.dge.iGABA) 
data.dge.iGlut <- calcNormFactors(data.dge.iGlut) 
data.dge.iAstro <- calcNormFactors(data.dge.iAstro) 


# shows how each library size was normalized
data.dge.iGABA$samples$norm.factors 
data.dge.iGlut$samples$norm.factors
data.dge.iAstro$samples$norm.factors

# --------------------------------------------------------------------------- DESIGN MATRIX
metadata.df.iGABA$Treatment <- factor(metadata.df.iGABA$Treatment)
#levels(metadata.df$Treatment) 
metadata.df.iGlut$Treatment <- factor(metadata.df.iGlut$Treatment)
metadata.df.iAstro$Treatment <- factor(metadata.df.iAstro$Treatment)

# refactor reference level (to compare against veh)
metadata.df.iGABA$Treatment <- relevel(metadata.df.iGABA$Treatment, ref = 'veh')
metadata.df.iGlut$Treatment <- relevel(metadata.df.iGlut$Treatment, ref = 'veh')
metadata.df.iAstro$Treatment <- relevel(metadata.df.iAstro$Treatment, ref = 'veh')


# define what to regress out for
design_matrix_iGABA <- model.matrix( ~ 0 + Treatment + Donor, data = metadata.df.iGABA) 
design_matrix_iGlut <- model.matrix(~ 0 + Treatment + Donor, data = metadata.df.iGlut) 
design_matrix_iAstro <- model.matrix(~ 0 + Treatment, data = metadata.df.iAstro) # only have 1 donor of astros


# --------------------------------------------------------------------------- VOOM TRANSFORM
# voom is a limma function that transforms the data into log2cpm which is what
# the limma package uses for RNA seq analysis. This allows your data to be
# in the right format for limma's math. 

data.iGABA.voom <- voom(counts=data.dge.iGABA, design=design_matrix_iGABA, plot=TRUE) # limma package function
data.iGlut.voom <- voom(counts=data.dge.iGlut, design=design_matrix_iGlut, plot=TRUE) # limma package function
data.iAstro.voom <- voom(counts=data.dge.iAstro, design=design_matrix_iAstro, plot=TRUE) # limma package function

# --------------------------------------------------------------------------- VARIANCE PARTITION
################# run variance partition to select covariates ################
# Check covariate corr w/ PCs 
# Variance partition 
# in how many genes is expression associated w/ a given covariate (include factors (metadata.df cols) that >1% var in 10% genes)
# make sure to include in model these covariates (>1% in 10% genes); if less than that threshold then can justify removing from subsequent model

library('variancePartition')

# select only cols of interest from metadata.df
meta <- subset(metadata.df, select=c('sampleName', 'Donor', 'Treatment', 'Well', 'cellType' ))

# load voom transformed expression (DGE) matrix
voom_matrix <- data.voom
# change voom mat sample names to match meta sample names
colnames(voom_matrix) <- gsub("^X", "", sub("_.*", "", colnames(voom_matrix)))

# model (all vars are categorical --> modeled as random effects)
form <- ~ (1|Donor) + (1|Treatment) + (1|Well) + (1|cellType)

# 1.  fit model and extract results
# since all vars are categorical, linear model is used 
# e/ entry  in results is a regression model fit on single gene

# 2. exact variance fractions from e/ model fit 
# for e/ gene, returns fraction of variance attributable to each variable

# INTERPRETATION
#The variance explained by e/ variable after correcting for all other variables
varPart <- fitExtractVarPartModel(voom_matrix, form, meta)

# sort vars by median fraction of variance explained
vp <- sortCols(varPart)

# write file
setwd('/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/')
write.table(vp, sprintf("variancePartition_output.txt"), sep="\t", quote=F, row.names=F)

################# visualize variance partition ################
library(ggplot2)
library(reshape2)

col =  c(ggColorHue(ncol(vp) - 1), "grey85")
vp2=(vp*100)
names=colnames(vp2)
t=apply(vp2, 2, FUN=median)
t=round(t,2)
t[t == 0.00] <- "<0.1"
textMatrix = paste(names, "\n(", t, "%)\n", sep = "");
melted=melt(vp2)

p=ggplot(melted, aes(x=variable, y=value,  color=variable)) +
  geom_violin(aes(x = variable, y = value), scale = "width") +
  geom_boxplot(aes(x = variable, y = value), width = .1,outlier.size = 0.4)+
  theme_bw()+ylim(0,100)+
  theme(#panel.border = element_blank(),
    axis.line.x = element_line(linewidth = 0.3, linetype = "solid", colour = "black"),
    axis.line.y = element_line(linewidth = 0.3, linetype = "solid", colour = "black"),
    axis.text.x = element_text(size=8, angle=45,colour="black",hjust=1),
    axis.text.y = element_text(size=7,colour="black"))+
  scale_x_discrete(labels=textMatrix) + ylab("Variance Explained (%)") 

png(sprintf("variancePartition.png"), res=300, units="in", width=10, height=7)
p

################# choose covariates ################
#alanna's criteria: continuous technical factors that explained ≥ 1% variation in ≥ 10% of genes

cat("Choosing Covariates\n")
test=colnames(vp)[1:length(colnames(vp))-1]
covs_to_adjust=c()
for(col_i in test) {
  print(col_i)
  if(length(which(vp[,col_i] >= 0.01))/(nrow(vp)) > 0.1) covs_to_adjust=c(covs_to_adjust,col_i)
}

write.table(covs_to_adjust, file=sprintf("lep_ghre_sema_allVarPartCovs.txt"), quote=F, sep="\n", col.names=T, row.names=F)

# --------------------------------------------------------------------------- LM FIT
# We have a general model that explains RNA expression and we want to fit our
# data to that model. This allows us to use limma to make contrasts

data.fit.iGABA <- lmFit(object=data.iGABA.voom, design=design_matrix_iGABA) # limma package function
data.fit.iGlut <- lmFit(object=data.iGlut.voom, design=design_matrix_iGlut) # limma package function
data.fit.iAstro <- lmFit(object=data.iAstro.voom, design=design_matrix_iAstro) # limma package function

# --------------------------------------------------------------------------- CONTRAST MATRIX
## make contrast to assess effects of different treatments (just start w/ lep, ghre, sema each v. veh, then do combinatorial/synergy)

cm_iGABA <-makeContrasts( 
  lep_v_veh = (Treatmentlep - Treatmentveh),
  ghre_v_veh = (Treatmentghre - Treatmentveh), # ghre isn't included - was it filtered out? (factor dropping?)
  sema_v_veh = (Treatmentsema - Treatmentveh),
  levels=design_matrix_iGABA)

cm_iGlut <-makeContrasts( 
  lep_v_veh = (Treatmentlep - Treatmentveh),
  ghre_v_veh = (Treatmentghre - Treatmentveh), # ghre isn't included - was it filtered out? (factor dropping?)
  sema_v_veh = (Treatmentsema - Treatmentveh),
  levels=design_matrix_iGlut)

cm_iAstro <-makeContrasts( 
  lep_v_veh = (Treatmentlep - Treatmentveh),
  ghre_v_veh = (Treatmentghre - Treatmentveh), # ghre isn't included - was it filtered out? (factor dropping?)
  sema_v_veh = (Treatmentsema - Treatmentveh),
  levels=design_matrix_iAstro)

# add contrasts to the Lmfit object. 
fit.iGABA <- contrasts.fit(fit=data.fit.iGABA, contrasts=cm_iGABA)
fit.iGlut <- contrasts.fit(fit=data.fit.iGlut, contrasts=cm_iGlut)
fit.iAstro <- contrasts.fit(fit=data.fit.iAstro, contrasts=cm_iAstro)


# perform bayes shrinkage to the analysis and estimate modified t and p values
fitDupCor.iGABA <- eBayes(fit.iGABA)
fitDupCor.iGlut <- eBayes(fit.iGlut)
fitDupCor.iAstro <- eBayes(fit.iAstro)


# --------------------------------------------------------------------------  FIND DEGs
# find most differentially expressed genes
DGE_lep_v_veh.iGABA <- topTable(fitDupCor.iGABA, coef="lep_v_veh", n=80000)
DGE_sema_v_veh.iGABA <- topTable(fitDupCor.iGABA, coef="sema_v_veh", n=80000) 
DGE_ghre_v_veh.iGABA <- topTable(fitDupCor.iGABA, coef="ghre_v_veh", n=80000) 

DGE_lep_v_veh.iGlut <- topTable(fitDupCor.iGlut, coef="lep_v_veh", n=80000)
DGE_sema_v_veh.iGlut <- topTable(fitDupCor.iGlut, coef="sema_v_veh", n=80000) 
DGE_ghre_v_veh.iGlut <- topTable(fitDupCor.iGlut, coef="ghre_v_veh", n=80000)

DGE_lep_v_veh.iAstro <- topTable(fitDupCor.iAstro, coef="lep_v_veh", n=80000)
DGE_sema_v_veh.iAstro <- topTable(fitDupCor.iAstro, coef="sema_v_veh", n=80000) 
DGE_ghre_v_veh.iAstro <- topTable(fitDupCor.iAstro, coef="ghre_v_veh", n=80000)

# add annotations to the table
DGE_lep_v_veh.iGABA <- merge(DGE_lep_v_veh.iGABA, annotations.df.iGABA, 
                             by.x="ensembl", by.y="ensembl")
DGE_sema_v_veh.iGABA <- merge(DGE_sema_v_veh.iGABA, annotations.df.iGABA, by.x="ensembl",
                              by.y="ensembl")
DGE_ghre_v_veh.iGABA <- merge(DGE_ghre_v_veh.iGABA, annotations.df.iGABA, by.x="ensembl",
                              by.y="ensembl")

DGE_lep_v_veh.iGlut <- merge(DGE_lep_v_veh.iGlut, annotations.df.iGlut, by.x="ensembl",
                             by.y="ensembl")
DGE_sema_v_veh.iGlut <- merge(DGE_sema_v_veh.iGlut, annotations.df.iGlut, by.x="ensembl",
                              by.y="ensembl")
DGE_ghre_v_veh.iGlut <- merge(DGE_ghre_v_veh.iGlut, annotations.df.iGlut, by.x="ensembl",
                              by.y="ensembl")

DGE_lep_v_veh.iAstro <- merge(DGE_lep_v_veh.iAstro, annotations.df.iAstro, by.x="ensembl",
                             by.y="ensembl")
DGE_sema_v_veh.iAstro <- merge(DGE_sema_v_veh.iAstro, annotations.df.iAstro, by.x="ensembl",
                              by.y="ensembl")
DGE_ghre_v_veh.iAstro <- merge(DGE_ghre_v_veh.iAstro, annotations.df.iAstro, by.x="ensembl",
                              by.y="ensembl")

# save DEGs
setwd('/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/3_DEG/')
write.table(DGE_lep_v_veh.iGABA, "DEG_response_iGABA_lep_v_veh.txt", sep="\t")
write.table(DGE_sema_v_veh.iGABA, "DEG_response_iGABA_sema_v_veh.txt", sep="\t")
write.table(DGE_ghre_v_veh.iGABA, "DEG_response_iGABA_ghre_v_veh.txt", sep="\t")

write.table(DGE_lep_v_veh.iGlut, "DEG_response_iGlut_lep_v_veh.txt", sep="\t")
write.table(DGE_sema_v_veh.iGlut, "DEG_response_iGlut_sema_v_veh.txt", sep="\t")
write.table(DGE_ghre_v_veh.iGlut, "DEG_response_iGlut_ghre_v_veh.txt", sep="\t")

write.table(DGE_lep_v_veh.iAstro, "DEG_response_iAstro_lep_v_veh.txt", sep="\t")
write.table(DGE_sema_v_veh.iAstro, "DEG_response_iAstro_sema_v_veh.txt", sep="\t")
write.table(DGE_ghre_v_veh.iAstro, "DEG_response_iAstro_ghre_v_veh.txt", sep="\t")   

# --------------------------------------------------------------------------- EXTRACT SIG DEGS
sigDEG_lep_v_veh.iGABA <- DGE_lep_v_veh.iGABA[DGE_lep_v_veh.iGABA$adj.P.Val <= 0.05, ]
sigDEG_sema_v_veh.iGABA <- DGE_sema_v_veh.iGABA[DGE_sema_v_veh.iGABA$adj.P.Val <= 0.05, ]
sigDEG_ghre_v_veh.iGABA <- DGE_ghre_v_veh.iGABA[DGE_ghre_v_veh.iGABA$adj.P.Val <= 0.05, ] 

sigDEG_lep_v_veh.iGlut <- DGE_lep_v_veh.iGlut[DGE_lep_v_veh.iGlut$adj.P.Val <= 0.05, ]
sigDEG_sema_v_veh.iGlut <- DGE_sema_v_veh.iGlut[DGE_sema_v_veh.iGlut$adj.P.Val <= 0.05, ]
sigDEG_ghre_v_veh.iGlut <- DGE_ghre_v_veh.iGlut[DGE_ghre_v_veh.iGlut$adj.P.Val <= 0.05, ] 

sigDEG_lep_v_veh.iAstro <- DGE_lep_v_veh.iAstro[DGE_lep_v_veh.iAstro$adj.P.Val <= 0.05, ]
sigDEG_sema_v_veh.iAstro <- DGE_sema_v_veh.iAstro[DGE_sema_v_veh.iAstro$adj.P.Val <= 0.05, ]
sigDEG_ghre_v_veh.iAstro <- DGE_ghre_v_veh.iAstro[DGE_ghre_v_veh.iAstro$adj.P.Val <= 0.05, ] 
