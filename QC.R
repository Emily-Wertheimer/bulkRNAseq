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

# --------------------------------------------------------------------------- IMPORT RAW COUTNS DATA
setwd('/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq')

# import data - opens csv and formats so that ensmbl ID is row names
library(data.table)
gene_counts.dt <- fread("/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/2_counts/count_matrices/merged/Count_Matrix_Annotated_FINAL.txt")
gene_counts.df <- as.data.frame(gene_counts.dt)

annotations.df <- read.csv("/gpfs/gibbs/pi/huckins/ekw28/resources/gene_annotations.csv")
# pick subset of annotation data that we want
annotations.df <- subset(annotations.df, select=c(ensembl, Chromosome, 
                                                  Gene_name, description))

gene_counts.df <- merge(gene_counts.df, annotations.df, by.x='Gene_name', by.y='Gene_name', all.x = TRUE)

# This adds a suffix to make row names unique
gene_counts.df$ensembl <- make.unique(gene_counts.df$ensembl)
# make the ensmbl ID into the row name
gene_counts.df <- gene_counts.df %>% remove_rownames() %>% column_to_rownames(var="ensembl")

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

# extract RIN 
new.vec.ID <- c(glut_gaba.ygca.qc$`Sample ID`, astro.ygca.qc$`Sample Name`)
new.vec.RIN <- c(glut_gaba.ygca.qc$RQN, astro.ygca.qc$RIN)

comb.rin.sampID$sampleName <- new.vec.ID 
comb.rin.sampID$RIN <- new.vec.RIN

metadata.df <- merge(metadata.df, comb.rin.sampID, by='sampleName')
metadata.df <- subset(metadata.df, select = -RIN.y)

# --------------------------------------------------------------------------- REMOVE LOW CPM GENES

# this takes away bias
#  if a gene is expressed in 1 transcript but doubles, it's not important but
#  the analysis will think it is. So we should remove it.
# THIS ANALYSIS: we want genes with at least 20 counts in at least 2 of our conditions

# make sure all data in gene_counts.df is numeric (ie no gene name, description, or chromosome cols); also remove GeneID
gene_counts.df <- subset(gene_counts.df, select = c(-Chromosome, -Gene_name, -description, -GeneID))

# convert to matrix
gene_counts_matrix <- as.matrix(gene_counts.df)

# run cpm 
gene_cpm.mtx <- cpm(gene_counts_matrix)

# high pass filter
genes_over_cpm_threshold.mtx <- gene_cpm.mtx > 0.5  # 0.5 is our threshold; this returns a matrix of booleans 

# how many of the samples have cpm >0.5 for each gene?
samples_over_thresh.int <- rowSums(genes_over_cpm_threshold.mtx)

# subset the matrix only for the genes that are >0.5 cpm in more than 2 samples
keep <- samples_over_thresh.int > 2
gene_counts_aftercpm.df <- gene_counts.df[keep,]

# check for validation that the 0.5 cpm corresponds to 20 counts
plot1 <- plot(gene_cpm.mtx[,8], gene_counts.df[,8], ylim=c(0,50), xlim=c(0,3))
plot1<-abline(v=0.5, h=20)
# --------------------------------------------------------------------------- DGE OBJECT CREATION

# The DGElist holds all of the data you want to analyze including:
#    counts, library size to normalize to, normalization factors, experimental
#    conditions, gene annotations

data.dge <-  DGEList(gene_counts_aftercpm.df) 

# add experimental condition metadata to the object
#data.dge$samples$group <- meta_data.df$mod.gene
metadata.df <- metadata.df %>% distinct()
data.dge$samples$group <- metadata.df$Treatment # is this my equivalent of "mod.gene"?

# add annotation data to the object
# first, match the annotations to the dge list object rownames
annotations.df <- annotations.df[match(rownames(data.dge), rownames(annotations.df)),]
data.dge$genes <- annotations.df

# --------------------------------------------------------------------------- TMM NORMALIZATION
  # TMM normalizes library sizes. If you do RNA seq on 2 samples and happen to get
  # better reads on one, it will have a bigger library size.
  # TMM accounts for this by defining the right library size, and then normalizing
  # all other library sizes to that one. So it may multiply libraries that are
  # too small by, say, 1.3, and those that are too large by, say, 0.8

data.dge <- calcNormFactors(data.dge) 

# shows how each library size was normalized
data.dge$samples$norm.factors 

# --------------------------------------------------------------------------- DESIGN MATRIX
# define what to regress out for
design_matrix <- model.matrix(~ 0 + cellType + Donor + Treatment + Well, data = metadata.df) # did not include sex (complete colinearity with line) or RIN (missing data)

# --------------------------------------------------------------------------- VOOM TRANSFORM
# voom is a limma function that transforms the data into log2cpm which is what
# the limma package uses for RNA seq analysis. This allows your data to be
# in the right format for limma's math. 

data.voom <- voom(counts=data.dge, design=design_matrix, plot=TRUE) # limma package function


# --------------------------------------------------------------------------- LM FIT
# We have a general model that explains RNA expression and we want to fit our
# data to that model. This allows us to use limma to make contrasts

data.fit <- lmFit(object=data.voom, design=design_matrix) # limma package function

# --------------------------------------------------------------------------- CONTRAST MATRIX
# make contrast to assess effects of different treatments (just start w/ lep, ghre, sema each v. veh)
cm <-makeContrasts( 
  lep_v_veh = (lep - veh),
  ghre_v_veh = (ghre - veh),
  sema_v_veh = (sema - veh),
  levels=design_matrix)

# add contrasts to the Lmfit object. 
fit <- contrasts.fit(fit=data.fit, contrasts=cm)

# perform bayes shrinkage to the analysis and estimate modified t and p values
fitDupCor <- eBayes(fit)


