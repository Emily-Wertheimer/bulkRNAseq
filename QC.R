################ install packages, set dir, import data #################
# install packages & dependencies
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("limma")
BiocManager::install("edgeR")
library(limma)
library(readr)

setwd('/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq')

# import data - opens csv and formats so that ensmbl ID is row names
library(data.table)
gene_counts.df <- fread("/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/2_counts/count_matrices/merged/Count_Matrix_Annotated_FINAL.txt")
#rownames(gene_counts.df) <- gene_counts.df[,1]

# Assuming the first column of gene_counts.df contains the identifiers you want as row names
# and ensuring it's converted from a tibble to a dataframe
gene_counts.df <- as.data.frame(gene_counts.df)

# Set the first column as row names correctly
#rownames(gene_counts.df) <- gene_counts.df[,1]

# Remove the first column properly
gene_counts.df <- gene_counts.df[ , -1, drop = FALSE]

# Ensure row names are unique if necessary
rownames(gene_counts.df) <- make.unique(rownames(gene_counts.df))

# Check the first few rows to verify the structure
head(gene_counts.df)



################ create csv for metadata #################

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
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 5) == "31825"] <- 'ghre + sema'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 5) == "31826"] <- 'lep + sema'

metadata.df$Treatment[substr(metadata.df$sampleName, 1, 6) == 'NGN2A1'] <- 'veh'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 6) == "NGN2A2"] <- 'ghre'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 6) == "NGN2A3"] <- 'lep'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 6) == "NGN2A4"] <- 'sema'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 6) == "NGN2A5"] <- 'ghre + sema'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 6) == "NGN2A6"] <- 'lep + sema'

metadata.df$Treatment[substr(metadata.df$sampleName, 1, 6) == "NGN2B1"] <- 'veh'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 6) == "NGN2B2"] <- 'ghre'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 6) == "NGN2B3"] <- 'lep'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 6) == "NGN2B4"] <- 'sema'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 6) == "NGN2B5"] <- 'ghre + sema'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 6) == "NGN2B6"] <- 'lep + sema'

metadata.df$Treatment[substr(metadata.df$sampleName, 1, 4) == "5531"] <- 'veh'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 4) == "5532"] <- 'ghre'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 4) == "5533"] <- 'lep'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 4) == "5534"] <- 'sema'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 4) == "5535"] <- 'ghre + sema'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 4) == "5536"] <- 'lep + sema'

metadata.df$Treatment[substr(metadata.df$sampleName, 1, 4) == "5As1"] <- 'veh'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 4) == "5As2"] <- 'ghre'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 4) == "5As3"] <- 'lep'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 4) == "5As4"] <- 'sema'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 4) == "5As5"] <- 'ghre + sema'
metadata.df$Treatment[substr(metadata.df$sampleName, 1, 4) == "5As6"] <- 'lep + sema'

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
metadata.df <- subset(metadata.df, select = -RIN.x)

############# annotations ##############
annotations.df <- read.csv("/gpfs/gibbs/pi/huckins/ekw28/resources/gene_annotations.csv")
# pick subset of annotation data that we want
annotations.df <- subset(annotations.df, select=c(ensembl, Chromosome, 
                                                  Gene_name, description))
# make the ensmbl ID into the row name
rownames(annotations.df) <- annotations.df[,1]

############## remove low cpm genes ##############
# remove low counts per million (cpm) genes - this takes away bias
#  if a gene is expressed in 1 transcript but doubles, it's not important but
#  the analysis will think it is. So we should remove it.
# THIS ANALYSIS: we want genes with at least 20 counts in at least 2 of our conditions

library(edgeR)
# Assuming gene_counts.df has gene names as row names and the rest of the data is numeric
gene_counts_matrix <- as.matrix(gene_counts.df)
gene_counts_matrix <- as.numeric(gene_counts_matrix)
gene_cpm.mtx <- cpm(gene_counts_matrix)  # using edgeR function
   #### error in cpm.default(gene_counts)matrix: library sizes should be finite and non-negativeß
