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
gene_counts.df <- read_csv("/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/2_counts/count_matrices/merged/Count_Matrix_Annotated_FINAL.txt")
#rownames(gene_counts.df) <- gene_counts.df[,1]

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
  Sex = rep(NA_character_, nSamples),
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


