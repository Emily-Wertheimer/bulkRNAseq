library(limma)
library(edgeR)

# import data - opens csv and formats so that ensmbl ID is row names

# COUNTS:
gene_counts.df <- read.csv(file="counts-multi.csv", header=TRUE, sep=",")
rownames(gene_counts.df) <- gene_counts.df[,1]

# META DATA:
meta_data.df <- read.csv(file="meta-multi.csv", header=TRUE, sep=",")

# ANNOTATIONS:
annotations.df <- read.csv(file="anno.csv", header=TRUE, sep=",")
# pick subset of annotation data that we want
annotations.df <- subset(annotations.df, select=c(ensembl, Chromosome, 
                                                  Gene_name, description))
# make the ensmbl ID into the row name
rownames(annotations.df) <- annotations.df[,1]


# ---------------------------------------------------------------------
# remove low counts per million (cpm) genes - this takes away bias
#  if a gene is expressed in 1 transcript but doubles, it's not important but
#  the analysis will think it is. So we should remove it.
# THIS ANALYSIS: we want genes with at least 20 counts in at least 2 of our 
#                conditions

gene_cpm.mtx <- cpm(gene_counts.df)  # using edgeR function
genes_over_cpm_threshold.mtx <- gene_cpm.mtx > 0.5  # 0.5 is our threshold
# this returns a matrix of booleans 

# how many of the samples have cpm >0.5 for each gene (34 samples total)
samples_over_thresh.int <- rowSums(genes_over_cpm_threshold.mtx)
keep <- samples_over_thresh.int > 2
# subset the matrix only for the genes that are >0.5 cpm in more than 2 samples
gene_counts_aftercpm.df <- gene_counts.df[keep,]

# check for validation that the 0.5 cpm corresponds to 20 counts
plot(gene_cpm.mtx[,1], gene_counts.df[,1], ylim=c(0,50), xlim=c(0,3))
abline(v=0.5, h=20)

# --------------------------------------------------------------------------
# Create DGElist object. 
# The DGElist holds all of the data you want to analyze including:
#    counts, library size to normalize to, normalization factors, experimental
#    conditions, gene annotations
data.dge <-  DGEList(gene_counts_aftercpm.df) # edgeR library function
# add experimental condition groups to the object
data.dge$samples$group <- meta_data.df$mod.gene
# add annotation data to the object
# first, match the annotations to the dge list object rownames
annotations.df <- annotations.df[match(rownames(data.dge), 
                                       rownames(annotations.df)),]
data.dge$genes <- annotations.df


# --------------------------------------------------------------------------
# TMM Normalization
# TMM normalizes library sizes. If you do RNA seq on 2 samples and happen to get
# better reads on one, it will have a bigger library size.
# TMM accounts for this by defining the right library size, and then normalizing
# all other library sizes to that one. So it may multiply libraries that are
# too small by, say, 1.3, and those that are too large by, say, 0.8
data.dge <- calcNormFactors(data.dge) # edgeR package function
data.dge$samples$norm.factors # this shows you how it normalized each library size

# ---------------------------------------------------------------------------
# Create design matrix
# define what to account for
design_matrix <- model.matrix(~ 0 + meta_data.df$dose + meta_data.df$donor
                              + meta_data.df$line) # basic R function

# ---------------------------------------------------------------------------
# Voom transform
# voom is a limma function that transforms the data into log2cpm which is what
# the limma package uses for RNA seq analysis. This allows your data to be
# in the right format for limma's math. 

data.voom <- voom(counts=data.dge, design=design_matrix, plot=TRUE) # limma package function

# ---------------------------------------------------------------------------
# Lmfit
# We have a general model that explains RNA expression and we want to fit our
# data to that model. This allows us to use limma to make contrasts

data.fit <- lmFit(object=data.voom, design=design_matrix) # limma package function

# ---------------------------------------------------------------------------
# Create contrast matrix

cm <-makeContrasts( #make a contrast to test the effect of different HCort doses
    Dose_100 = (Dose100 - Dose0),
    Dose_1000 = (Dose1000 - Dose0),
levels=design)

# add contrasts to the Lmfit object. 
fit <- contrasts.fit(fit=data.fit, contrasts=cm)

# perform bayes shrinkage to the analysis and estimate modified t and p values
fitDupCor <- eBayes(fit)

# --------------------------------------------------------------------------  
# find most differentially expressed genes
DGE_NGN2_100nm <- topTable(fitDupCor, coef="Dose_100", n=80000)
DGE_NGN2_1000nm <- topTable(fitDupCor, coef="Dose_1000", n=80000)

# add annotations to the table
DGE_NGN2_100nm <- merge(DGE_NGN2_100nm, annotations.df, by.x="row.names",
                               y.y="ensembl")
DGE_NGN2_1000nm <- merge(DGE_NGN2_1000nm, annotations.df, by.x="row.names",
                               y.y="ensembl")

# save DEGs
write.table(DGE_NGN2_100nm, "DEG_response_100nM.txt", sep="\t")
write.table(DGE_NGN2_1000nm, "DEG_response_1000.txt", sep="\t")

