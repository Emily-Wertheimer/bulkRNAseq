#--------------------
# title:"edgeR & limma"
# author: "Joy Lee"
# date: "10/3/2022"
​
#--------------------
### library
library(limma)
library(edgeR)
library(tidyverse)
# for gradient color
library(grDevices)
# for heatmap
library(pheatmap)
​
#--------------------
### data import
## set working directory
setwd("/Users/cl2375/OneDrive - Yale University/PhD study/Labs/Brennand Lab/RNA-seq/RNAseq_iMGLpilot/filter_count/")
## import count data
gene_counts.df <- read.table(file = "3182_filtered_genes_counts.txt", head = TRUE, sep = "\t", row.names = 1)
## clean up column & row name
colnames(gene_counts.df) <- gsub("X", "", colnames(gene_counts.df))
rownames(gene_counts.df) <- str_split_fixed(rownames(gene_counts.df), "\\.", 2)[,1]
## data column reordering
sample_orders <- c("iPS", "HPC", "D8", "D16", "D24", "M5", "M10")
sample_date <- c(0,12,20,28,36,41,46)
column_orders<-paste0("3182_",rep(sample_orders,each=3),"_",c(1:3))
gene_counts.df <- gene_counts.df[,column_orders]
## import annotation
annotations.df <- read.csv(file="anno.csv", header=TRUE, sep=",")
# pick subset of annotation data that we want
annotations.df <- subset(annotations.df, select=c(ensembl, Chromosome, Gene_name, description))
# make the ensmbl ID into the row name
rownames(annotations.df) <- annotations.df[,1]
​
### clean up data
## remove low counts per million (cpm) genes, retain genes with at least 20 counts in at least 3 replicates
gene_cpm.mtx <- cpm(gene_counts.df)
genes_over_cpm_threshold.mtx <- gene_cpm.mtx > 0.5 # threhold 0.5 and return booleans matrix
samples_over_thresh.int <- rowSums(genes_over_cpm_threshold.mtx)
keep <- samples_over_thresh.int > 3
## subset the matrix only for the genes that are >0.5 cpm in more than 2 samples
gene_counts_aftercpm.df <- gene_counts.df[keep,]
## check for validation that the 0.5 cpm corresponds to 20 counts
plot(gene_cpm.mtx[,1], gene_counts.df[,1], ylim=c(0,50), xlim=c(0,3))
abline(v=0.5, h=20)
​
#--------------------
### create DEGlist object
## The DGElist holds all of the data you want to analyze including: counts, library size to normalize to, normalization factors, experimental conditions, gene annotations
##  group is the experimental condition
data.dge <-  DGEList(gene_counts_aftercpm.df,group=rep(sample_orders,each=3))
## add in genes annotation
# check the percentage of genes covered by annotation list
# sum(rownames(gene_counts_aftercpm.df) %in% rownames(annotations.df)) / nrow(gene_counts_aftercpm.df)
# subset the annotation needed
annotations.df <- annotations.df[match(rownames(data.dge), rownames(annotations.df)),]
# add annotation information into DEGlist
data.dge$genes <- annotations.df
​
#--------------------
### TMM normalization
## TMM normalizes library sizes. TMM accounts for this by defining the right library size, and then normalizing all other library sizes to that one with multiplying factors
data.dge <- calcNormFactors(data.dge) 
## add the normalized factor of each library size into metadata
data.dge$samples$norm.factors 
​
### --------------------
### want to see the distribution of all sample
## simple plot
plotMDS(data.dge, labels=data.dge$samples$group)
# focus the plot on the iMGL development
plotMDS(data.dge, labels=data.dge$samples$group, xlim=c(1.5,2.2))
​
## colored plot
# create a color palette
colfunc <- colorRampPalette(c("light blue", "dark blue"))
# plot
par(mar=c(5,5,4,1)+.1)
plotMDS(data.dge, pch = 19, col =rep(colfunc(7),each=3), cex = 4.5,
        cex.axis = 2,cex.lab = 2)
plotMDS(data.dge, pch = 19, col =rep(colfunc(7),each=3), cex = 4.5,
        cex.axis = 2,cex.lab = 2, ylim=c(-1,1), xlim=c(1.4,2.2))
# legend
plot(rep(1,21),col =rep(colfunc(7),each=3), pch = 19, cex = 2)
​
#--------------------
### Design matrix
# include all covariates in design
design_matrix <- model.matrix(~ 0 + data.dge$samples$group)
colnames(design_matrix) <- c("D16","D24","D8","HPC","iPS","M10","M5")
​
#--------------------
### Voom transform
## voom is a limma function that transforms the data into log2cpm for RNA seq analysis in limma.
data.voom <- voom(counts=data.dge, design=design_matrix, plot=TRUE)
​
# ---------------------------------------------------------------------------
### Lmfit
# We have a general model that explains RNA expression and we want to fit our data to that model.
data.fit <- lmFit(object=data.voom, design=design_matrix) # limma package function
​
​
# ---------------------------------------------------------------------------
### time-series design (from DEGlist)
data.dge.time <-  DGEList(gene_counts_aftercpm.df,group=rep(sample_date,each=3))
annotations.df <- annotations.df[match(rownames(data.dge.time), rownames(annotations.df)),]
data.dge.time$genes <- annotations.df
data.dge.time <- calcNormFactors(data.dge.time)
data.dge.time$samples$norm.factors
# polynomial with 3 degrees of freedom
X <- poly(rep(sample_date, each=3), degree=3)
design_time <- model.matrix(~X)
# estimating the NB dispersion
data.dge.time <- estimateDisp(data.dge.time,design_time)
sqrt(data.dge.time$common.dispersion)
plotBCV(data.dge.time)
# estimating the QL dispersion
fit <- glmQLFit(data.dge.time,design_time, robust=TRUE)
plotQLDisp(fit)
#--------------------
### Time course trend analysis
##  looking for genes that change expression level over time
##  the design matrix uses 3 natural spline basis vectors to model smooth changes over time
##  test for a trend by conducting F-tests on 3 df for each gene
# run model
fit <-glmQLFTest(fit, coef=2:4)
# extract top set of genes with most significant time effects
tab <- as.data.frame(topTags(fit, n=30))
## visualize fitted spline curves fr the top four genes
# computing observed and fitted log-CPM value for each gene
logCPM.obs <- cpm(data.dge.time, log=TRUE, prior.count=fit$prior.count)
logCPM.fit <- cpm(fit, log=TRUE)
# loop through the first four genes and plot
par(mfrow=c(2,2))
for(i in 1:4) {
  FlybaseID <- row.names(tab)[i]
  Symbol <- tab$Gene_name[i]
  logCPM.obs.i <- logCPM.obs[FlybaseID,]
  logCPM.fit.i <- logCPM.fit[FlybaseID,]
  plot(rep(sample_date, each=3), logCPM.obs.i, ylab="log-CPM", main=Symbol, pch=16, xlab="Days")
  lines(rep(sample_date, each=3), logCPM.fit.i, col="red", lwd=2)}
# clear the figure setting: par(mfrow=c(1,1))
# ---------------------------------------------------------------------------
### time-series design (from DEGlist)
## temp: remove iPSC stage
gene_counts_aftercpm_noiPSC.df <- gene_counts_aftercpm.df[,-1:-3]
data.dge.time <-  DGEList(gene_counts_aftercpm_noiPSC.df,group=rep(sample_date[-1],each=3))
annotations.df <- annotations.df[match(rownames(data.dge.time), rownames(annotations.df)),]
data.dge.time$genes <- annotations.df
data.dge.time <- calcNormFactors(data.dge.time)
data.dge.time$samples$norm.factors
# polynomial with 3 degrees of freedom
X <- poly(rep(sample_date[-1], each=3), degree=3)
design_time <- model.matrix(~X)
# estimating the NB dispersion
data.dge.time <- estimateDisp(data.dge.time,design_time)
sqrt(data.dge.time$common.dispersion)
plotBCV(data.dge.time)
# estimating the QL dispersion
fit <- glmQLFit(data.dge.time,design_time, robust=TRUE)
plotQLDisp(fit)
#--------------------
### Time course trend analysis (exclude iPSC stage)
##  looking for genes that change expression level over time
##  the design matrix uses 3 natural spline basis vectors to model smooth changes over time
##  test for a trend by conducting F-tests on 3 df for each gene
# run model
fit <-glmQLFTest(fit, coef=2:4)
# extract top set of genes with most significant time effects
tab <- as.data.frame(topTags(fit, n=30))
## visualize fitted spline curves fr the top four genes
# computing observed and fitted log-CPM value for each gene
logCPM.obs <- cpm(data.dge.time, log=TRUE, prior.count=fit$prior.count)
logCPM.fit <- cpm(fit, log=TRUE)
# loop through the first four genes and plot
par(mfrow=c(2,2))
for(i in 1:4) {
  FlybaseID <- row.names(tab)[i]
  Symbol <- tab$Gene_name[i]
  logCPM.obs.i <- logCPM.obs[FlybaseID,]
  logCPM.fit.i <- logCPM.fit[FlybaseID,]
  plot(rep(sample_date[-1], each=3), logCPM.obs.i, ylab="log-CPM", main=Symbol, pch=16, xlab="Days")
  lines(rep(sample_date[-1], each=3), logCPM.fit.i, col="red", lwd=2)}
# clear the figure setting: par(mfrow=c(1,1))
​
​
​
# ---------------------------------------------------------------------------
# Create contrast matrix
#make a contrast to test the effect of different developmental timepoints
cm <-makeContrasts( 
  iPS_HPC = (HPC - iPS),
  HPC_iMG = (D24 - HPC),
  M10_M5 = (M10 - M5),
  HPC_D8 = (D8 - HPC),
  D8_D16 = (D16 - D8),
  D16_D24 = (D24 -D16),
  D24_M10 = (M10 - D24),
  levels=design_matrix)
​
# add contrasts to the Lmfit object. 
fit <- contrasts.fit(fit=data.fit, contrasts=cm)
​
# perform bayes shrinkage to the analysis and estimate modified t and p values
fitDupCor <- eBayes(fit)
​
# --------------------------------------------------------------------------  
# find most differentially expressed genes
DGE_iPS_HPC <- topTable(fitDupCor, coef="iPS_HPC", n=80000)
DGE_HPC_iMGL <- topTable(fitDupCor, coef="HPC_iMG", n=80000)
DGE_iMGL_maturation <- topTable(fitDupCor, coef="M10_M5", n=80000)
# get all DEG
subset_iPS_HPC <- na.omit(subset(DGE_iPS_HPC, DGE_iPS_HPC$adj.P.Val<0.01 & abs(DGE_iPS_HPC$logFC)>=1))
subset_HPC_iMG <- na.omit(subset(DGE_HPC_iMGL, DGE_HPC_iMGL$adj.P.Val<0.01 & abs(DGE_HPC_iMGL$logFC)>=1))
subset_iMGL_maturation <- na.omit(subset(DGE_iMGL_maturation, DGE_iMGL_maturation$adj.P.Val<0.01 & abs(DGE_iMGL_maturation$logFC)>=1))
# get upregulated DEG
subset_iPS_HPC_P <- na.omit(subset(DGE_iPS_HPC, DGE_iPS_HPC$adj.P.Val<0.01 & DGE_iPS_HPC$logFC>=1))
subset_HPC_iMG_P <- na.omit(subset(DGE_HPC_iMGL, DGE_HPC_iMGL$adj.P.Val<0.01 & DGE_HPC_iMGL$logFC>=1))
subset_iMGL_maturation_P <- na.omit(subset(DGE_iMGL_maturation, DGE_iMGL_maturation$adj.P.Val<0.01 & DGE_iMGL_maturation$logFC>=1))
​
# Output upregulated DEG table
write.table(subset_iPS_HPC_P$Gene_name, "top_DE_iPStoHPC_genes_edgeR_U.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(subset_HPC_iMG_P$Gene_name, "top_DE_HPCtoiMGL_genes_edgeR_U.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(subset_iMGL_maturation_P$Gene_name, "top_DE_matureiMGL_genes_edgeR_U.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
​
​
# --------------------------------------------------------------------------
#comparison results with DESeq2
DESeq_iPS_HPC <- read.table(file = "/Users/cl2375/OneDrive - Yale University/PhD study/Labs/Brennand Lab/RNA-seq/RNAseq_iMGLpilot/DEG/DESeq2/cell stages/top_DE_iPStoHPC_genes_U.txt", head = FALSE)
DESeq_HPC_iMG <- read.table(file = "/Users/cl2375/OneDrive - Yale University/PhD study/Labs/Brennand Lab/RNA-seq/RNAseq_iMGLpilot/DEG/DESeq2/cell stages/top_DE_HPCtoiMGL_genes_U.txt", head = FALSE)
# intersact result from both DESeq2 and edgeR, extract edgeR result of union
DGE_iPS_HPC_both <- subset_iPS_HPC_P[rownames(subset_iPS_HPC_P) %in% intersect(subset_iPS_HPC_P$ensembl,DESeq_iPS_HPC$V1),]
DGE_HPC_iMGL_both <- subset_HPC_iMG[rownames(subset_HPC_iMG) %in% intersect(subset_HPC_iMG$ensembl,DESeq_HPC_iMG$V1),]
intersect(head(subset_iPS_HPC$ensembl,100),head(DESeq_iPS_HPC$V1,100))
​
# add annotations to the table
# DGE_iPS_HPC <- merge(DGE_iPS_HPC, annotations.df, by.x="row.names", by.y="ensembl", all=FALSE)
# DGE_HPC_iMGL <- merge(DGE_HPC_iMGL, annotations.df, by.x="row.names", by.y="ensembl")
​
# save DEGs
# write.table(DGE_NGN2_100nm, "DEG_response_100nM.txt", sep="\t")
# write.table(DGE_NGN2_1000nm, "DEG_response_1000.txt", sep="\t")
​
​
# do a heatmap!!
​
​
