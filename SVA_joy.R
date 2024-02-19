#--------------------
# title:"SVA & RUV for batch effect"
# author: "Joy Lee"
# date: "9/28/2023"

#--------------------
### library
# basic operation
library(tidyverse)
# DEG analysis
library(limma)
library(edgeR)
# for plot
library(ggplot2)
library(ggfortify) # PCA
# for variance partition
library(variancePartition)
library(pheatmap)
# for SVA
library(sva)
# for RUV
library(RUVSeq)
library(RColorBrewer)




## function for variance partitioning
varPartAnalysis <- function(dge, formula) {
  
  ## Ignore genes with variance 0
  genes_var_zero <- which(apply(log(dge$counts+0.5), 1, var) == 0)
  if (length(genes_var_zero) > 0) {
    dge <- dge[-genes_var_zero, ]
  }
  ## Loop over each gene to fit model and extract variance explained by each variable
  varPart <- fitExtractVarPartModel(log(dge$counts+0.5), formula, dge$samples)
  # Sort variables by median fraction of variance explained
  vp <- sortCols(varPart)
  p <- plotVarPart(vp)
  return(list(p, vp))
}







#--------------------
### SVA to detact batch effect and unknown noise (method 2)
# extract counts table, use gene_counts_aftercpm.df / expression_matrix_res
# add all known variable back into metadata
meta.df$lib.size <- data.dge[["samples"]][["lib.size"]]/1000000
meta.df$norm.factors <- data.dge[["samples"]][["norm.factors"]]
# create two models: FULL and NULL
# make sure the scale of all variables are not too far off
mod <- model.matrix(~ group + batch + lib.size + norm.factors, meta.df)
mod0 <- model.matrix(~ batch + lib.size + norm.factors, meta.df)
## run sva, can try two methods 
n.sv = num.sv(gene_counts_aftercpm.df, mod, method="leek", B = 20)
head(n.sv)
# 5 surrogate variables. Now run sva to estimate the variable and make sure your counts is matrix
svobj <- svaseq(as.matrix(gene_counts_aftercpm.df), mod, mod0, n.sv=4) # WHY CAN I TRY 5???????
# include SV into data
svobj$sv <- as.data.frame(svobj$sv)
meta_sv.df <- as.data.frame(cbind(meta.df,svobj$sv))
## plot correlation matrix with surrogate variable to decide whether to include
# reset formula with SV
form <- ~ group + batch + lib.size + norm.factors + V1 + V2 + V3 + V4
# plot
C <- canCorPairs(form, meta_sv.df)
graphics.off()
pheatmap(
  C, ## data
  color = hcl.colors(50, "YlOrRd", rev = TRUE), ## color scale
  fontsize = 8, ## text size
  border_color = "black", ## border color for heatmap cells
  cellwidth = unit(0.4, "cm"), ## height of cells
  cellheight = unit(0.4, "cm") ## width of cells
)
## rerun variance partition to see if residuals are smaller
data.dge.sva <- data.dge
data.dge.sva$samples$V4 <- meta_sv.df$V4
form <- ~ (1|group) + (1|batch) + lib.size + norm.factors + V4
varPart <- varPartAnalysis(data.dge.sva, form)
varPart[[1]]



### Design matrix (with SVA)
# include all possible covariates in design, ie, plate, RIN, sex, etc
design_matrix_sva <- model.matrix(~ 0 + group + batch + V4, meta_sv.df)
# estimate the dispersion
data.dge.sva <- estimateDisp(data.dge.sva,design_matrix_sva)
### Voom transformation
## voom is a limma function that transforms the data into log2cpm based on normalization factors
data.sva.voom <- voom(counts=data.dge.sva, design=design_matrix_sva, plot=TRUE)
### Lmfit
# We have a general model that explains RNA expression and we want to fit our data to that model
data.sva.fit <- lmFit(object=data.sva.voom, design=design_matrix_sva) # limma package function
# check each coeficient in the model
head (coef(data.sva.fit))



### calculate the contrasts
# Create contrast matrix
#make a contrast to test the effect of different developmental timepoints
cm.sva <-makeContrasts( 
  TNF_on_all = (groupHCP_TNF + groupSCZ_TNF - groupHCP_NS - groupSCZ_NS),
  TNF_on_SCZ = (groupSCZ_TNF - groupSCZ_NS),
  TNF_on_HCP = (groupHCP_TNF - groupHCP_NS),
  SCZ_baseline = (groupSCZ_NS - groupHCP_NS),
  levels=design_matrix_sva)
# add contrasts to the Lmfit object. 
fit.sva <- contrasts.fit(fit=data.sva.fit, contrasts=cm.sva)
# perform bayes shrinkage to the analysis and estimate modified t and p values
fitDupCor.sva <- eBayes(fit.sva)



### find most differentially expressed genes
DGE_TNF_sva <- topTable(fitDupCor.sva, coef="TNF_on_all", n=80000)
DGE_TNF_HCP_sva <- topTable(fitDupCor.sva, coef="TNF_on_HCP", n=80000)
DGE_TNF_SCZ_sva <- topTable(fitDupCor.sva, coef="TNF_on_SCZ", n=80000)
DGE_SCZ_sva <- topTable(fitDupCor.sva, coef="SCZ_baseline", n=80000)
subset_HCPvsSCZ_sva <- na.omit(subset(DGE_SCZ_sva, DGE_SCZ_sva$adj.P.Val<0.05 & abs(DGE_SCZ_sva$logFC)>=1))
subset_TNFonHCP_sva <- na.omit(subset(DGE_TNF_HCP_sva, DGE_TNF_HCP_sva$adj.P.Val<0.05 & abs(DGE_TNF_HCP_sva$logFC)>=1))
subset_TNFonSCZ_sva <- na.omit(subset(DGE_TNF_SCZ_sva, DGE_TNF_SCZ_sva$adj.P.Val<0.05 & abs(DGE_TNF_SCZ_sva$logFC)>=1))
subset_TNFonALL_sva <- na.omit(subset(DGE_TNF_sva, DGE_TNF_sva$adj.P.Val<0.05 & abs(DGE_TNF_sva$logFC)>=1))

write.table(subset_TNFonHCP_sva, "topDEG_TNFonHCP_sva.txt", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(subset_TNFonSCZ_sva, "topDEG_TNFonSCZ_sva.txt", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(subset_TNFonALL_sva, "topDEG_TNFonALL_sva.txt", quote = FALSE, col.names = TRUE, row.names = FALSE)






#--------------------
### RUV using cpm data
## plottting to look at data
# plot the data to look at distribution
x <- as.factor(meta.df$group)
set <- newSeqExpressionSet(as.matrix(gene_counts_aftercpm.df),
                           phenoData = data.frame(meta.df, row.names=rownames(meta.df)))
colors <- brewer.pal(3, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x], las=2)
plotPCA(set, col=colors[x], cex=1.2)
# normalized data using UQ normalization
set <- betweenLaneNormalization(set, which="upper")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)



## RUV with 5 variables
## RUVs - method 1 with replicates
# contruct matrix specifying the replicates
differences <- makeGroups(x)
genes <- rownames(counts(set))
set1 <- RUVs(set, genes, k=5, differences)
## check correlation between SUV variables and known variable
# include W data into metadata
meta.df$lib.size <- data.dge[["samples"]][["lib.size"]]
meta.df$norm.factors <- data.dge[["samples"]][["norm.factors"]]
meta_ruv.df <- meta.df
meta_ruv.df$W1 <- set1$W_1
meta_ruv.df$W2 <- set1$W_2
meta_ruv.df$W3 <- set1$W_3
meta_ruv.df$W4 <- set1$W_4
meta_ruv.df$W5 <- set1$W_5
# model with RUV variables
form <- ~ group + batch + lib.size + norm.factors + W1 + W2 + W3 + W4 + W5
# plot
C <- canCorPairs(form, meta_ruv.df)
graphics.off()
pheatmap(
  C, ## data
  color = hcl.colors(50, "YlOrRd", rev = TRUE), ## color scale
  fontsize = 8, ## text size
  border_color = "black", ## border color for heatmap cells
  cellwidth = unit(0.4, "cm"), ## height of cells
  cellheight = unit(0.4, "cm") ## width of cells
)
# W1 = batch


## rerun variance partition to see if residuals are smaller
data.dge.ruv <- data.dge
data.dge.ruv$samples$W2 <- meta_ruv.df$W2
data.dge.ruv$samples$W3 <- meta_ruv.df$W3
data.dge.ruv$samples$W4 <- meta_ruv.df$W4
data.dge.ruv$samples$W5 <- meta_ruv.df$W5
form <- ~ (1|group) + (1|batch) + lib.size + norm.factors + W2 + W3 + W4 + W5
varPart <- varPartAnalysis(data.dge.ruv, form)
varPart[[1]]


### Design matrix (with RUV-method1)
design_matrix_ruv <- model.matrix(~ 0 + group + batch + W2 + W3 + W4 + W5, meta_ruv.df)
# estimate the dispersion
data.dge.ruv <- estimateDisp(data.dge.ruv,design_matrix_ruv)
### Voom transformation
## voom is a limma function that transforms the data into log2cpm based on normalization factors
data.ruv.voom <- voom(counts=data.dge.ruv, design=design_matrix_ruv, plot=TRUE)
### Lmfit
# We have a general model that explains RNA expression and we want to fit our data to that model
data.ruv.fit <- lmFit(object=data.ruv.voom, design=design_matrix_ruv) # limma package function
# check each coeficient in the model
head (coef(data.ruv.fit))



### calculate the contrasts
# Create contrast matrix
#make a contrast to test the effect of different developmental timepoints
cm.ruv <-makeContrasts( 
  TNF_on_all = (groupHCP_TNF + groupSCZ_TNF - groupHCP_NS - groupSCZ_NS),
  TNF_on_SCZ = (groupSCZ_TNF - groupSCZ_NS),
  TNF_on_HCP = (groupHCP_TNF - groupHCP_NS),
  SCZ_baseline = (groupSCZ_NS - groupHCP_NS),
  levels=design_matrix_ruv)
# add contrasts to the Lmfit object. 
fit.ruv <- contrasts.fit(fit=data.ruv.fit, contrasts=cm.ruv)
# perform bayes shrinkage to the analysis and estimate modified t and p values
fitDupCor.ruv <- eBayes(fit.ruv)



### find most differentially expressed genes
DGE_TNF_ruv <- topTable(fitDupCor.ruv, coef="TNF_on_all", n=80000)
DGE_TNF_HCP_ruv <- topTable(fitDupCor.ruv, coef="TNF_on_HCP", n=80000)
DGE_TNF_SCZ_ruv <- topTable(fitDupCor.ruv, coef="TNF_on_SCZ", n=80000)
DGE_SCZ_ruv <- topTable(fitDupCor.ruv, coef="SCZ_baseline", n=80000)
subset_HCPvsSCZ_ruv <- na.omit(subset(DGE_SCZ_ruv, DGE_SCZ_ruv$adj.P.Val<0.05 & abs(DGE_SCZ_ruv$logFC)>=1))
subset_TNFonHCP_ruv <- na.omit(subset(DGE_TNF_HCP_ruv, DGE_TNF_HCP_ruv$adj.P.Val<0.05 & abs(DGE_TNF_HCP_ruv$logFC)>=1))
subset_TNFonSCZ_ruv <- na.omit(subset(DGE_TNF_SCZ_ruv, DGE_TNF_SCZ_ruv$adj.P.Val<0.05 & abs(DGE_TNF_SCZ_ruv$logFC)>=1))
subset_TNFonALL_ruv <- na.omit(subset(DGE_TNF_ruv, DGE_TNF_ruv$adj.P.Val<0.05 & abs(DGE_TNF_ruv$logFC)>=1))


## RUVr
design_matrix_ruv <- model.matrix(~ 0 + group + batch, meta.df)



mod <- model.matrix(~ group + batch + lib.size + norm.factors, meta.df)
