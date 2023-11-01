## packages needed in this script
# install.packages("dplyr")
# install.packages("tibble")
# install.packages("stringr")
# install.packages("BiocManager")
# library(BiocManager)
# BiocManager::install("DESeq2")
### Bioconductor competibility with Macbook Air M1 https://support.bioconductor.org/p/9137290/
# install.packages("pheatmap")
# install.packages("svglite")
​
​
## load library
# library(DESeq2)
library(dplyr)
library(tibble)
library(stringr)
library(DESeq2)
library(pheatmap)
library(svglite)
​
## set working directory
setwd("/home/cl2375/project/RNAseq/")
​
## prepare input for DESeq2
countData <- read.table(file = "3182_filtered_genes_counts.txt", head = TRUE, sep = "\t", row.names = 1)
## clean up column & row name (row name specifically for DAVID)
colnames(countData) <- gsub("X", "", colnames(countData))
rownames(countData) <- str_split_fixed(rownames(countData), "\\.", 2)[,1]
​
## extract analysis targets
sampleNames <- c("iPSC", "HPC", "iMGL_D8", "iMGL_D16", "iMGL_D24", "iMGL_M5", "iMGL_M10")
​
sampleNames_iPSC <- c("3182_iPS_1", "3182_iPS_2", "3182_iPS_3")
sampleNames_HPC <- c("3182_HPC_1", "3182_HPC_2", "3182_HPC_3")
sampleNames_iMGL_D24 <- c("3182_D24_1", "3182_D24_2", "3182_D24_3")
sampleNames_iMGL_M10 <- c("3182_M10_1", "3182_M10_2", "3182_M10_3")
​
countData_iPSvHPC <- countData[,c(sampleNames_iPSC,sampleNames_HPC)]
countData_HPCviMGL <- countData[,c(sampleNames_HPC,sampleNames_iMGL_D24)]
countData_iMGLvmiMGL <- countData[,c(sampleNames_iMGL_D24,sampleNames_iMGL_M10)]
​
Conditions_iPSvHPC <- c("iPS","iPS","iPS","HPC","HPC","HPC")
Conditions_HPCviMGL <- c("HPC","HPC","HPC","iMGL","iMGL","iMGL")
Conditions_iMGLvmiMGL <- c("homeostatic","homeostatic","homeostatic","mature","mature","mature")
​
countDataTable_iPSvHPC <- data.frame(sampleName = c(sampleNames_iPSC,sampleNames_HPC), condition = Conditions_iPSvHPC)
countDataTable_HPCviMGL <- data.frame(sampleName = c(sampleNames_HPC,sampleNames_iMGL_D24), condition = Conditions_HPCviMGL)
countDataTable_iMGLvmiMGL <- data.frame(sampleName = c(sampleNames_iMGL_D24,sampleNames_iMGL_M10), condition = Conditions_iMGLvmiMGL)
​
## create DESeq dataset
DDS_iPSvHPC <- DESeqDataSetFromMatrix(countData = countData_iPSvHPC, colData = countDataTable_iPSvHPC, design = ~ condition)
DDS_HPCviMGL <- DESeqDataSetFromMatrix(countData = countData_HPCviMGL, colData = countDataTable_HPCviMGL, design = ~ condition)
DDS_iMGLvmiMGL <- DESeqDataSetFromMatrix(countData = countData_iMGLvmiMGL, colData = countDataTable_iMGLvmiMGL, design = ~ condition)
DDS_iPSvHPC$condition <- relevel(DDS_iPSvHPC$condition, ref = "iPS")
DDS_HPCviMGL$condition <- relevel(DDS_HPCviMGL$condition, ref = "HPC")
DDS_iMGLvmiMGL$condition <- relevel(DDS_iMGLvmiMGL$condition, ref = "homeostatic")
​
## Run DESeq
DDS_iPSvHPC <- DESeq(DDS_iPSvHPC)
DDS_HPCviMGL <- DESeq(DDS_HPCviMGL)
DDS_iMGLvmiMGL <- DESeq(DDS_iMGLvmiMGL)
​
## Retrieve results
res_iPSvHPC <- results(DDS_iPSvHPC, pAdjustMethod = "fdr")
res_iPSvHPC <- as.data.frame(res_iPSvHPC) %>% rownames_to_column("gene") %>% arrange(padj)
res_HPCviMGL <- results(DDS_HPCviMGL, pAdjustMethod = "fdr")
res_HPCviMGL <- as.data.frame(res_HPCviMGL) %>% rownames_to_column("gene") %>% arrange(padj)
res_iMGLvmiMGL <- results(DDS_iMGLvmiMGL, pAdjustMethod = "fdr")
res_iMGLvmiMGL <- as.data.frame(res_iMGLvmiMGL) %>% rownames_to_column("gene") %>% arrange(padj)
​
# Count adjusted p-value less than 0.01 and print it
print(sum(res_iPSvHPC$padj < 0.01, na.rm=TRUE))
print(sum(res_HPCviMGL$padj < 0.01, na.rm=TRUE))
print(sum(res_iMGLvmiMGL$padj < 0.01, na.rm=TRUE))
​
# Extract the significantly upregulated genes
top_DE_iPSvHPC <- na.omit(res_iPSvHPC [res_iPSvHPC $padj < 0.01 & res_iPSvHPC $log2FoldChange >= 1,])
top_DE_HPCviMGL <- na.omit(res_HPCviMGL [res_HPCviMGL $padj < 0.01 & res_HPCviMGL $log2FoldChange >= 1,])
top_DE_iMGLvmiMGL <- na.omit(res_iMGLvmiMGL [res_iMGLvmiMGL $padj < 0.01 & res_iMGLvmiMGL $log2FoldChange >= 1,])
​
# Output table
write.table(top_DE_iPSvHPC$gene, "top_DE_iPStoHPC_genes.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(top_DE_HPCviMGL$gene, "top_DE_HPCtoiMGL_genes.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(top_DE_iMGLvmiMGL$gene, "top_DE_matureiMGL_genes.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
​
# Extract only the top 50 significantly upregulated genes for heatmap generation
top_DE_iPSvHPC_50 <- top_DE_iPSvHPC[1:50,]
top_DE_HPCviMGL_50 <- top_DE_HPCviMGL[1:50,]
top_DE_iMGLvmiMGL_50 <- top_DE_iMGLvmiMGL[1:50,]
​
# Make a matrix of the counts for the top significantly upregulated genes
counts_top_DE_iPSvHPC_50 <- as.matrix(assay(DDS_iPSvHPC[intersect(top_DE_iPSvHPC_50$gene, rownames(assay(DDS_iPSvHPC))),]))
counts_top_DE_HPCviMGL_50 <- as.matrix(assay(DDS_HPCviMGL[intersect(top_DE_HPCviMGL_50$gene, rownames(assay(DDS_HPCviMGL))),]))
counts_top_DE_iMGLvmiMGL_50 <- as.matrix(assay(DDS_iMGLvmiMGL[intersect(top_DE_iMGLvmiMGL_50$gene, rownames(assay(DDS_iMGLvmiMGL))),]))
​
# Generate the heatmap
fig1 <- pheatmap(scale(counts_top_DE_iPSvHPC_50), fontsize_row = 9, main = "Heatmap of Top Differentially Expressed Genes from iPSC to HPC")
fig2 <- pheatmap(scale(counts_top_DE_HPCviMGL_50), fontsize_row = 9, main = "Heatmap of Top Differentially Expressed Genes from HPC to microglia")
fig3 <- pheatmap(scale(counts_top_DE_iMGLvmiMGL_50), fontsize_row = 9, main = "Heatmap of Top Differentially Expressed Genes from homeostatic to mature microglia")
​
# Save the heatmap into a pdf file
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(fig1, "iPSC_HPC.pdf")
save_pheatmap_pdf(fig2, "HPC_iMGL.pdf")
save_pheatmap_pdf(fig3, "iMGL_maturation.pdf")
​
