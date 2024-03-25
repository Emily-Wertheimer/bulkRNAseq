### Running fgsea ###

## Set up
library(fgsea)
library(ggplot2)
library(clusterProfiler)
library(DOSE)
library(dplyr)

setwd("~/Library/CloudStorage/OneDrive-YaleUniversity/Lab/Analyses/20230505_hormones_RNAseq/DE_ds1_alcohol_analysis")

## Starting from limma res list
res.list <- readRDS("scripts/alcohol_RNAseq_reslist.rds") # This is the results object from limma


## Reading in Nadine's gene set files
files <- list.files(path = "./genesets", pattern = ".gmt")

## Running analysis and saving dataframes to txt, also saving significant pathways figure

for (j in 1:length(files)) {
  fgsea.res.list <- list()
  gmt <- read.gmt(paste0("genesets/", files[j]))
  gmt.set <- files[[j]]
  for (i in 1:length(res.list)) { # b/c my res list contains many contrasts, looping this
    pathways <- gmtPathways(paste0("genesets/", gmt.set))
    str(head(pathways))
    
    ranks <- res.list[[i]][,c("hgnc_symbol","logFC")]
    ranks <- setNames(ranks$logFC, ranks$hgnc_symbol)
    str(ranks)
    
    fgsea.res.list[[i]] <- fgsea(pathways,
                                 ranks)
    fgsea.res.list[[i]] <- fgsea.res.list[[i]][order(pval), ]
  }
  names(fgsea.res.list) <- names(res.list)
  
  p <- list()
  for (i in 1:length(fgsea.res.list)) {
    name <- names(fgsea.res.list)[i]
    fgseaRes <- fgsea.res.list[[i]]
    data.table::fwrite(fgseaRes, file = paste0("results/", gmt.set, "_", name, "_fgsea.txt"), sep = "\t", sep2 = c("", " ", ""))
    topPathwaysUp <- fgseaRes[ES > 0][padj <= 0.05, pathway]
    topPathwaysDown <- fgseaRes[ES < 0][padj <= 0.05, pathway]
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    p[[i]] <- plotGseaTable(pathways[topPathways], ranks, fgseaRes, 
                            gseaParam=0.5)
  }
  
  pdf(paste0("figures/", gmt.set, "_", name, "_fgsea.pdf"))
  for (i in 1:length(p)) {
    print(p[[i]])
  }
  dev.off()
}
