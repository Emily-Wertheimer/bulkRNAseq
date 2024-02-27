#tenemos los archivos de arrays en zip
#los cargamos en partek flow de 2 en 2 (for y rev) en un proyecto nuevo [import data_transform automatically]
# create samples

#analyses --> unsigned reads --> prealignment QA/QC
#unsigned reads-->aligners-->star (mouse)

#aligner reads--> quantification
#quantify to annotation model (Partek E/M)
#Gene/Feature annotation (refseq transcript 98-2021-05-05)
#DESCARGAR LOS ARCHIVOS--> GENE COUNTS


####PODEMOS HACER EL ANALISIS GLOBAL, PERO LO VAMOS A HACER PRIMERO POR PARES
#PARA ELLO: unzip --> abrimos excel --> seleccionamos las filas de los que necesitamos y lo trasponermos en otro excel (matrix)--> guardar CSV
#PONEMOS LOS GRUPOS en otro excel (GROUPS)--> GUARDAR CSV
#files(more:set as working directory)

###cuando sale de partekflow 

library(readr)
library(readxl)
library(tibble)
#First, we upload the data matrix with the genes and the samples#
matrix_Tmem86b<-read.csv("Tmem86b_PCSK9_WD_Quantify_to_annotation_model_Gene_counts.csv",header=TRUE, sep=",")
View(matrix_Tmem86b)

matrix_Tmem86b<-as.data.frame(matrix_Tmem86b)
row.names(matrix_Tmem86b)<-matrix_Tmem86b$Gene_Symbol
matrix_Tmem86b<-matrix_Tmem86b[-1]
head(matrix_Tmem86b)

#Then, we upload the data matrix with the information of the samples (genotype,...). This is the group list#
group_list<-read.csv("GROUPS_WD.csv",header = TRUE,sep=",")
row.names(group_list)<-group_list$X
group_list<-group_list[-1]
group_list$Genotype <- as.factor(group_list$Genotype)
str(group_list)



#####If we don't have a matrix with the information of the sample we can do this. This is even good for future actions with sample names.#####
# Get metadata
#SampleID<-gsub("(.*)\\()",colnames(matrix_Tmem86b))
#GroupVar<-gsub("(.*)\\_(.*)\\_(.*)\\_(.*)\\_(.*)\\_(.*)","\\1",colnames(matrix_Tmem86b))
GroupVar<-gsub("(.*?)\\d(.*)","\\1",colnames(matrix_Tmem86b))
#Replicate<-gsub("(.*)\\_(.*)\\_(.*)\\_(.*)\\_(.*)\\_(.*)","\\2",colnames(matrix_Tmem86b))
Replicate<-gsub("(.*)\\_(.*)","\\2",colnames(matrix_Tmem86b))

##Change names
GroupVar<-sub("TKO","Tmem86b",GroupVar)
GroupVar<-sub("WT","WT",GroupVar)

colnames(matrix_Tmem86b)<-paste(GroupVar,Replicate,sep="_")
metadata<-data.frame(cbind(sample_id=colnames(matrix_Tmem86b),GroupVar,Replicate))

log_traza(metadata)

# Reorder sampleID's to match featureCounts column order. 
metadata <- metadata[match(colnames(matrix_Tmem86b), metadata$sample_id), ]

# Fix this to set logFC by default order 2/1 in this case DEGS1/SCRAM
metadata$GroupVar<-factor(metadata$GroupVar, levels=c("WT","Tmem86b")) 




#Now, we start to build the matrix and the package for the DESeq2 analysis#
library(DESeq2)
#If this does not work:
#install.packages("DESeq2")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")


######This is one option to build the matrix for the DESeq2 analysis######
# - countData : count dataframe
# - colData : sample metadata in the dataframe with row names as sampleID's
# - design : The design of the comparisons to use. 
#            Use (~) before the name of the column variable to compare
ddsMat <- DESeqDataSetFromMatrix(countData = round(matrix_Tmem86b), ## non-negative integer only 
                                 colData = metadata,
                                 design = ~GroupVar)


######This is another option  to do it######


#Contruimos la matriz para DESeq2
DEXTmem86bWD<-DESeqDataSetFromMatrix(countData = round(matrix_Tmem86b), colData = group_list, design = ~Genotype)
head(DEXTmem86bWD)



#Remove some conditions TODO: make options

unique(ddsMat$sample_id)
unique(ddsMat$GroupVar)

#dds<-dds[,which(!dds$sample_id %in% c("ChowAnimal4","ChowAnimal5","ChowAnimal6"))]

## Remove condition which is outlier clearly. I've seen which ones quit in the following representations.

ddsMat_WO_O<-ddsMat[,which(!ddsMat$sample_id %in% c("WT_S5","Tmem86b_S10"))] # Pero en analisis de foldchange incluir
ddsMat_WO_O02<- ddsMat[,which(!ddsMat$sample_id %in% c("WT_S4","Tmem86b_S10"))]
# Transform counts for data visualization
rld <- rlog(ddsMat, blind=TRUE)
rld02 <- rlog(ddsMat_WO_O, blind = TRUE)
rld03 <- rlog(ddsMat_WO_O02, blind = TRUE)

# Plot PCA
colData(ddsMat)
DESeq2::plotPCA(rld, ntop = 500, intgroup = "GroupVar")
DESeq2::plotPCA(rld, ntop = 500, intgroup = "Replicate")
DESeq2::plotPCA(rld, ntop = 500, intgroup = "sample_id")

DESeq2::plotPCA(rld02, ntop = 500, intgroup = "GroupVar")
DESeq2::plotPCA(rld02, ntop = 500, intgroup = "Replicate")
DESeq2::plotPCA(rld02, ntop = 500, intgroup = "sample_id")

DESeq2::plotPCA(rld03, ntop = 500, intgroup = "GroupVar")
DESeq2::plotPCA(rld03, ntop = 500, intgroup = "Replicate")
DESeq2::plotPCA(rld03, ntop = 500, intgroup = "sample_id")

QC_pca_plot<- DESeq2::plotPCA(rld03, ntop = 500, intgroup = "GroupVar",)+ 
  geom_point(size=5,shape=21,color="black")+
  geom_text(label=rld03$sample_id,check_overlap = T,
            nudge_x = 1, nudge_y =1,
            color="black",
            size=3)+
  theme_classic()

#ggtitle(paste0(names(dds_ls)[i],", cells= ",sum(rld$cell_count)))+
#lims(y=c(-40,40)) +
#lims(x=c(-40,40)) +


# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

rld_mat02 <- assay(rld02)
rld_cor02 <- cor(rld_mat02)

rld_mat03 <- assay(rld03)
rld_cor03 <- cor(rld_mat03)

# Plot heatmap
rownames(metadata)<-metadata$sample_id
pheatmap(rld_cor03, annotation = metadata[, c("GroupVar"), drop=F])
#rownames(metadata)<-paste0(metadata$GroupVar,"_",metadata$Replicate)
#pheatmap(rld_cor, annotation = rownames(metadata))

#Pulimos la matriz quitando los genes que se hayan detectado con 10 o menos cuentas
keep_HExpG <- rowSums(counts(ddsMat_WO_O02)) >=10

ddsMat_WO_O_H02 <- ddsMat_WO_O02[keep_HExpG,]


# Run DESeq2 differential expression analysis
#dds <- DESeq(dds)

dds_diff<-DESeq(ddsMat_WO_O_H02,test = 'Wald',
           fitType = 'parametric',
           sfType = 'poscounts',
           betaPrior = F)
# Plot dispersion estimates
plotDispEsts(dds_diff)

# Check the coefficients for the comparison
resultsNames(dds_diff)

# Generate results object
resultados <- results(dds_diff, 
               name = resultsNames(dds)[2],
               alpha = 0.05)

resultados

# Shrink the log2 fold changes to be more appropriate using the apeglm method - should cite [paper]() when using this method
resultados <- lfcShrink(dds_diff, 
                 coef = resultsNames(dds)[2],
                 res=resultados,
                 type = "apeglm")


res_tbl_diff <- resultados %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  arrange(padj)

##Force Scientific notation in some columns
# res_tbl$pvalue<-formatC(res_tbl$pvalue, format = "e")
# res_tbl$padj<-formatC(res_tbl$padj, format = "e") 



###--------------        Annotate results      -----------------
## -------------------------------------------------------------
library(org.Mm.eg.db)
library(org.Hs.eg.db)
#results_annotated<-res_tbl
rownames(res_tbl_diff)<-res_tbl_diff$gene
# Add gene full name
res_tbl_diff$description <- mapIds(x = org.Mm.eg.db,
                              keys = rownames(res_tbl_diff),
                              column = "GENENAME",
                              keytype = "SYMBOL",
                              multiVals = "first")

# Add gene symbol
#res_tbl$symbol <- res_tbl$gene
rownames(res_tbl_diff)<-res_tbl_diff$gene
# Add ENTREZ ID
res_tbl_diff$entrez <- mapIds(x = org.Mm.eg.db,
                         keys = rownames(res_tbl_diff),
                         column = "ENTREZID",
                         keytype = "SYMBOL",
                         multiVals = "first")

# Add ENSEMBL
rownames(res_tbl_diff)<-res_tbl_diff$gene
res_tbl_diff$ensembl <- mapIds(x = org.Mm.eg.db,
                          keys = row.names(res_tbl_diff),
                          column = "ENSEMBL",
                          keytype = "SYMBOL",
                          multiVals = "first")
#Remove any genes that do not have any entrez identifiers

res_tbl_diff <- subset(res_tbl_diff, is.na(entrez) == FALSE)

res_tbl_diff.na<-res_tbl_diff[is.na(res_tbl_diff$log2FoldChange),"gene"]
na_list_diff<-res_tbl_diff[is.na(res_tbl_diff$log2FoldChange),"gene"]
#normalized_counts.na<-normalized_counts[which(rownames(normalized_counts) %in% na_list$gene),]

## Remove genes with zero=NA form results

res_tbl_diff <- subset(res_tbl_diff, is.na(log2FoldChange) == FALSE)



###Normalize and write the counts###
normCounts_dds_diff<-counts(dds_diff, normalized=TRUE)

write.csv(normCounts_dds_diff,"DEXTmem86bWD_NormCounts_dds_diff.csv")


write.csv(resultados, "DEXTmem86bWD_results_dds.csv")
write.csv(res_tbl_diff, "DEXTmem86bWD_results_dds_complete_diff.csv")


plotCounts(dds, gene = "Tmem164",intgroup = "GroupVar")

plotCounts(dds, gene = "Il6ra",intgroup = "GroupVar")

MAplot_res <- plotMA(dds, alpha = 0.05, ylim = c(-6,6))


summary(resultados)




#Pulimos la matriz quitando los genes que se hayan detectado con 10 o menos cuentas
keep_HExpG <- rowSums(counts(DEXTmem86bWD)) >=10

DEXTmem86bWD_H <- DEXTmem86bWD[keep_HExpG,]

DEXTmem86bWD<-DESeq(DEXTmem86bWD_H)

normCounts<-counts(DEXTmem86bWD, normalized=TRUE)

write.csv(normCounts,"DEXTmem86bWD_NormCounts.csv")

resultsNames(DEXTmem86bWD)

res<- results(DEXTmem86bWD)

res02 <- results(DEXTmem86bWD, contrast = list("Genotype_WT_vs_KO"))
res02_padj05 <- results(DEXTmem86bWD, contrast = list("Genotype_WT_vs_KO"), alpha=0.05)

res02 <- results(DEXTmem86bWD, contrast = c("Genotype","KO","WT"))
res02_padj05 <- results(DEXTmem86bWD, contrast = c("Genotype","KO","WT"), alpha=0.05)


summary(res02)
summary(res02_padj05)

write.csv(res02, "DEXTmem86bWD_results02.csv")
write.csv(res02_padj05, "DEXTmem86bWD_results_a0.0502.csv")

DEXTmem86bWD_log <- rlog(DEXTmem86bWD)
#When we have too many samples (more than 50) we can use vst() function, which is much faster#
DEXDEXTmem86b_vst <- vst(DEXTmem86b)

plotPCA(DEXTmem86bWD_log, intgroup = "Genotype")


DEXTmem86bWD_var <- varianceStabilizingTransformation(DEXTmem86bWD)
plotPCA(DEXTmem86bWD_var, intgroup = "Genotype")


MAplot_res <- plotMA(DEXTmem86bWD, alpha = 0.05, ylim = c(-6,6))
MAplot_res2 <- plotMA(DEXTmem86bWD, contrast = list("Genotype"), alpha= 0.05, ylim = c(-6,6))

plotCounts(DEXTmem86bWD, gene = "Il6",intgroup = "Genotype")


d2 <- normCounts_dds


#Para representar heatmaps. TRAS SCRIPT PASO 1


library(pheatmap)
library(edgeR)
install.packages("edgeR")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")


cpm_converted_DEXTmem86bWD<-cpm(dds, log = TRUE)
head(cpm_converted_DEXTmem86bWD)
n<-t(scale(t(cpm_converted_DEXTmem86bWD)))
x<-scale(cpm_converted_DEXTmem86bWD)
n[n>2]=2
n[n<-2]=-2
n[1:4,1:4]

#HEATMAPS

#fibrosis new all
c<-n[row.names(n)%in%c("Fn1","Lama5","Lamb1","Lamb2","Lamc1","Lamc2","Vtn","Thbs1","Thbs3","Col3a1","Col1a1","Col1a2","Acta1","Acta2","Tgfb1","Tgfb2","Tgfb3","Tgfa","Pdgfa","Pdgfb","Pdgfc","Pdgfd","Map4","Lama1","Lama2","Lama3","Lama4","Lamb3","Lamb4","Lamc3","Lamc4","Egflam","Lamb2p1","Col17a1","Ctgf","Vegfa","Vegfb","Vegfc","Vegfd","Mmp2","Mmp9","Igf1","Timp1","Timp2","Timp3","Mmp13","Mmp12"),]
#fibrosis new sig
c<-n[row.names(n)%in%c("Fn1","Lama4","Lamb2","Mmp2","Pdgfa","Pdgfb","Tgfb1","Tgfb2","Tgfb3","Thbs1","Timp2","Vtn","Krt18","Fgb","Fga","Fgg"),]
#fibrosis new sig
c<-n[row.names(n)%in%c("Mmp2","Pdgfb","Tgfb2","Pdgfa","Thbs1","Timp2","Lama4","Tgfb1","Lamb2","Tgfb3","Fga","Fgb","Fgg","Krt18","Fn1","Vtn"),]

#inflammation sig
c<-n[row.names(n)%in%c("Cd101", "Itgb3", "Aif1","Aim2","Aimp1","Alox15","Aoah","Arel1","C1qa","C3ar1","Casp12","Ccl24","Cd200r3","Cd5l","Cebpa","Cebpb","Celf1","Crp","Csf1","Csrp3","Cxcl1","Cxcl13","Dhx9","Ets1","F12","F7","Fabp4","Fgfr1","Fpr1","Grn","Hdac5","Hps1","Ighg1","Ighg2b","Il1a","Il1b","Il1r1","Irf5","Itga2","Itgam","Itgb6","Lbp","Lilra5","Mapk14","Metrnl","Naip1","Naip2","Naip5","Ndfip1","Nfkbia","Nlrx1","Nr1d1","Nr1h4","Nrros","Pdcd4","Per1","Pik3ap1","Pik3cd","Plaa","Pld3","Psen2","Ptgs1","Rabgef1","Rarres2","Rora","Saa3","Serpinb1a","Siglece","Sirpa","Smad3","Socs3","Stap1","Sucnr1","Ticam2","Tlr4","Tlr6","Tnf","Xcl1","Ywhaz","Zfp36"),]
#inflammation sig
c<-n[row.names(n)%in%c("Cd86","Cxcl1","Il1b","Itgam","Ltb","Mmp9","Socs3","Tnf"),]

#Apoptosis sig
c<-n[row.names(n)%in%c("Bcl2l2","Bnip3l","Casp2","Cflar","Cycs","Ikbkg","Irf1","Irf2","Irf5","Map3k1","Myc","Nfkbia","Tnf","Tnfrsf25","Tnfsf10"),]

#Lipogenesis sig
c<-n[row.names(n)%in%c("Acacb","Scd1","Acot1","Pecr","Hsd17b10","Echdc3","Acss2","Hadh","Mogat1","Tecr","Acsl5","Acot2","Acsl4","Hsd17b6","Ech1","Echs1","Hsd17b4","Decr1","Srebf1"),]

#FAO sig
c<-n[row.names(n)%in%c("Acox1","Ehhadh","Acss2","Hadh","Lpl","Acad9","Acadl","Gcdh","Ivd","Acad10","Acadvl","Hadha","Acsl5","Hadhb","Acsl4","Cpt1a","Ech1","Rsad1","Aldh2","Echs1","Hsd17b4","Acadm","Acad8","Decr1","Chkb","Slc25a20","Ppargc1a"),]

#Chol metab sig
c<-n[row.names(n)%in%c("Hmgcr","Lss","Fdps","Cyp51","Dhcr7","Sc5d","Mvd","Srebf2"),]

#Electron transport chain sig
c<-n[row.names(n)%in%c("Ndufa8","Atpif1","Ndufab1","Uqcr10","Sdhb","Ndufs4","Cox7b","Sdha"),]

#Lipoprotein metabolism sig
c<-n[row.names(n)%in%c("Lipg","Lipc","Mttp","Apoa1","Apoc1","Lpl","Apoa4","Apob","Apoe","Apoc4","Apoa5"),]

#FA transporters sig
c<-n[row.names(n)%in%c("Fabp1","Dbi","Acsl5","Acsl4","Cd36","Slc27a4","Fabp2","Slc27a5","Slc27a2"),]

#Plasmalogen metab sig
c<-n[row.names(n)%in%c("Tmem86b","Far1","Pla2g15","Pla2g6","Pla2g7","Plcg2","Selenoi","Plpp1"),]

#OX. Stress sig
c<-n[row.names(n)%in%c("Hmox1","Mt1","Fos","Cyp4a14","Nfe2l2","Mgst1","Cyp2b9","Gclc","Gsr","Cyba","Gpx1","Mapk14","Nfkb1","Sod2","Nqo1","Junb","Txnrd2"),]

#ER. Stress sig
c<-n[row.names(n)%in%c("Creb3","Creb3l3","Hspa5","Crebrf","Eif2s2","Atf6"),]

#EIF2 response genes sig
c<-n[row.names(n)%in%c("Myof","Iqgap1","Peg3","Dock7","Myo1b","Atf5","Mthfd1","Filip1l","Aldh18a1","Eprs","Stat3","Hadha","Safb","Crybg3","Hsph1","Aars","Top2a","Gbf1","Ubr2","Hectd1","Trip12","Ganab","Sec23ip","Helb","Dhx9","Pi4ka","Mical1","Nup155","Polr2a","Ylpm1","Ap3b1","Dhx36","Baz1b","Eif3a","Ttc17","Tanc2","Cdk5rap2","Ubr1","Snd1"),]


#Acute phase response POSITIVE regulators genes sig
c<-n[row.names(n)%in%c("Crp","Orm1","Lbp","Serpina1a","C3","C4b","Hamp","Fga","Fgb","Fgg","Vtn","C6"),]

#Acute phase response NEGATIVE regulators genes sig
c<-n[row.names(n)%in%c("Alb","Ttr","Trf","Serpina6"),]



pheatmap(c, cluster_cols = F, cluster_rows = T,  border_color = 'black',
         main = 'APR Neg_GENES', fontsize_row = 8, fontsize_col = 6)




####volcano######

install.packages("EnhancedVolcano")
library(EnhancedVolcano)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")

library(ggplot2)
library(ggrepel)


dds2 <- DESeq(dds, betaPrior=FALSE)
res2 <- results(DEXTmem86bWD,
                contrast = c('dex','trt','untrt'))
res2 <- lfcShrink(DEXTmem86bWD,
                  contrast = c('dex','trt','untrt'), res=res, type = 'normal')


topT <- as.data.frame(res_tbl)
rownames(topT) <- topT$gene
?EnhancedVolcano
EnhancedVolcano(topT,
                lab = rownames(topT),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                ylim = c(0,30),
                xlim = c(-5,5),
                col = c('grey', 'lightskyblue3','lightblue','red'))





#######################GSEA analysis###################

library(tidyverse)
library(msigdbr)
library(fgsea)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(logr)

##We load the data to analyze and de databases we want to compare with
normCounts02 <- read.csv("DEXTmem86bWD_NormCounts_WO_WT5KO9_copy.csv")
normCounts_dds02 <- read.csv("DEXTmem86bWD_NormCounts_dds.csv")
row.names(normCounts02)<-normCounts02$X
normCounts02<-normCounts02[-1]

row.names(normCounts_dds02)<-normCounts_dds02$X
normCounts_dds02<-normCounts_dds02[-1]

H<-msigdbr(species = "Mus musculus", category = "H")
M<-msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP")
M.a<-msigdbr(species = "Mus musculus", category = "C2")
class(H)
head(H)

write.csv(M,"Mouse_gene_dataset_M.csv")
write.csv(M.a,"Mouse_gene_dataset_M.a.csv")
write.csv(Mo,"Mouse_gene_dataset_Mo.csv")

### Check some duplicated genes, genes not annotated, without ensembl

ifelse(sum(duplicated(as.character(res_tb$gene)))>0,"CHECK FOR DUPLICATED GENES !!","GREAT, NOT DUPLICATED GENES !!")
ifelse(sum(duplicated(as.character(res_tb$ensembl)))>0,"CHECK FOR DUPLICATED Ensmbl GENES !!","GREAT, NOT DUPLICATED GENES !!")

#look for genes without ensembl, Check lattest version of org.Mm.eg.db
model.results.NA<-res_tbl[is.na(res_tbl$ensembl),]
# genes starting with 'Gm-' or ending with 'Rik' include protein-coding RNAs, long non-coding RNAs, and antisense transcripts.
model.results.NA_not_RikGm<-model.results.NA[!grepl("Rik|Gm",model.results.NA$gene),]
#length(model.results.NA_not_RikGm$gene)
## Number of genes Rik or Gm
model.results.Rik_Gm<-res_tbl[grepl("Rik|Gm",res_tbl$gene),]

log_print(paste0("Number of genes without annotation in ensembl = ",length(model.results.NA$gene)))
log_print(paste0("Not annotated and not .*Rik or Gm.* = ",length(model.results.NA_not_RikGm$gene)))
log_print(paste0("total.*Rik or Gm.* in model.results = ",length(model.results.Rik_Gm$gene)))


### Select Molecular Signature to Explore
MolSignatures<- c("H","C1","C2","C3","C4","C5","C6","C7","C8")
MolSignatures<- c("H","C1","C2")
names(MolSignatures)<-c("Hallmark",
                        "Positional",
                        "Curated gen sets",
                        "Regulatory",
                        "Computational",
                        "Ontologies",
                        "Oncogenic",
                        "Immunologic",
                        "Cell type sets"
)

names(MolSignatures)<-c("Hallmark",
                        "Positional",
                        "Curated gen sets"
                        
)


MolSignatures <-selectOption(MolSignatures,"Select Molecular Signatures to expore")
#M <- msigdbr(species = "Mus musculus", category = "C2", subcategory = c("CP:KEGG"))
Mo <-bind_rows(lapply(MolSignatures,FUN=function(x){msigdbr(species = "Mus musculus", category = x)}),.id = "column_label")[,-c(1)] 
Mo$gs_subcat[Mo$gs_subcat==""]<-"no subcat"

# Select which Subcats
MolSubcats<-unique(Mo$gs_subcat)
MolSubcats <-selectOption(MolSubcats,"Select Molecular Sub-Categories (TO REMOVE)",toRemove = T)
Mo <-Mo[which(Mo$gs_subcat %in% MolSubcats),]

print(unique(Mo[,c(1,2)]))

## Filter specific description terms

# repeat{
#   statements...
#   if(condition){
#     break
#   }
# }

list_MolTerms<- c("List of Terms to filter separated by ,"="phospholip,plasmalog,collagen,stellate,macrophage,monocyte,lipase,lipoprotein,cholesterol,inflamma,fibro,apopto,ROS,oxidat,mitochon,fatty,endoplasmic_reticulum_st,NFKB,ECM")
fix(list_MolTerms)

list_MolTerms<-str_split_1(list_MolTerms,pattern=",")


log_traza(list_MolTerms)

MolTerms<-paste0(".*(",
                 paste0(paste(paste0(toupper(gsub("(.).*","[\\1|",list_MolTerms)),gsub("(.).*","\\1]",list_MolTerms),gsub("(.)(.*)","\\2",list_MolTerms)),collapse="|"),
                        #paste(toupper(paste0(toupper(gsub("(.).*","[\\1|",list_MolTerms)),gsub("(.).*","\\1]",list_MolTerms),gsub("(.)(.*)","\\2",list_MolTerms))),collapse="|")
                        #),
                        ")"))

MolTerms
#MolTerms<-".*([S|s]phingo|[C|c]eramide|[M|m]acrophage|[M|m]onocyte|[S|s]tellate|[C|c]ollagen)|[L|l]ipase|[L|l]ipoprotein|[C|c]holesterol)"
#names(MolTerms)<-"Specific name terms to filter by: * (if all)"
#fix(MolTerms)

M<-Mo[grepl(MolTerms,Mo$gs_description),]
M_pathways_selected<-M %>% dplyr:: mutate(pathway=paste(gs_subcat,gs_name,sep="_"))%>% select(pathway,gs_description) %>% distinct()
print(M_pathways_selected,n=50)
print(paste0("You selected a list of N = ",length(M_pathways_selected$pathway)," Pathways and n_genes = ", length(unique(M$gene_symbol))))





#con la siguiente linea de código lo que hace únicamente coger las columnas gs_name y gene_symbol del data set H. Ahora los va a agrupar por gs_name y va a hacer que se junten todos los genes que compartan gs_name en un vector para cada gs_name. Ademas con deframe lo que hace es convertir todo en una lista.
H.Genes.ls_new <- H %>%
  select(gs_name, gene_symbol) %>%
  group_by(gs_name) %>%
  summarise(all.genes = list(unique(gene_symbol))) %>%
  deframe()

M_Genes.ls_new02 <- M %>%
  select(gs_name, gene_symbol) %>%
  group_by(gs_name) %>%
  summarise(all.genes = list(unique(gene_symbol))) %>%
  deframe()

H.Genes.ls_new03 <- Mo %>%
  select(gs_name, gene_symbol) %>%
  group_by(gs_name) %>%
  summarise(all.genes = list(unique(gene_symbol))) %>%
  deframe()

str(normCounts02)
str(normCounts_dds)
normCounts_dds<-as.data.frame(normCounts_dds)


##We calculate the means and the diference between both genotypes (KOvsWT) and create a matrix with the gene names and the delta, which is the diference between KO vs WT.
FC <- as.data.frame(normCounts_dds02) %>%
  rownames_to_column("gene_symbol") %>%
  pivot_longer(-gene_symbol, 
               names_to = "libID", values_to = "expression") %>%
  
  separate(libID, into = c("genotype","sample_number"), sep = "_") %>%
 
  
  pivot_wider(names_from = genotype, values_from = expression) %>%
  
  mutate(delta = Tmem86b-WT) %>%
  
  group_by(gene_symbol) %>%
  summarise(mean.delta = mean(delta, na.rm = TRUE)) %>%
  
  arrange(desc(mean.delta))

##El package necesita que convirtamos la matriz a esta estructura##
FC.vec <- FC$mean.delta
names(FC.vec) <- FC$gene_symbol

str(FC.vec)


min(FC.vec)
max(FC.vec)
scoreType <- "std"

##We run the GSEA analysis with the fgseaSimple function.
gsea.H <- fgseaSimple(pathways = H.Genes.ls_new, 
                      stats = FC.vec,
                      scoreType = scoreType,
                      nperm = 6000)

gsea.H02 <- fgseaSimple(pathways = M_Genes.ls_new02, 
                        stats = FC.vec,
                        scoreType = scoreType,
                        nperm = 1000)

gsea.H03 <- fgseaSimple(pathways = H.Genes.ls_new03, 
                        stats = FC.vec,
                        scoreType = scoreType,
                        nperm = 1000)


class(gsea.H)

gsea.H03 %>%
  filter(pval <= 0.005) %>%
  mutate(pathway = gsub ("HOFFMANN","", pathway),
         pathway = gsub ("HUMMERICH","", pathway),
         pathway = gsub ("MATZUK","", pathway),
         pathway = gsub ("REACTOME","", pathway),
         pathway = gsub ("BOYLAN","", pathway),
         pathway = gsub ("BLANCO_MELO","", pathway),
         pathway = gsub ("NABA","", pathway),
         pathway = gsub ("WP","", pathway),
         pathway = gsub ("PID","", pathway),
         pathway = gsub ("ZHENG","", pathway),
         pathway = gsub ("MILICIC","", pathway),
         pathway = gsub ("WONG","", pathway),
         pathway = gsub ("GAUSSMAN","", pathway),
         pathway = gsub ("LIANG","", pathway),
         pathway = gsub ("LIU","", pathway),
         pathway = gsub ("XU","", pathway),
         pathway = gsub ("WANG","", pathway),
         pathway = gsub ("MOHANKUMAR","", pathway),
         pathway = gsub ("SHETH","", pathway),
         pathway = gsub ("ODONNELL","", pathway),
         pathway = gsub ("FUJIWARA","", pathway),
         pathway = gsub ("PIONTEK","", pathway),
         pathway = gsub ("BURTON","", pathway),
         pathway = gsub ("SCHOEN","", pathway),
         pathway = gsub ("OUELLET","", pathway),
         pathway = gsub ("PEDERSEN","", pathway),
         pathway = gsub ("GOBP","", pathway),
         pathway = gsub ("GOMF","", pathway),
         pathway = gsub ("NEMETH","", pathway),
         pathway = gsub ("BROWNE","", pathway),
         pathway = gsub ("GOCC","", pathway),
         pathway = gsub ("JACKSON","", pathway),
         pathway = gsub("_"," ", pathway),) %>%
  
  ggplot(aes(x = reorder(pathway, NES),
             y = NES)) +
  geom_col() +
  theme_classic() +
  lims(y = c(-2,2)) +
  
  coord_flip() +
  
  labs(y = "Normalizaed enrichment score (NES)",
       x = "Gene set",
       title = "Hallamrk GSEA (Pval<0.01)\n Down in +Mtb <--          --> Up in +Mtb")



########################### dot plot with bars in ggplot2#########################
install.packages("tidyverse")
install.packages("ddply")
library(ggplot2)
library(tidyverse)
library(ggsci)
library(ggpubr)
library(plyr)
library(dplyr)


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ddply")

normCounts_WD <- read.csv("DEXTmem86bWD_NormCounts_dds.csv")
row.names(normCounts_WD)<-normCounts_WD$X
normCounts_WD<-normCounts_WD[-1]



FC02 <- as.data.frame(normCounts_WD) %>%
  rownames_to_column("gene_symbol") %>%
  pivot_longer(-gene_symbol, names_to = "libID", values_to = "expression") %>%
  
  separate(libID, into = c("genotype","sample_number"), sep = "_") %>%
  
  group_by(gene_symbol, genotype) %>%
  
  summarise(mean.exp = mean(expression, na.rm = TRUE)) 



FC02 <- as.data.frame(normCounts_WD) %>%
  rownames_to_column("gene_symbol") %>%
  pivot_longer(-gene_symbol, names_to = "libID", values_to = "expression") %>%
  
  separate(libID, into = c("genotype","sample_number"), sep = "_") %>%
  
  group_by(gene_symbol, genotype) %>%
  
  summarise(mean.exp = mean(expression, na.rm = TRUE), sd_exp = sd(expression))%>%
  
  pivot_wider(names_from = genotype, values_from = c(mean.exp, sd_exp))%>%
  
  mutate(fold_change = mean.exp_Tmem86b/mean.exp_WT)
  
row.names(FC02)<-FC02$gene_symbol
FC02<-FC02%>%select(-gene_symbol)


  
str(FC02)














#group_by(gene_symbol) %>%





#Los anteriores genes de fibrosis que ponía
#"Fn1","Lamb1","Lamb2","Vtn","Thbs1","Thbs2", "Thbs3","Col3a1","Col1a1","Col1a2","Acta1","Acta2","Tgfb1","Tgfb2","Tgfb3","Tgfa","Pdgfa","Pdgfb","Pdgfc","Pdgfd","Map4","Lamb3","Lamb4","Egflam","Col17a1","Ctgf","Vegfa","Vegfb","Vegfc","Vegfd","Mmp2","Mmp9","Igf1","Timp1","Timp2","Timp3","Mmp13","Mmp12"


#fibrosis new all

d<-normCounts_WD[row.names(normCounts_WD)%in%c("Hamp","Orm1","Serpina1a","C3","Vtn","Crp","Fga","Lbp","C4b","C6","Fgb","Fgg","Il6ra","Il1r1","Plg"),]
d<-normCounts_WD[row.names(normCounts_WD)%in%c("Aifm2", "Gpx4", "Slc7a11", "Slc39a7", "Tmem164", "Hmox1", "Ptgs2", "Egfr", "Cav1", "Tlr4", "Cdkn2a", "Rasal3", "Ptgs2", "Cox2", "Tfrc", "Acsl4", "Chac1", "Arntl", "Nfe2l2", "Vdac2", "Vdac3"),]
str(d)
head(d)
d<-FC02[row.names(FC02)%in%c("Ccn2", "Ccnd1", "Cebpb", "Cybb", "Fos", "Fzd8", "Il1a", "Il1r1", "Itga2", "Itgam", "Itgb3", "Mras", "Pdgfb", "Prok1", "Rac2", "Sucnr1", "Tlr4"),]
row.names(d)<-d$gene_symbol
"Adipogenesis" "Bmp7","Cebpb","Cebpd","Egr2","Fgfr3","Fzd8","Lpin1","Lpl","Smad9","Srebf1","Plin5"

"NRF2-mediated Oxidative Stress Response""Abcc1","Cyp2a12","Cyp2a22","Cyp2a6","Cyp2c40","Cyp2c8","Cyp2g1","Cyp2u1","Cyp3a25","Cyp3a5","Cyp4a11","Cyp4a14","Dnajb11","Fkbp5","Fos","Gsta3","Gsta5","Gstm5","Gstm6","Herpud1","Mgst3","Mras"

"Hepatic Fibrosis / Hepatic Stellate Cell Activation"	"Ccn2","Csf1","Cxcl3","Il1a","Il1r1","Il6r","Mmp2","Pdgfb","Prok1","Tlr4"

"Hepatic Fibrosis" "Ccn2","Ccnd1","Cebpb","Cybb","Fos","Fzd8","Il1a","Il1r1","Itga2","Itgam","Itgb3","Mras","Pdgfb","Prok1","Rac2","Sucnr1","Tlr4"

"Hepatic cholestasis" "Abcb1a","Abcc1","Cyp3a5","Il1a","Il1r1","Slc10a2","Srebf1","Tlr4","Tnfsf10"

"p53 Signaling"	"Ccnd1","Cdkn1a","E2f1","Gadd45a","Gadd45g","Serpine2"

"Granulocyte Adhesion and Diapedesis"	"Ccl2","Cxcl10","Cxcl14","Cxcl3","Il1a","Il1r1","Itgam","Mmp2"


"LXR/RXR Activation"	"A1bg","Abcg1","Apoa4","Fga","Il1a","Il1r1","Lpl","Lyz","Pltp","Scd1", "Scd2","Srebf1","Tlr4"

"Apoptosis Signaling"		"Bcl2a1a","Bcl2a1b","Bcl2a1d","Capn11","Capn8","Mras","Naip","Naip1"

"CLEAR Signaling Pathway"		"Atp6v0d2","Ddit4","Fgf21","Fgfr3","Hexb","Mras","Ntrk2","Pdgfb","Tlr1","Tlr2","Tlr4","Tlr6"


"Acute Phase Response Signaling"		"Cebpb","Fga","Fgb","Fgg","Fos","Il1a","Il1r1","Il6ra","Itih3","Mras"

"Coagulation System"	"Fga","Fgb","Fgg","Plau"

"Senescence Pathway"		"Atf3","Capn11","Capn8","Ccnd1","Cdkn1a","Cebpb","E2f1","Gadd45a","Gadd45g","Il1a","Mras","Smad9","Tlr2"

"Ferroptosis" "Aifm2", "Gpx4", "Slc7a11", "Slc39a7", "Tmem164", "Hmox1", "Ptgs2", "Egfr", "Cav1", "Tlr4", "Mapk8", "Cdkn2a", "Rasal3", "Hmgb1", "Ptgs2", "Cox2", "Tfrc", "Acsl4", "Chac1", "Arntl", "Nfe2l2", "Vdac2", "Vdac3", "Kras"


f<-normCounts_all02[row.names(normCounts_all02)%in%c("Fn1","Lama5","Lamb1","Lamb2","Lamc1","Lamc2","Vtn","Thbs1","Thbs3","Col3a1","Col1a1","Col1a2","Acta1","Acta2","Tgfb1","Tgfb2","Tgfb3","Tgfa","Pdgfa","Pdgfb","Pdgfc","Pdgfd","Map4","Lama1","Lama2","Lama3","Lama4","Lamb3","Lamb4","Lamc3","Lamc4","Egflam","Lamb2p1","Col17a1","Ctgf","Vegfa","Vegfb","Vegfc","Vegfd","Mmp2","Mmp9","Igf1","Timp1","Timp2","Timp3","Mmp13","Mmp12"),]
str(d)
head(d)

d<-as.data.frame(d)
f<-as.data.frame(f)

#Pulimos los datos para representar el fold change y que no se vea cada gen por su lado y para ponerlo en buen formato
d <- d %>%
  mutate(mean.WT = (WT_1+WT_2+WT_3+WT_4)/4) %>%
  mutate(mean.TKO = (Tmem86b_1+Tmem86b_2+Tmem86b_3+Tmem86b_4)/4) %>%
  mutate(WT_1bis = WT_1/mean.WT, WT_2bis = WT_2/mean.WT, WT_3bis = WT_3/mean.WT, WT_4bis = WT_4/mean.WT) %>%
  mutate(Tmem86b_1bis = Tmem86b_1/mean.WT, Tmem86b_2bis = Tmem86b_2/mean.WT, Tmem86b_3bis = Tmem86b_3/mean.WT, Tmem86b_4bis = Tmem86b_4/mean.WT) 
  
  mutate(sd.WT = lapply(d[,c("WT_1bis","WT_2bis","WT_3bis","WT_4bis")],sd)) %>%
  mutate(sd.TKO = lapply(d[,c("TKO_1bis","TKO_2bis","TKO_3bis","TKO_4bis")],sd))


d <- d[-c(1,2,3,4,5,6,7,8,9,10)]

d_tris <- d[-c(1,2,3,4,5,6,7,8,11,12,13,14,15,16,17,18)]
d.mean <- d_tris[-c(3,4)]
d.SD <- d_tris[-c(1,2)]

Plot <- as.data.frame(d) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, 
               names_to = "ID", values_to = "expression") %>%
  separate(ID, into = c("genotype","sample_number"), sep = "_") %>%
 
  mutate(genotype_good = factor(genotype, levels = c("WT","Tmem86b")))




Plot02<-Plot %>%
  mutate(genotype_good = factor(genotype, levels = c("WT","Tmem86b")))


#Hago lo mismo pero con un dataframe que tiene todas las muestras a ver que tal
f <- f %>%
  mutate(mean.WT = (WT_1+WT_2+WT_3+WT_4+WT_5)/5) %>%
  mutate(mean.TKO = (TKO_1+TKO_2+TKO_3+TKO_4+TKO_5)/5) %>%
  mutate(WT_1bis = WT_1/mean.WT, WT_2bis = WT_2/mean.WT, WT_3bis = WT_3/mean.WT, WT_4bis = WT_4/mean.WT, WT_5bis = WT_5/mean.WT) %>%
  mutate(TKO_1bis = TKO_1/mean.WT, TKO_2bis = TKO_2/mean.WT, TKO_3bis = TKO_3/mean.WT, TKO_4bis = TKO_4/mean.WT, TKO_5bis = TKO_5/mean.WT)

f_bis <- f[-c(1,2,3,4,5,6,7,8,9,10,11,12)]

Plot_f <- as.data.frame(f_bis) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, 
               names_to = "ID", values_to = "expression") %>%
  separate(ID, into = c("genotype","sample_number"), sep = "_")

Plot02_f<-Plot_f %>%
  mutate(genotype_good = factor(genotype, levels = c("WT","TKO")))




##We represent the data with ggplot
ggplot(data = Plot) + aes(y = expression, x = Gene, fill = genotype_good) + 
  geom_boxplot(position = position_dodge(width = 0.7), outlier.shape = 17, outlier.colour = "red", outlier.size = 0.1, aes(fill = genotype_good)) + 
  stat_summary(geom = 'pointrange', width = 0.9, fatten = 2, color = 'black',
               fun.data = function(x){
                 return(c(
                   y = mean(x),
                   ymin = mean(x),
                   ymax = mean(x)))}) +
  geom_jitter( position = position_jitterdodge(0.5), color = "black") +
  labs(y = "Expression", x = "Gene set",
       title = "Apoptosis")

##este es el que ha hecho Mar y el que sale bien!##
ggplot(data = Plot02, aes(y = expression, x = Gene, fill = genotype_good)) + 
  geom_bar(stat="summary", fun = "mean", position=position_dodge(0.9), color = "black")  +
  geom_jitter( position = position_jitterdodge(0.1), color = "black") +
  theme_bw(base_size=15, base_family = "", base_line_size = NA, base_rect_size = 2) +
  scale_fill_manual(name=NULL, breaks = c("WT","TKO"), values = c("snow2", "#F8766D")) +
  labs(y = "Expression", x = "Gene set",title = "Apoptosis") 

geom_errorbar(data = Plot.Mean, aes(ymin = as.double(expression-SD), ymax = as.double(expression+SD), width = 0.2, position=position_dodge(0.9)))

##este es el que modifico yo para varias cosas y sale bien también##
##He metido las error bar de las STANDARD DEVIATION y modificado como aparecen los puntos##
ggplot(data = Plot02, aes(y = expression, x = Gene, fill = genotype_good)) + 
  geom_bar(stat="summary", fun = "mean", position=position_dodge(0.9), color = "black")  +
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, position = position_dodge(0.9)) +
  theme_bw(base_size=15, base_family = "", base_line_size = NA, base_rect_size = 2) +
  scale_fill_manual(name=NULL, breaks = c("WT","Tmem86b"), values = c("snow2", "#F8766D")) +
  labs(y = "Expression", x = "Gene set",title = "Apoptosis") +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
               geom="errorbar", 
               color = "black", 
               position = position_dodge(0.9))


##esta es una forma de poner el INTERVALO DE CONFIANZA AL 95%##
ggplot(data = Plot_new, aes(y = expression, x = Gene, fill = ID_good)) + 
  geom_bar(stat="summary", fun = "mean", position=position_dodge(0.9), color = "black")  +
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, position = position_dodge(0.9)) +
  theme_bw(base_size=15, base_family = "", base_line_size = NA, base_rect_size = 2) +
  scale_fill_manual(name=NULL, breaks = c("WT","TKO"), values = c("snow2", "#F8766D")) +
  labs(y = "Expression", x = "Gene set",title = "Fibrosis") +
  stat_summary(fun.data="mean_cl_normal",
               geom="errorbar", 
               color = "black", 
               position = position_dodge(0.9),
               width = 0.7)

##esta es una forma de poner el STANDARD ERROR.##
ggplot(data = Plot, aes(y = expression, x = Gene, fill = genotype_good)) + 
  geom_bar(stat="summary", fun = "mean", position=position_dodge(0.9), color = "black", width = 0.8)  +
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.9, position = position_dodge(0.9)) +
  theme_bw(base_size=15, base_family = "", base_line_size = NA, base_rect_size = 2) +
  scale_fill_manual(name=NULL, breaks = c("WT","Tmem86b"), values = c("snow2", "#F8766D")) +
  labs(y = "Expression", x = "Gene set",title = "Acute Pahse Response Downregulated Genes") +
  stat_summary(fun = "mean",
               geom = "errorbar",
               position = position_dodge(0.9),
               fun.max = function(expression) mean(expression) + sd(expression) / sqrt(length(expression)),
               fun.min = function(expression) mean(expression) - sd(expression) / sqrt(length(expression))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18), aspect.ratio = 2/6)





##############################GO analysis###########################

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationDbi")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Mm.eg.db")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tidyverse")

library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)

results = resultados
results<- as.data.frame(results)
results <- subset(results, is.na(padj) == FALSE)
write.csv(results,"resultsTmem86bWD.csv")
new_results<-read.csv("DEXTmem86bWD_results_dds_complete02.csv")

results_sig <- results[results$padj < 0.05,]
t<-results_sig$X
rownames(results_sig)=t
duplicated(t)

dupli<-duplicated(t)
dupli_names<-t[dupli]
print(dupli_names)

genes_to_test_up
genes_to_test_down

genes_to_test_up <- rownames(results_sig[results_sig$log2FoldChange > 0.5,])


GO_results_UP <- enrichGO(gene = genes_to_test_up, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
#Si por lo que fuera a la hora de correr esto en vez de simbolos de genes tuviera los codigos de genes de Ensembl tendría que poner "ENSEMBL", si tubiera entrez ID tendría que poner "ENTREZ"

GO_results_dataframe_UP <- as.data.frame(GO_results_UP)

fit_UP <- plot(barplot(GO_results_UP, showCategory = 50))

fit_UP02 <- plot(barplot(GO_results_UP, x = "pvalue", showCategory = 20))



genes_to_test_down <- rownames(results_sig[results_sig$log2FoldChange < -0.5,])

GO_results_DOWN <- enrichGO(gene = genes_to_test_down, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
#Si por lo que fuera a la hora de correr esto en vez de simbolos de genes tuviera los codigos de genes de Ensembl tendría que poner "ENSEMBL", si tubiera entrez ID tendría que poner "ENTREZ"

GO_results_dataframe_DOWN <- as.data.frame(GO_results_DOWN)

fit_DOWN <- plot(barplot(GO_results_DOWN, showCategory = 50))

fit_DOWN02 <- plot(barplot(GO_results_DOWN, , x = "pvalue", showCategory = 10))



ggplot(GO_results_dataframe_DOWN, aes(x = reorder(Description, pvalue),
                                             y = pvalue)) +
  geom_col() +
  theme_classic() +
  coord_flip() +
  labs(y = "pvalue",
       x = "Gene set",
       title = "Hallamrk GO Downregulated Genes (Pval<0.05)")

library(ggplot2)









write_csv(GO_results_dataframe_UP, "GO_Upregulated_genes Tmem86bWD.csv")
write_csv(GO_results_dataframe_DOWN, "GO_Downregulated_genes Tmem86bWD.csv")













#####################Another GO analysis####################
install.packages("msigdbr")

library(tidyverse)
library(msigdbr)
library(fgsea)
packageVersion("msigdbr")


new_results<-read.csv("resultsTmem86bHFD.csv")

#Esto es para la base de datos que uso
H<-msigdbr(species = "Mus musculus", category = c("H","C1","C2","C3","C4","C5","C6","C7","C8"))
class(H)
H<-msigdbr(species = "Mus musculus", category = "H")
M<-msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP")
M.a<-msigdbr(species = "Mus musculus", category = "C2")

ggplot(results, aes(x=pvalue))+
  geom_histogram(bins=4000) +
  theme_classic()+
  lims(x=c(0,0.05))


table(results$pvalue == 0)


result_sig02 <- results %>%
  filter(pvalue <= 0.05)

result_sig03 <- results %>%
  filter(pvalue <= 0.05) %>%
  filter(log2FoldChange > 0.5)

result_sig04 <- results %>%
  filter(pvalue <= 0.05) %>%
  filter(log2FoldChange < -0.5)


sig_geneNames <- rownames(result_sig02)
sig_geneNamesU <- rownames(result_sig03)
sig_geneNamesD <- rownames(result_sig04)
Ma.genes <- select(M.a,gs_name, gene_symbol)



enrich.H <- enricher(gene = sig_geneNames, TERM2GENE = Ma.genes)
enrich.H_U <- enricher(gene = sig_geneNamesU, TERM2GENE = Ma.genes)
enrich.H_D <- enricher(gene = sig_geneNamesD, TERM2GENE = Ma.genes)

class(enrich.H)

enrich.H@result
##Aquí lo que hago es por un lado, de los resultados, separar BgRatio y GenRatio por / y dar a cada columna un nombre distinto. Luego les da funcion de número, para que piueda operar con ello. Con lo ultimo lo que hace es crear una nueva columna donde hace una operación con los datos que tiene en las dos columnas nuevas creadas.
enrich.H.df <- enrich.H@result %>%
  
  separate(BgRatio, into = c("size.term", "size.category"), sep = "/") %>%
  separate(GeneRatio, into = c("size.overlap.term", "size.overlap.category"), sep = "/") %>%
  
  mutate_at(vars("size.term", "size.category", "size.overlap.term", "size.overlap.category"), as.numeric) %>%
  
  mutate("k.K" = size.overlap.term/size.term)

enrich.H.df_UP <- enrich.H_U@result %>%
  
  separate(BgRatio, into = c("size.term", "size.category"), sep = "/") %>%
  separate(GeneRatio, into = c("size.overlap.term", "size.overlap.category"), sep = "/") %>%
  
  mutate_at(vars("size.term", "size.category", "size.overlap.term", "size.overlap.category"), as.numeric) %>%
  
  mutate("k.K" = size.overlap.term/size.term)

enrich.H.df_DOWN <- enrich.H_D@result %>%
  
  separate(BgRatio, into = c("size.term", "size.category"), sep = "/") %>%
  separate(GeneRatio, into = c("size.overlap.term", "size.overlap.category"), sep = "/") %>%
  
  mutate_at(vars("size.term", "size.category", "size.overlap.term", "size.overlap.category"), as.numeric) %>%
  
  mutate("k.K" = size.overlap.term/size.term)


##Ahora se representan los datos


enrich.H.df %>%
  filter(p.adjust <= 0.00005) %>%
  
  mutate(Description = gsub ("HALLMARK_","", Description),
         Description = gsub("_"," ", Description)) %>%
  
  ggplot(aes(x = reorder(Description, k.K),
             y = k.K)) +
  geom_col() +
  theme_classic() +
  
  coord_flip() +
  
  labs(y = "Significant genes in set / Total genes in set \nk/K",
       x = "Gene set",
       tittle = "Mtb significant genes (Pval<0.05) \ Enriched in Hallmark gene sets (Pval<0.05)")


enrich.H.df_UP %>%
  filter(p.adjust <= 0.005) %>%
  
  mutate(Description = gsub ("HALLMARK_","", Description),
         Description = gsub("_"," ", Description)) %>%
  
  ggplot(aes(x = reorder(Description, k.K),
             y = k.K)) +
  geom_col() +
  theme_classic() +
  
  coord_flip() +
  
  labs(y = "Significant genes in set / Total genes in set \nk/K",
       x = "Gene set",
       tittle = "Mtb significant genes (Pval<0.05) \ Enriched in Hallmark gene sets (Pval<0.05)")


enrich.H.df_DOWN %>%
  filter(p.adjust <= 0.000005) %>%
  
  mutate(Description = gsub ("HALLMARK_","", Description),
         Description = gsub("_"," ", Description)) %>%
  
  ggplot(aes(x = reorder(Description, k.K),
             y = k.K)) +
  geom_col() +
  theme_classic() +
  
  coord_flip() +
  
  labs(y = "Significant genes in set / Total genes in set \nk/K",
       x = "Gene set",
       title = "Mtb significant genes (Pval<0.05) \ Enriched in Hallmark gene sets (Pval<0.05)")











####cosas####


myc<-grep("HALLMARK_",names(Mo$gs_name))





