#!/usr/bin/env Rscript
#run variance partition by sample
library(tidyverse)
library(data.table)
library(variancePartition)
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
library(reshape)

setwd("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/projects/trauma_DEGs/coexpression/QC")


args = commandArgs(trailingOnly=TRUE)
tissue=args[1]
#tissue="DLPFC"


# 1. Variance Partitioning for each region
# this tells you whether confounders are a significant contributor to variance 

#get phenotype/demographic data
meta <- data.frame(read.delim("Girgenti_brains_demographic_data_clean_with_trauma_cols_by_sample.txt"))

#filter by region
meta_region <- meta %>% filter(Region==tissue)
rownames(meta_region) <- meta_region$BrNum

#collect RNums
samples_region <- meta_region %>% pull(RNum)
samples_BrNums <- meta_region %>% pull(BrNum)

#load expression matrix
voom_matrix <- data.frame(fread("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/projects/trauma_DEGs/coexpression/data/postmortem_brain_voom_matrix_chr1_22.txt"))
voom_matrix <- voom_matrix %>% column_to_rownames("gene_id")

voom_matrix_region <- voom_matrix[,c(samples_region)] #grab just the samples for reigon of interest
colnames(voom_matrix_region) <- c(samples_BrNums) #rename by brnums

#make sure the order of samples is the same
colnames(voom_matrix_region)==rownames(meta_region)

voom_matrix_region <- voom_matrix_region[,which(colnames(voom_matrix_region) %in% rownames(meta_region))] # make sure same people
meta_region <- subset(meta_region, rownames(meta_region) %in% colnames(voom_matrix_region))


meta_region$Sex <- as.factor(meta_region$Sex)
meta_region$Race <- as.factor(meta_region$Race)
meta_region$plate <- as.factor(meta_region$plate)
meta_region$Group <- as.factor(meta_region$Group)
meta_region$Smoking <- as.factor(meta_region$Smoking)
meta_region$GERD <- as.factor(meta_region$GERD)
meta_region$Atherosclerotic_CVD <- as.factor(meta_region$Atherosclerotic_CVD)
meta_region$Hypertension <- as.factor(meta_region$Hypertension)
meta_region$Migraine <- as.factor(meta_region$Migraine)
meta_region$Seizures <- as.factor(meta_region$Seizures)
meta_region$Chronic_Pain <- as.factor(meta_region$Chronic_Pain)
meta_region$Lifetime.Antipsych <- as.factor(meta_region$Lifetime.Antipsych)
meta_region$Lifetime.Anticonvulsant <- as.factor(meta_region$Lifetime.Anticonvulsant)
meta_region$Lifetime.Antidepress <- as.factor(meta_region$Lifetime.Antidepress)
meta_region$Lifetime.Lithium <- as.factor(meta_region$Lifetime.Lithium)
meta_region$Past.SI <- as.factor(meta_region$Past.SI)
meta_region$Past.Self.Mutilation <- as.factor(meta_region$Past.Self.Mutilation)
meta_region$Psychosis <- as.factor(meta_region$Psychosis)
meta_region$Hallucinations <- as.factor(meta_region$Hallucinations)
meta_region$ECT <- as.factor(meta_region$ECT)
meta_region$Delusions <- as.factor(meta_region$Delusions)
meta_region$trauma_any <- as.factor(meta_region$trauma_any)

attach(meta_region)

#set up formula

#form <- ~  RIN + PMI + (1|plate) + AgeDeath + (1|Sex) + (1|Race) + (1|Smoking) + Education_years + BMI..calculated. + (1|Smoking) + (1|GERD) + (1|Atherosclerotic_CVD) + (1|Hypertension) + (1|Migraine) + (1|Seizures) + (1|Chronic_Pain) + (1|Lifetime.Antipsych) + (1|Lifetime.Anticonvulsant) + (1|Lifetime.Antidepress) + (1|Lifetime.Lithium) + Past.Suicide.Attempts + (1|Past.SI) + (1|Past.Self.Mutilation) + (1|Psychosis) + (1|Hallucinations) + (1|ECT) + (1|Delusions) + (1|Group)

## Smoking, BMI..calculated., Seizures, Lifetime.Antipsych, Lifetime.Anticonvulsant, Lifetime.Antidepress, Lifetime.Lithium, Past.Suicide.Attempts, Past.SI, Past.Self.Mutilation, ECT, Delusions all have missing values remove those

form <- ~  RIN + PMI + (1|plate) + AgeDeath + (1|Sex) + (1|Race) + (1|Smoking) + Education_years + (1|GERD) + (1|Atherosclerotic_CVD) + (1|Hypertension) + (1|Migraine) + (1|Chronic_Pain) + (1|Psychosis) + (1|Hallucinations) + (1|Group)

varPart <- fitExtractVarPartModel(voom_matrix_region, form, meta_region )
vp <- sortCols(varPart)

setwd("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/projects/trauma_DEGs/coexpression/QC/variance_partition")

write.table(vp, sprintf("variancePartition_output_%s.txt",tissue), sep="\t", quote=F, row.names=F)


#vp = data.frame(fread(sprintf("variancePartition_output_%s.txt",tissue)))

#plotting to vizualize
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

png(sprintf("variancePartition_%s.png",tissue), res=300, units="in", width=10, height=7)
p
dev.off()


setwd("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/projects/trauma_DEGs/coexpression/QC/residuals")

## choose covariates 
#alanna's criteria: continuous technical factors that explained ≥ 1% variation in ≥ 10% of genes
cat("Choosing Covariates\n")
test=colnames(vp)[1:length(colnames(vp))-1]
covs_to_adjust=c()
for(col_i in test) {
	print(col_i)
	if(length(which(vp[,col_i] >= 0.01))/(nrow(vp)) > 0.1) covs_to_adjust=c(covs_to_adjust,col_i)
}

write.table(covs_to_adjust, file=sprintf("Girgenti_expression_matrix_autosomes_%s_allVarPartCovs.txt",tissue), quote=F, sep="\n", col.names=T, row.names=F)


##### were going to make several residual matrices #### 

#a helper function

get_residuals <- function(my_covs_to_adjust) {
	adjust_function <- paste(my_covs_to_adjust, collapse=" + ")
	adjust_model <- reformulate(adjust_function)
	mod = model.matrix(adjust_model, data=meta_region)

	#double check samples line up
	voom_matrix_region <- voom_matrix_region[, which(colnames(voom_matrix_region) %in% rownames(mod))]
	mod <- mod[which(rownames(mod) %in% colnames(voom_matrix_region)),]
	all(colnames(voom_matrix_region)==rownames(mod))

	#Generate residuals
	fit=lmFit(voom_matrix_region, mod)
	residuals = residuals.MArrayLM(fit, voom_matrix_region)
	residualsdf <- as.data.frame(residuals)
	residual_tosave <- residualsdf %>% rownames_to_column('gene_id')
	return(residual_tosave)
}


#####  mod1 = residualize all above covariates
#####  mod2 = residualize all above covariates + trauma any
#####  mod3 = residualize all above covariates + trauma count
#####  mod4 = residualize all above covariates that are technical/demographic covs
#####  mod5 = residualize all above covariates that are technical/demographic covs + trauma any
#####  mod6 = residualize all above covariates that are technical/demographic covs + trauma count


#####  mod1 = residualize all above covariates
#get residuals
residual_tosave1 <- get_residuals(covs_to_adjust)
write.table(residual_tosave1, file=sprintf("Girgenti_expression_matrix_autosomes_%s_correct_allVarPartCovs.residuals",tissue), quote=F, sep="\t", col.names=T, row.names=F)


#####  mod2 = residualize all above covariates + trauma any
covs_to_adjust_tr <- c(covs_to_adjust,"trauma_any")
#get residuals
residual_tosave2 <- get_residuals(covs_to_adjust_tr)
#save
write.table(residual_tosave2, file=sprintf("Girgenti_expression_matrix_autosomes_%s_correct_allVarPartCovs_plus_trauma_any.residuals",tissue), quote=F, sep="\t", col.names=T, row.names=F)

#####  mod3 = residualize all above covariates + trauma any
covs_to_adjust_tr_count <- c(covs_to_adjust,"trauma_count")
#get residuals
residual_tosave3 <- get_residuals(covs_to_adjust_tr_count)
#save
write.table(residual_tosave3, file=sprintf("Girgenti_expression_matrix_autosomes_%s_correct_allVarPartCovs_plus_trauma_count.residuals",tissue), quote=F, sep="\t", col.names=T, row.names=F)


#####  mod4 = residualize all above covariates that are technical adn demographic covs
technical_covs=c("plate", "RIN", "PMI", "Sex", "AgeDeath", "Race")
covs_to_adjust_tech <- intersect(covs_to_adjust,technical_covs)
#save
write.table(covs_to_adjust_tech, file=sprintf("Girgenti_expression_matrix_autosomes_%s_TechDemoVarPartCovs.txt",tissue), quote=F, sep="\n", col.names=T, row.names=F)
#get residuals
residual_tosave4 <- get_residuals(covs_to_adjust_tech)

#save residual matrix 1 
write.table(residual_tosave4, file=sprintf("Girgenti_expression_matrix_autosomes_%s_correct_TechDemoVarPartCovars.residuals",tissue), quote=F, sep="\t", col.names=T, row.names=F)



#####  mod5 = residualize all above covariates that are technical adn demographic covs + trauma any
covs_to_adjust_tech_tr <- c(covs_to_adjust_tech,"trauma_any")
#get residuals
residual_tosave5 <- get_residuals(covs_to_adjust_tech_tr)
#save
write.table(residual_tosave5, file=sprintf("Girgenti_expression_matrix_autosomes_%s_correct_TechDemoVarPartCovars_plus_trauma_any.residuals",tissue), quote=F, sep="\t", col.names=T, row.names=F)


#####  mod6 = residualize all above covariates that are technical adn demographic covs + trauma any
covs_to_adjust_tech_tr_count <- c(covs_to_adjust_tech,"trauma_count")
#get residuals
residual_tosave6 <- get_residuals(covs_to_adjust_tech_tr_count)
#save
write.table(residual_tosave6, file=sprintf("Girgenti_expression_matrix_autosomes_%s_correct_TechDemoVarPartCovars_plus_trauma_count.residuals",tissue), quote=F, sep="\t", col.names=T, row.names=F)



#softlink file to files_to_use
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/projects/trauma_DEGs/coexpression/files_to_use


#### to run this script ###
## cd /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/projects/trauma_DEGs/coexpression/QC/variance_partition

# for tissue in DLPFC dACC MedialAmyg BasoAmyg
# do
# 	bsub -P acc_psychgen -q premium -W 12:00 -n 8 -R rusage[mem=5000] -R span[hosts=1] -o var_part_${tissue}.o -e var_part_${tissue}.e /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/projects/trauma_DEGs/coexpression/scripts/QC_variance_partition.R $tissue
# done
