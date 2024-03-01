#Create two models: FULL and NULL
preserve_function <- paste(colnames(phenodata), collapse=" + ")
full_model <- reformulate(preserve_function)
null_model <- reformulate('1')

cat("Print models\n")
print(full_model)
print(null_model)

phenodata <- phenodata[order(rownames(phenodata)),,drop=F]
assayData <- assayData[,order(colnames(assayData))]


cat("Creating full model: variables to adjust and preserve\n")
mod = model.matrix(full_model, data=phenodata)

cat("Creating null model: only adjustment\n")
mod0 <- model.matrix(null_model , data=phenodata)


assayData <- assayData[, which(colnames(assayData) %in% rownames(mod))]
mod <- mod[which(rownames(mod) %in% colnames(assayData)),]
mod0 <- mod0[which(rownames(mod0) %in% colnames(assayData)),]


all(colnames(assayData)==rownames(mod))

# remove infinite values
assayData <- assayData[is.finite(rowSums(assayData)),]

cat("Calculate number of SVs and printing\n")
n.sv = num.sv(assayData,mod,method="be")
print(n.sv) # BasoAmyg-25, dACC-21, DLPFC-22, MedialAmyg-27

cat("Applying SVA to estimate SVs\n")
svobj = sva(assayData,mod,mod0,n.sv=n.sv)

cat("Saving Outputs\n")

cat("Surrogate Variables\n")
#Save SVs in full model
colnames(svobj$sv) <- paste0("SV",1:n.sv)
modSv = cbind(mod,svobj$sv)
modSv_save <- as.data.frame(modSv) 
modSv_save <- rownames_to_column(modSv_save, 'sampleid')
write.table(modSv_save, file=paste0("/sc/arion/projects/psychgen/collab/postmortem_brain/expression/sva_preserve_trauma/",tissue,"_SVA_preserve_total_trauma.SVs"), quote=F, sep="\t", col.names=T, row.names=F)

cat("Saving Residuals\n")
#Generate residuals
#Save SVs in null model
mod0Sv = cbind(mod0,svobj$sv)
fit=lmFit(assayData, mod0Sv)
residuals = residuals.MArrayLM(fit, assayData)
residualsdf <- as.data.frame(residuals)
residual_tosave <- residualsdf %>% rownames_to_column('gene_id')
write.table(residual_tosave, file=paste0("/sc/arion/projects/psychgen/collab/postmortem_brain/expression/sva_preserve_trauma/",tissue,"_SVA_preserve_total_trauma.residuals"), quote=F, sep="\t", col.names=T, row.names=F)

cat("Saving posterior probabilities\n")
pg <- cbind(rownames(assayData),svobj$pprob.gam)
colnames(pg) <- c('Name', 'pprob.gam')
write.table(pg, file=paste0("/sc/arion/projects/psychgen/collab/postmortem_brain/expression/sva_preserve_trauma/",tissue,"_SVA_preserve_total_trauma.pprobgam"), quote=F, sep="\t", col.names=T, row.names=F)

pb <- cbind(rownames(assayData),svobj$pprob.b)
colnames(pb) <- c('Name', 'pprob.b')
write.table(pb, file=paste0("/sc/arion/projects/psychgen/collab/postmortem_brain/expression/sva_preserve_trauma/",tissue,"_SVA_preserve_total_trauma.pprobb"), quote=F, sep="\t", col.names=T, row.names=F)

cat("DONE!\n")
