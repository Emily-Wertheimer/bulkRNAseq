library(ggplot2)
library(dplyr)
install.packages("ggrepel")
library(ggrepel)

# --------------------------------------------------------------------------- IMPORT DEG tables
setwd('/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/3_DEG')

iAstro_ghre_v_veh <- read.delim("/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/3_DEG/DEG_response_iAstro_ghre_v_veh.txt")
iAstro_lep_v_veh <- read.delim('/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/3_DEG/DEG_response_iAstro_lep_v_veh.txt')
iAstro_sema_v_veh <- read.delim('/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/3_DEG/DEG_response_iAstro_sema_v_veh.txt')

iGlut_ghre_v_veh <- read.delim("/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/3_DEG/DEG_response_iGlut_ghre_v_veh.txt")
iGlut_lep_v_veh <- read.delim('/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/3_DEG/DEG_response_iGlut_lep_v_veh.txt')
iGlut_sema_v_veh <- read.delim('/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/3_DEG/DEG_response_iGlut_sema_v_veh.txt')

iGABA_ghre_v_veh <- read.delim("/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/3_DEG/DEG_response_iGABA_ghre_v_veh.txt")
iGABA_lep_v_veh <- read.delim('/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/3_DEG/DEG_response_iGABA_lep_v_veh.txt')
iGABA_sema_v_veh <- read.delim('/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/3_DEG/DEG_response_iGABA_sema_v_veh.txt')

# --------------------------------------------------------------------------- EXTRACT DEG
# TOP 5 DEG/CONDITION

iAstro_ghre_v_veh_sig <- iAstro_ghre_v_veh[order(iAstro_ghre_v_veh$adj.P.Val),] %>% head(5) 
iAstro_lep_v_veh_sig <- iAstro_lep_v_veh[order(iAstro_lep_v_veh$adj.P.Val),] %>% head(5) 
iAstro_sema_v_veh_sig <- iAstro_sema_v_veh[order(iAstro_sema_v_veh$adj.P.Val),] %>% head(5) 

iGlut_ghre_v_veh_sig <- iGlut_ghre_v_veh[order(iGlut_ghre_v_veh$adj.P.Val),] %>% head(5) 
iGlut_lep_v_veh_sig <- iGlut_lep_v_veh[order(iGlut_lep_v_veh$adj.P.Val),] %>% head(5) 
iGlut_sema_v_veh_sig <- iGlut_sema_v_veh[order(iGlut_sema_v_veh$adj.P.Val),] %>% head(5) 

iGABA_ghre_v_veh_sig <- iGABA_ghre_v_veh[order(iGABA_ghre_v_veh$adj.P.Val),] %>% head(5) 
iGABA_lep_v_veh_sig <- iGABA_lep_v_veh[order(iGABA_lep_v_veh$adj.P.Val),] %>% head(5) 
iGABA_sema_v_veh_sig <- iGABA_sema_v_veh[order(iGABA_sema_v_veh$adj.P.Val),] %>% head(5) 

# all adj pval 
iAstro_ghre_v_veh_all <- iAstro_ghre_v_veh[order(iAstro_ghre_v_veh$adj.P.Val),] 
iAstro_lep_v_veh_all <- iAstro_lep_v_veh[order(iAstro_lep_v_veh$adj.P.Val),] 
iAstro_sema_v_veh_all <- iAstro_sema_v_veh[order(iAstro_sema_v_veh$adj.P.Val),] 

iGlut_ghre_v_veh_all <- iGlut_ghre_v_veh[order(iGlut_ghre_v_veh$adj.P.Val),] 
iGlut_lep_v_veh_all <- iGlut_lep_v_veh[order(iGlut_lep_v_veh$adj.P.Val),] 
iGlut_sema_v_veh_all <- iGlut_sema_v_veh[order(iGlut_sema_v_veh$adj.P.Val),] 

iGABA_ghre_v_veh_all <- iGABA_ghre_v_veh[order(iGABA_ghre_v_veh$adj.P.Val),]
iGABA_lep_v_veh_all <- iGABA_lep_v_veh[order(iGABA_lep_v_veh$adj.P.Val),] 
iGABA_sema_v_veh_all <- iGABA_sema_v_veh[order(iGABA_sema_v_veh$adj.P.Val),] 

# all reg. pval (make sure to color dots by bonf sig)
iAstro_ghre_v_veh_all_nonAdj <- iAstro_ghre_v_veh[order(iAstro_ghre_v_veh$P.Value),] 
iAstro_lep_v_veh_all_nonAdj <- iAstro_lep_v_veh[order(iAstro_lep_v_veh$P.Value),] 
iAstro_sema_v_veh_all <- iAstro_sema_v_veh[order(iAstro_sema_v_veh$P.Value),] 

iGlut_ghre_v_veh_all_nonAdj <- iGlut_ghre_v_veh[order(iGlut_ghre_v_veh$P.Value),] 
iGlut_lep_v_veh_all_nonAdj <- iGlut_lep_v_veh[order(iGlut_lep_v_veh$P.Value),] 
iGlut_sema_v_veh_all_nonAdj <- iGlut_sema_v_veh[order(iGlut_sema_v_veh$P.Value),] 

iGABA_ghre_v_veh_all_nonAdj <- iGABA_ghre_v_veh[order(iGABA_ghre_v_veh$P.Value),]
iGABA_lep_v_veh_all_nonAdj <- iGABA_lep_v_veh[order(iGABA_lep_v_veh$P.Value),] 
iGABA_sema_v_veh_all_nonAdj <- iGABA_sema_v_veh[order(iGABA_sema_v_veh$P.Value),] 

# --------------------------------------------------------------------------- PLOT F(X)
makeVolcanoPlotAndSave <- function(df, title = "Volcano Plot", filename = "volcano_plot.png") {
  p <- df %>% 
    ggplot(aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(alpha = 0.5) +
    geom_text_repel(aes(label = ifelse(adj.P.Val < 0.05, as.character(Gene), '')),
                    box.padding = 0.35, 
                    point.padding = 0.5,
                    segment.color = 'grey50') +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
    theme_minimal()
  
  # Save the plot
  setwd('/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/3_DEG/figures')
  ggsave(filename, plot = p, width = 10, height = 8, dpi = 300)
}


# --------------------------------------------------------------------------- SAVE
setwd('/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/3_DEG/figures')

# write out
makeVolcanoPlotAndSave(iAstro_ghre_v_veh_all_nonAdj, "iAstro: Ghrelin (100nm) v. Vehicle", "iAstro_ghre_v_veh_all.png")
makeVolcanoPlotAndSave(iAstro_lep_v_veh_all_nonAdj, "iAstro: Leptin (100nm) v. Vehicle", "iAstro_lep_v_veh_all.png")
makeVolcanoPlotAndSave(iAstro_sema_v_veh_all_nonAdj, "iAstro: Sema (100nm) v. Vehicle", "iAstro_sema_v_veh_all.png")

makeVolcanoPlotAndSave(iGlut_ghre_v_veh_all, "iGlut: Ghrelin (100nm) v. Vehicle", "iGlut_ghre_v_veh_all.png")
makeVolcanoPlotAndSave(iGlut_lep_v_veh_all, "iGlut: Leptin (100nm) v. Vehicle", "iGlut_lep_v_veh_all.png")
makeVolcanoPlotAndSave(iGlut_sema_v_veh_all, "iGlut: Sema (100nm) v. Vehicle", "iGlut_sema_v_veh_all.png")

makeVolcanoPlotAndSave(iGABA_ghre_v_veh_all, "iGABA: Ghrelin (100nm) v. Vehicle", "iGABA_ghre_v_veh_all.png")
makeVolcanoPlotAndSave(iGABA_lep_v_veh_all, "iGABA: Leptin (100nm) v. Vehicle", "iGABA_lep_v_veh_all.png")
makeVolcanoPlotAndSave(iGABA_sema_v_veh_all, "iGABA: Sema (100nm) v. Vehicle", "iGABA_sema_v_veh_all.png")

