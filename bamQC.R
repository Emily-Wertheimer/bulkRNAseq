## adapted from Jasmine & Kayla

#####in shell ######
## find unique .bam files & move to new dir
cd /gpfs/gibbs/pi/huckins/ekw28/ghre_lep_sema_RNAseq/bam_aligned_files
mkdir -p /gpfs/gibbs/pi/huckins/ekw28/ghre_lep_sema_RNAseq/bam_aligned_files/dotBAM_only
find cd /gpfs/gibbs/pi/huckins/ekw28/ghre_lep_sema_RNAseq/bam_aligned_files -type f -regex ".*\.bam$" -exec cp {} /gpfs/gibbs/pi/huckins/ekw28/ghre_lep_sema_RNAseq/bam_aligned_files/dotBAM_only \;





--
ls
find /gpfs/gibbs/pi/huckins/ekw28/ghre_lep_sema_RNAseq/bam_aligned_files -type f -name "*.bam*" | sort | uniq

# make new wd w/ only unique .bam files
mkdir /path/to/new_directory
find /path/to/search -type f -name "*pattern*" -exec cp {} /path/to/new_directory \;



##### in R ######
setwd("/gpfs/gibbs/pi/huckins/ekw28/ghre_lep_sema_RNAseq/bam_aligned_files")





# Check quality scores of reads based on the BAM files
library(Rsubread)
bam.files <- list.files(path = "./data/combined_treatment/bam/merged",
                       pattern = "_merged.BAM$", full.names = TRUE)

for (x in bam.files) {
	qs <- qualityScores(filename = x,
			nreads = 100,
			input_format="BAM")
	write.csv(qs, paste("./results/combined_treatment/quality_scores/quality_", basename(x), ".csv"))
}
