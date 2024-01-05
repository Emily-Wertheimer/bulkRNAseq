## modified from jasmine

## merge fastq files from same read/sample run in different lanes (bash) ##
# load R in bash
ml R/4.3.0-foss-2020b
R

# install packages/libraries, set wd
library(Rsubread)
setwd("/gpfs/gibbs/pi/huckins/ekw28/ghre_lep_sema_RNAseq/")

# make objects with filenames
read1.fastq.files <- list.files(path = "./sample_dir/merged_all_lanes", pattern = "R1_all_lanes.fastq.gz$", full.names = TRUE)
read2.fastq.files <- list.files(path = "./sample_dir/merged_all_lanes", pattern = "R2_all_lanes.fastq.gz$", full.names = TRUE)

#paths
output_dir <- "/gpfs/gibbs/pi/huckins/ekw28/ghre_lep_sema_RNAseq/bam_aligned_files"
index_path <- "/gpfs/gibbs/pi/huckins/ekw28/ghre_lep_sema_RNAseq/hg38_index_files/hg38_index"

# Assuming read1 and read2 files are in the same order and correspond to each other
for(i in 1:length(read1.fastq.files)) {
    sample_name <- gsub(".*\\/(.*)_R1_all_lanes.fastq.gz", "\\1", read1.fastq.files[i])
    output_file_path <- paste0(output_dir, sample_name, ".bam")

    align(index = index_path, read1.fastq.files[i], read2.fastq.files[i],
          type = "rna", input_format = "gzFASTQ", output_format = "BAM",
          output_file = output_file_path,  phredOffset = 33, nsubreads = 10, 
      TH1 = 3, TH2 = 1, maxMismatches = 3, unique = FALSE, nBestLocations = 1, 
      indels = 5, complexIndels = FALSE, nTrim5 = 0, nTrim3 = 0, 
      minFragLength = 50, maxFragLength = 600, PE_orientation = "fr", 
      nthreads = 1, readGroupID = NULL, readGroup = NULL, keepReadOrder = FALSE, 
      sortReadsByCoordinates = FALSE, color2base = FALSE, DP_GapOpenPenalty = -1, 
      DP_GapExtPenalty = 0, DP_MismatchPenalty = 0, DP_MatchScore = 2, 
      detectSV = FALSE, useAnnotation = FALSE, annot.inbuilt = "hg38", 
      annot.ext = NULL, isGTF = FALSE, GTF.featureType = "exon", 
      GTF.attrType = "gene_id", chrAliases = NULL) 
          
          
}

#####
# align
align(index = "./hg38_index_files/hg38_index", read1.fastq.files, read2.fastq.files,
      type = "rna", input_format = "gzFASTQ", output_format = "BAM",
      output_file = paste(output_dir,read1.fastq.files, "subread", "BAM", sep = "."),
      phredOffset = 33, nsubreads = 10, 
      TH1 = 3, TH2 = 1, maxMismatches = 3, unique = FALSE, nBestLocations = 1, 
      indels = 5, complexIndels = FALSE, nTrim5 = 0, nTrim3 = 0, 
      minFragLength = 50, maxFragLength = 600, PE_orientation = "fr", 
      nthreads = 1, readGroupID = NULL, readGroup = NULL, keepReadOrder = FALSE, 
      sortReadsByCoordinates = FALSE, color2base = FALSE, DP_GapOpenPenalty = -1, 
      DP_GapExtPenalty = 0, DP_MismatchPenalty = 0, DP_MatchScore = 2, 
      detectSV = FALSE, useAnnotation = FALSE, annot.inbuilt = "hg38", 
      annot.ext = NULL, isGTF = FALSE, GTF.featureType = "exon", 
      GTF.attrType = "gene_id", chrAliases = NULL) 
