## modified from jasmine

library(Rsubread)
setwd("/gpfs/gibbs/pi/huckins/ekw28/ghre_lep_sema_RNAseq")

# make objects with filenames
read1.fastq.files <- list.files(path = "./sample_dir", pattern = "L001_R1_001.fastq.gz$", full.names = TRUE)
read2.fastq.files <- list.files(path = "./sample_dir", pattern = "L001_R2_001.fastq.gz$", full.names = TRUE)
index_path <- "/gpfs/gibbs/pi/huckins/ekw28/ghre_lep_sema_RNAseq/hg38_index_files/hg38_index"

# check if order for R1 and R2 are the same
read1.fastq.files
read2.fastq.files

# align
align(index = "./hg38_index_files/hg38_index", read1.fastq.files, read2.fastq.files,
      type = "rna", input_format = "gzFASTQ", output_format = "BAM",
      output_file = paste(read1.fastq.files, "subread", "BAM", sep = "."),
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
