## align fastQ files to reference genome in parallel using mcapply

# install packages/libraries, set wd
library(Rsubread)
library(parallel)
setwd("/gpfs/gibbs/pi/huckins/ekw28/ghre_lep_sema_RNAseq/")

output_dir <- '/gpfs/gibbs/pi/huckins/ekw28/ghre_lep_sema_RNAseq/bam_aligned_files'
fastq_dir <- '/gpfs/gibbs/pi/huckins/ekw28/ghre_lep_sema_RNAseq/sample_dir/merged_all_lanes_old/553gaba_NGN2_553astro'
index_path <- '/gpfs/gibbs/pi/huckins/ekw28/ghre_lep_sema_RNAseq/hg38_index_files/hg38_index'

fastq_files <- list.files(fastq_dir, pattern = "_R1.*_all_lanes\\.fastq\\.gz$", full.names = TRUE)
sample_names <- gsub("_R1.*_all_lanes\\.fastq\\.gz$", "", basename(fastq_files)) # Use fastq_files to define sample_names

# Create a function that takes a sample name and performs the alignment using align().
align_and_save <- function(sample_names) {
  readfile1 <- file.path(fastq_dir, paste0(sample_names, "_R1_all_lanes.fastq.gz"))
  readfile2 <- file.path(fastq_dir, paste0(sample_names, "_R2_all_lanes.fastq.gz"))
  output_file <- file.path(output_dir, paste0(sample_names, ".bam"))

  # Run the alignment
  align(index=index_path,
        readfile1=readfile1,
        readfile2=readfile2,
        output_file=output_file,
        type='rna',
        minFragLength=50,
        maxFragLength=600
       )
       
    return(paste0("Alignment completed for ", sample_name))
}


# Detect the number of available cores
num_cores <- detectCores() - 1

# Run the alignment in parallel & save output
results <- mclapply(sample_names, align_and_save, mc.cores = num_cores)
