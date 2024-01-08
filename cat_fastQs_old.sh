### merge fastQs for reads of same sample run across multiple lanes and save as new file (paired end sequencing)

##save as .sh##
#!/bin/bash

# Directory where your FASTQ files are located
FASTQ_DIR="/gpfs/gibbs/pi/huckins/ekw28/ghre_lep_sema_RNAseq/sample_dir"

# Navigate to the FASTQ directory
cd $FASTQ_DIR

# Loop through each sample
for sample in `ls *R1_001.fastq.gz | sed 's/_L00[0-9]_R1_001.fastq.gz//' | uniq`
do
    # Define a unique name for the concatenated output file for R1
    output_R1="${sample}_R1_all_lanes.fastq.gz"

    # Define a unique name for the concatenated output file for R2
    output_R2="${sample}_R2_all_lanes.fastq.gz"

    # Concatenate all R1 files for the sample
    cat ${sample}_L00*_R1_001.fastq.gz > $output_R1

    # Concatenate all R2 files for the sample (if paired-end)
    cat ${sample}_L00*_R2_001.fastq.gz > $output_R2
done
##

## in command line
# make executable
chmod +x cat_fastQs.sh

# run script
./cat_fastQs.sh
