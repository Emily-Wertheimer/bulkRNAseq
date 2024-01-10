# install seqtk environment 
git clone https://github.com/lh3/seqtk.git;
cd seqtk; make

####### shell script #########
#!/bin/bash
module load seqtk

# Directory containing FASTQ files
FASTQ_DIR="/gpfs/gibbs/pi/huckins/ekw28/ghre_lep_sema_RNAseq/seqtk/sample_dir/merged_all_lanes"

# Directory to output paired reads
PAIRED_DIR="/gpfs/gibbs/pi/huckins/ekw28/ghre_lep_sema_RNAseq/seqtk/sample_dir/merged_all_lanes_paired"

# Directory to output unpaired reads
SINGLE_DIR="/gpfs/gibbs/pi/huckins/ekw28/ghre_lep_sema_RNAseq/seqtk/sample_dir/merged_all_lanes_UNpaired"

# Loop through R1 files
for r1_file in $FASTQ_DIR/*_R1_*.fastq.gz; do
    # Corresponding R2 file
    r2_file="${r1_file/_R1_/_R2_}"

    # Output file names
    r1_paired="${PAIRED_DIR}/$(basename $r1_file)"
    r2_paired="${PAIRED_DIR}/$(basename $r2_file)"
    singles="${SINGLE_DIR}/$(basename ${r1_file/_R1_/_singles_})"

    # Run seqtk
    seqtk pairfq -t m1 -f $r1_file -r $r2_file -1 $r1_paired -2 $r2_paired -s $singles
done
#############

## make executable


