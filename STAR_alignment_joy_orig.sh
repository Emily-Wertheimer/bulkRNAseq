#!/bin/bash
#SBATCH --partition day
#SBATCH --output STAR-%A_%1a-%N.out
#SBATCH --job-name STAR_alignment
#SBATCH --mem-per-cpu 10G -c 10 -t 1-00:00:00 --mail-type ALL
​
module load STAR
​
for sample in /vast/palmer/scratch/huckins/cl2375/RNAseq_LentiInflam_Abeta/rawdata/*; do
    name=$(echo ${sample##*/})
    cd $sample/Unaligned/
    STAR --runThreadN 4 \
    --runMode alignReads \
    --genomeDir /gpfs/gibbs/project/huckins/cl2375/hg38/STAR_index/ \
    --readFilesIn *_R1_001.fastq.gz *_R2_001.fastq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix /vast/palmer/scratch/huckins/cl2375/RNAseq_LentiInflam_Abeta/alignment/${name}_ \
    --outSAMtype BAM SortedByCoordinate \
    --outFilterMismatchNmax 2
done
