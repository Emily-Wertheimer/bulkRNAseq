#!/bin/bash
#SBATCH --output featureCounts-%A_%1a-%N.out
#SBATCH --job-name featureCounts
#SBATCH --mem-per-cpu 10G -c 8 -t 12:00:00 --mail-type ALL
​
## determine RNA strandedness tool
# http://rseqc.sourceforge.net/
​
​
## featureCounts for-loops
# go to directory containing bam files
cd /vast/palmer/scratch/huckins/cl2375/RNAseq_LentiInflam_Abeta/alignment
# featureCounts is part of Subread package (https://bioconductor.org/packages/release/bioc/vignettes/Rsubread/inst/doc/SubreadUsersGuide.pdf)
module load Subread
for sample in Ctrl AB D4 D7; do
    featureCounts -p -g gene_id -s 2 -t exon \ # -p for pair-end bam files, -g for features to count, -s 2 meant reverse-stranded
    -a /gpfs/gibbs/project/huckins/cl2375/hg38/gencode.v40.primary_assembly.annotation.gtf \ # annotation GTF file
    -o /vast/palmer/scratch/huckins/cl2375/RNAseq_LentiInflam_Abeta/readcounts/${sample}_featureCounts.txt \ # out put txt name
    ${sample}_*_Aligned.sortedByCoord.out.bam # input bam files
done
​
# make a meta file for filtercounts
featureCounts -p -B -T 4 -s 2 \
-a /gpfs/gibbs/project/huckins/cl2375/hg38/gencode.v40.primary_assembly.annotation.gtf \
-o /vast/palmer/scratch/huckins/cl2375/RNAseq_LentiInflam_Abeta/readcounts/all_featureCounts.txt \
*_Aligned.sortedByCoord.out.bam
