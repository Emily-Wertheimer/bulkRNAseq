### code to cat fastQs from different lanes (from Carina Seah, PhD)

# navigate to sample directory

# find fastQ files 
find ./ -type f -name "*.fastq.gz" 

# extract sample names
| while read F; do basename $F | cut -d _ -f 1,3 ; done 

# sort & find unique identifiers
| sort | uniq | 

# concatenate files for each identifer (R1)
while read P; do find ./ -type f -name "${P}_*R1*.fastq.gz"  -exec cat '{}' ';'  > ${P}.R1.merged.fq.gz ; done

# concatenate files for each identifer (R2)
while read P; do find ./ -type f -name "${P}_*R2*.fastq.gz"  -exec cat '{}' ';'  > ${P}.R2.merged.fq.gz ; done

## all together, the above code can be run on the command line like this: 

# R1
find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | cut -d _ -f 1,3 ; done | sort | uniq | while read P; do find ./ -type f -name "${P}_*R1*.fastq.gz"  -exec cat '{}' ';'  > ${P}.R1.merged.fq.gz ; done

# R2
find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | cut -d _ -f 1,3 ; done | sort | uniq | while read P; do find ./ -type f -name "${P}_*R2*.fastq.gz"  -exec cat '{}' ';'  > ${P}.R2.merged.fq.gz ; done

## relocate merged files to their own, new directory
mkdir -p /gpfs/gibbs/pi/huckins/ekw28/ghre_lep_sema_RNAseq/sample_dir/merged_all_lanes_NEW
mv *.{R1,R2}.merged.fq.gz /gpfs/gibbs/pi/huckins/ekw28/ghre_lep_sema_RNAseq/sample_dir/merged_all_lanes_NEW

