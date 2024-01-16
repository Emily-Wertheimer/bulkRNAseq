#!/bin/sh
#Goal: This script checks the qualitiy of our fastq files and performs an alignment to the human cDNA transcriptome reference with Kallisto.
#Meilin Fernandez Garcia 
#03-07-2023

#To run this 'shell script' you will need to open your terminal and navigate to the directory where this script resides on your computer.
# This should be the same directory where you fastq files and reference fasta file are found.
# Then type './readMapping.sh' (without the quotes) at the prompt.  
# This will begin the process of running each line of code in the shell script.


#load modules:
ml kallisto/0.45.0
ml FastQC/0.11.9-Java-11
ml MultiQC/1.10.1-foss-2020b-Python-3.8.6

# first use fastqc to check the quality of our fastq files:
fastqc *.gz -t 4

#Build and index: 
#building an index takes a long time so I downloaded the kallisto index built with ENSEMBL from :
wget https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/94/Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx.gz
#This index file was produced using kallisto version 0.45.0 and Ensembl Transcriptomes v94
#unzip index
gunzip Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx.gz

# now make  map reads to the indexed reference host transcriptome
# use as many 'threads' as your machine will allow in order to speed up the read mapping process.
# note that we're also including the '&>' at the end of each line
# this takes the information that would've been printed to our terminal, and outputs this in a log file that is saved in /data/course_data

kallisto quant -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx -o 3AsBSA1 -t 8 3182iAstro_BSA_1_096_097_S1_L001_R1_001.fastq.gz 3182iAstro_BSA_1_096_097_S1_L001_R2_001.fastq.gz &> 3AsBSA1.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx -o 3AsBSA2 -t 8 3182iAstro_BSA_2_084_109_S2_L001_R1_001.fastq.gz 3182iAstro_BSA_2_084_109_S2_L001_R2_001.fastq.gz &> 3AsBSA2.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx -o 3AsBSA3 -t 8 3182iAstro_BSA_3_072_121_S3_L001_R1_001.fastq.gz 3182iAstro_BSA_3_072_121_S3_L001_R2_001.fastq.gz &> 3AsBSA3.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx -o 3AsA11 -t 8 3182iAstro_A1_1_060_133_S4_L001_R1_001.fastq.gz 3182iAstro_A1_1_060_133_S4_L001_R2_001.fastq.gz &> 3AsA11.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx -o 3AsA12 -t 8 3182iAstro_A1_2_048_145_S5_L001_R1_001.fastq.gz 3182iAstro_A1_2_048_145_S5_L001_R2_001.fastq.gz &> 3AsA12.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx -o 3AsA13 -t 8 3182iAstro_A1_3_036_157_S6_L001_R1_001.fastq.gz 3182iAstro_A1_3_036_157_S6_L001_R2_001.fastq.gz &> 3AsA13.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx -o 3AsA21 -t 8 3182iAstro_A2_1_024_169_S7_L001_R1_001.fastq.gz 3182iAstro_A2_1_024_169_S7_L001_R2_001.fastq.gz &> 3AsA21.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx -o 3AsA22 -t 8 3182iAstro_A2_2_012_181_S8_L001_R1_001.fastq.gz 3182iAstro_A2_2_012_181_S8_L001_R2_001.fastq.gz &> 3AsA22.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx -o 3AsA23 -t 8 3182iAstro_A2_3_095_098_S9_L001_R1_001.fastq.gz 3182iAstro_A2_3_095_098_S9_L001_R2_001.fastq.gz &> 3AsA23.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx -o 5AsBSA1 -t 8 553iAstro_BSA_1_083_110_S10_L001_R1_001.fastq.gz 553iAstro_BSA_1_083_110_S10_L001_R2_001.fastq.gz &> 5AsBSA1.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx -o 5AsBSA2 -t 8 553iAstro_BSA_2_071_122_S11_L001_R1_001.fastq.gz 553iAstro_BSA_2_071_122_S11_L001_R2_001.fastq.gz &> 5AsBSA2.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx -o 5AsBSA3 -t 8 553iAstro_BSA_3_059_134_S12_L001_R1_001.fastq.gz 553iAstro_BSA_3_059_134_S12_L001_R2_001.fastq.gz &> 5AsBSA3.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx -o 5AsA11 -t 8 553iAstro_A1_1_047_146_S13_L001_R1_001.fastq.gz 553iAstro_A1_1_047_146_S13_L001_R2_001.fastq.gz &> 5AsA11.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx -o 5AsA12 -t 8 553iAstro_A1_2_035_158_S14_L001_R1_001.fastq.gz 553iAstro_A1_2_035_158_S14_L001_R2_001.fastq.gz &> 5AsA12.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx -o 5AsA13 -t 8 553iAstro_A1_3_023_170_S15_L001_R1_001.fastq.gz 553iAstro_A1_3_023_170_S15_L001_R2_001.fastq.gz &> 5AsA13.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx -o 5AsA21 -t 8 553iAstro_A2_1_011_182_S16_L001_R1_001.fastq.gz 553iAstro_A2_1_011_182_S16_L001_R2_001.fastq.gz &> 5AsA21.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx -o 5AsA22 -t 8 553iAstro_A2_2_094_099_S17_L001_R1_001.fastq.gz 553iAstro_A2_2_094_099_S17_L001_R2_001.fastq.gz &> 5AsA22.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx -o 5AsA23 -t 8 553iAstro_A2_3_082_111_S18_L001_R1_001.fastq.gz 553iAstro_A2_3_082_111_S18_L001_R2_001.fastq.gz &> 5AsA23.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx -o uHABSA1 -t 8 unknown_HA_BSA_1_070_123_S19_L001_R1_001.fastq.gz unknown_HA_BSA_1_070_123_S19_L001_R2_001.fastq.gz &> uHABSA1.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx -o uHABSA2 -t 8 unknown_HA_BSA_2_058_135_S20_L001_R1_001.fastq.gz unknown_HA_BSA_2_058_135_S20_L001_R2_001.fastq.gz &> uHABSA2.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx -o uHABSA3 -t 8 unknown_HA_BSA_3_046_147_S21_L001_R1_001.fastq.gz unknown_HA_BSA_3_046_147_S21_L001_R2_001.fastq.gz &> uHABSA3.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx -o uHAA11 -t 8 unknown_HA_A1_1_034_159_S22_L001_R1_001.fastq.gz unknown_HA_A1_1_034_159_S22_L001_R2_001.fastq.gz &> uHAA11.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx -o uHAA12 -t 8 unknown_HA_A1_2_022_171_S23_L001_R1_001.fastq.gz unknown_HA_A1_2_022_171_S23_L001_R2_001.fastq.gz &> uHAA12.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx -o uHAA13 -t 8 unknown_HA_A1_3_010_183_S24_L001_R1_001.fastq.gz unknown_HA_A1_3_010_183_S24_L001_R2_001.fastq.gz &> uHAA13.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx -o uHAA21 -t 8 unknown_HA_A2_1_093_100_S25_L001_R1_001.fastq.gz unknown_HA_A2_1_093_100_S25_L001_R2_001.fastq.gz &> uHAA21.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx -o uHAA22 -t 8 unknown_HA_A2_2_081_112_S26_L001_R1_001.fastq.gz unknown_HA_A2_2_081_112_S26_L001_R2_001.fastq.gz &> uHAA22.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx -o uHAA23 -t 8 unknown_HA_A2_3_069_124_S27_L001_R1_001.fastq.gz unknown_HA_A2_3_069_124_S27_L001_R2_001.fastq.gz &> uHAA23.log

# summarize fastqc and kallisto mapping results in a single summary html using MultiQC
multiqc -d .

echo "Finished"

exit 
