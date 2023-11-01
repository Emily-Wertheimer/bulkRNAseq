#! /usr/bin/python
​
'''
Sep 20 2023
Chia-Yi Joy Lee (cl2375)
filter counts >20
'''
'''
module load Python/3.8.6-GCCcore-10.2.0
module load SciPy-bundle/2020.11-foss-2020b (for pandas)
python filter_counts.py -I /home/cl2375/palmer_scratch/RNAseq_LentiInflam_Abeta/readcounts/ -O /home/cl2375/palmer_scratch/RNAseq_LentiInflam_Abeta/filtercounts/
'''
import argparse
import os
import pandas as pd
import numpy as np
​
# take in the input directory and output directory
parser = argparse.ArgumentParser()
parser.add_argument("-I", type=str, help="Input directory")
parser.add_argument("-O", type=str, help="output directory")
# set up the input directory and output directory
args = parser.parse_args()
input_dir = args.I
output_dir = args.O
​
# all timepoints = ["iPS","HPC", "D4","D8","D16","D24","M8"]
# set up the function
def filter(input_dir):
​
	# read in meta file from featureCounts
	table = pd.read_csv(input_dir+"all_featureCounts.txt", sep = "\t", skiprows = 1, index_col = 0)
	# rename columns in meta file (2607_iPS_1_Aligned.sortedByCoord.out.bam --> 2607_iPS_1)
	name_list = []
	# rearrange column names
	for col in table.columns:
		if col.endswith("_Aligned.sortedByCoord.out.bam"):
			col_new = col.replace("_Aligned.sortedByCoord.out.bam","")
			name_list.append(col_new)
		else:
			name_list.append(col)
	# make the columns of table into rearranged names
	table.columns = name_list
	# select the column for expression and check if bigger than 20
	table_select = np.sum(table.iloc[:,5:] > 20, axis=1) > 0
	# make the dataframe containing boolean information into a list
	table_index = table_select.tolist() 
	# use the list as an index to obtain the subset of data > 20
	table_sub = table.loc[table_index,:]
​
	# output the file of gene id and gene counts
	table_sub.index.to_series().to_csv(output_dir + "/" + "filtered_genes.txt", sep='\t', index=False)
	table_sub.iloc[:,5:].to_csv(output_dir + "/" + "filtered_genes_counts.txt", sep = '\t')
	
# call the function
filter(input_dir)
