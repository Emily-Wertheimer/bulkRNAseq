## job submission for align parallel script

#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --mem=200g
#SBATCH -t 23:00:00
#SBATCH --cpus-per-task=10
#SBATCH -o alignParallel.log

ml R/4.3.0-foss-2020b
Rscript alignParallel.R
