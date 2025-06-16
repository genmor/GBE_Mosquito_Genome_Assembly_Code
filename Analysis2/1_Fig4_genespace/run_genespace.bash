#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --time=24:00:00
#SBATCH --mem=90G

####### Set environment variables ###############
module load lib/zlib/1.2.11
module load gcc/10.2.0
module load perl/5.32.0
PATH=/work/soghigian_lab/apps/conda/envs/R/bin:$PATH
PATH=/work/soghigian_lab/apps/conda/envs/orthofinder/bin:$PATH
#conda activate /work/soghigian_lab/apps/conda/envs/orthofinder

Rscript example_genespace.r

echo Finished...Consider opening the RData file in a local RStudio instance and editing the grobs/gtables