#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --time=100:00:00
#SBATCH --mem=75G


#this job used 51G/150G of RAM 

module load gcc/10.2.0
PATH=/work/soghigian_lab/apps/conda/envs/RepeatModeler2/bin:$PATH

#activate the environment first, then submit the job!
#conda activate /work/soghigian_lab/apps/conda/envs/RepeatModeler2/

#rep is the path + the prefix to the database created in step1
#fasta must be in the currrent working directory
#-ucsctools_dir is the directory where these tools are located
asm=assembly.fasta
rep=./rep_db/${asm%.fa*}
ucsc=/work/soghigian_lab/apps/conda/envs/RepeatModeler2/ucsc_utils #not sure this does anything...
RepeatModeler -threads 25 -database ${rep} -LTRStruct -ucsctools_dir ${ucsc}