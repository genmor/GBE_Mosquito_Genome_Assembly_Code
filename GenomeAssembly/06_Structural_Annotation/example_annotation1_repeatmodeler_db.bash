#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=100:00:00
#SBATCH --mem=5G

module load gcc/10.2.0
PATH=/work/soghigian_lab/apps/conda/envs/RepeatModeler2/bin:$PATH

#activate the environment first, then submit the job!

#conda activate /work/soghigian_lab/apps/conda/envs/RepeatModeler2/

#ensure that the genome assembly has "simple" fasta headers.
#fasta headers with lots of spaces seems to mess up pipeline

# sed 's/ .*//' GCF_000001735.4_TAIR10.1_genomic.fna >tmp.fna
# mv tmp.fna GCF_000001735.4_TAIR10.1_genomic.fna

mkdir rep_db
cp PATH/TO/GENOME/ASSEMBLY .
asm=assembly.fasta
sed -i '/^>/s/ .*//' assembly.fasta


BuildDatabase -name ./rep_db/${asm%.fa*} ${asm}