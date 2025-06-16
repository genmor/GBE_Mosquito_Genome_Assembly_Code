#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --time=100:00:00
#SBATCH --mem=250G

####### Set environment variables ###############
module load gcc/10.2.0
module load bioconda/conda3
source activate /work/soghigian_lab/apps/conda/envs/flye
export PATH="/work/soghigian_lab/apps/conda/envs/flye/bin:$PATH"

####### Run your script #########################
cd /scratch/${SLURM_JOB_ID}
echo check /scratch/${SLURM_JOB_ID}


###set taxon name
spp=Bf05

###set location of hifi reads
reads=/work/soghigian_lab/data/pacbio/maryland1/BF05/deepcon/filtered/*.fastq.gz

flye --pacbio-hifi $reads \
--genome-size 1.3g \
-o ${spp}_flye \
-t 40 \