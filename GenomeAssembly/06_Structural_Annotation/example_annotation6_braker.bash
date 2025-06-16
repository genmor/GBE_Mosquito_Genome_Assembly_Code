#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=35
#SBATCH --time=50:00:00
#SBATCH --mem=80G


##this worked!
export braker=/work/soghigian_lab/apps/braker3/braker3_latest.sif
export GENEMARK_PATH=/work/soghigian_lab/apps/braker3/GeneMark-ETP/bin

module load perl/5.32.0
module load mysql/8.0.30
module load lib/boost/1.70.0-mpi

#move to scratch
cd /scratch/${SLURM_JOB_ID}
echo check /scratch/${SLURM_JOB_ID}

echo if the analysis ran without errors copy: /scratch/${SLURM_JOB_ID}/braker_out/braker
#create output directory for braker and go there
mkdir -p braker_out
cd braker_out

#copy genome, the protein database, and config files for AUGUSTUS
cp /work/soghigian_lab/gen.morinaga/mosquito_hifi_new2/MadSpNov/MadSpNov_canu/ph/strut_annotation/rep_mask_combined/MadSpNov_rep_mask_combined_masked.fasta .
cp /work/soghigian_lab/gen.morinaga/mosquito_hifi_new2/Arthropoda.fa .
cp -r /work/soghigian_lab/apps/Augustus-3.4.0/config config


export TMPDIR=$(pwd)/braker_tmp 
    mkdir -p ${TMPDIR}

export genome="./MadSpNov_rep_mask_combined_masked.fasta"
export proteins="./Arthropoda.fa"
export AUGUSTUS_CONFIG_PATH="/scratch/${SLURM_JOB_ID}/braker_out/config"

# run braker
# paths need to be defined within apptainer, hence -B followed by directory variables


apptainer exec -B ${PWD}:${PWD},${HOME},${GENEMARK_PATH},${AUGUSTUS_CONFIG_PATH} ${braker} braker.pl --genome=${genome} --prot_seq=${proteins} --DIAMOND_PATH=/work/soghigian_lab/apps/diamond --GENEMARK_PATH=${GENEMARK_PATH} --AUGUSTUS_CONFIG_PATH=${AUGUSTUS_CONFIG_PATH} --softmasking --threads 35 --gff3 --species MadSpNov
