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

mkdir repeat_masker
cd repeat_masker
ln -s /work/soghigian_lab/gen.morinaga/mosquito_hifi_new2/TH21/annotation/rep_db/rep_classifier/classified_repeats_for_repeatmasker.fa
replib=classified_repeats_for_repeatmasker.fa

ln -s PATH/TO/GENOME/ASSEMBLY #previous scripts assembly.fasta

#Round 1 mask simple repeats first. Takes about an hour; used less than 10 Gb memory

asm_r1=assembly.fasta
RepeatMasker -pa 25 -a -dir rep_mask_r1 -noint -noisy -xsmall ${asm_r1}
ln -s ./rep_mask_r1/${asm_r1}.masked

#Round 2 mask elements founds in Diptera. Takes less than an hour; used less than 10 Gb memory

asm_r2=${asm_r1}.masked
RepeatMasker -pa 25 -a -dir rep_mask_r2 -nolow -species diptera -noisy -xsmall ${asm_r2}
ln -s ./rep_mask_r2/${asm_r2}.masked

#Round 3 mask elements from combined de novo and TEfam repeat libraries
replib=classified_repeats_for_repeatmasker.fa
asm_r3=${asm_r2}.masked
RepeatMasker -pa 25 -a -dir rep_mask_r3 -nolow  -lib ${replib} -noisy -xsmall ${asm_r3}
