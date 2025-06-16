#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=12:00:00
#SBATCH --mem=80G

####### Set environment variables ###############
module load lib/zlib/1.2.11
module load gcc/10.2.0
module load perl/5.32.0
PATH=/work/soghigian_lab/apps/minimap2-2.23_x64-linux:$PATH
export PATH=/work/soghigian_lab/apps/R-4.1.2/bin:$PATH
#conda activate /work/soghigian_lab/apps/conda/envs/purge_haplotigs

asm=Bf05.contigs.fasta
#run this script where you have all the files generated from step 1
#sets the cutoffs homo/heterozygous coverage -l is low, -m middle, -h high
purge_haplotigs contigcov -i ${asm%.fa*}_aligned.bam.gencov -l 5 -m 31 -h 110

#the purge thing step -b is bam file, -g assembly, -o is output name, -c the read coverage stats file 
#-d generates dotplots -a % cutoff for to identify contig as haplotig
purge_haplotigs purge -b ${asm%.fa*}_aligned.bam -g \
${asm} \
-o ${asm%.fa*}_ph -c coverage_stats.csv -a 60