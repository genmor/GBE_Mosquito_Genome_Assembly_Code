#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --time=12:00:00
#SBATCH --mem=80G

####### Set environment variables ###############
module load lib/zlib/1.2.11
module load gcc/10.2.0
module load perl/5.32.0
export PATH=/work/soghigian_lab/apps/conda/envs/ins/Inspector:$PATH

#conda activate /work/soghigian_lab/apps/conda/envs/ins/

####change working directory to scratch
cd /scratch/${SLURM_JOB_ID}
echo check /scratch/${SLURM_JOB_ID}

####copy assenmbly to scratch
ln -s /work/soghigian_lab/data/pacbio/maryland1/BF05/deepcon/filtered/BF05_A01.filt.fastq.gz
ln -s /work/soghigian_lab/data/pacbio/maryland1/BF05/deepcon/filtered/BF05_B01.filt.fastq.gz
ln -s /work/soghigian_lab/gen.morinaga/mosquito_hifi_new2/Bf05/flye/Bf05_flye.fasta

####set variable names. 
#asm is the assembly to inspect
#prefix is the prefix given to output directory
asm=Bf05_flye.fasta
reads=*.fastq.gz

####this is the "inspection" part of inspector. 
#-c is the assembly to be inspected
#-r is the reads
#-d is how the reads were generated e.g., hifi, ont, clr...
#-o is the output directory name
#-t is the number of threads
python2 /work/soghigian_lab/apps/conda/envs/ins/Inspector/inspector.py -c ${asm} \
-r $reads \
-d hifi \
-o ${asm%.fa*} \
-t 30

mkdir ./${asm%.fa*}/out
cp ./${asm%.fa*}/{*.bed,summary_statistics,contig_length_info} ./${asm%.fa*}/out
####this is the "correction" part of inspector.
#-i is the output directory that the "inspection" part created
#--datatype is hoow the reads were generated e.g., hifi, ont, clr
#-o is the name of the output directory. here, differentiated by "-cor" from the "inspection" part
#-t is the number of threads
python2 /work/soghigian_lab/apps/conda/envs/ins/Inspector/inspector-correct.py -i ${asm%.fa*} \
--datatype pacbio-hifi \
-o ${asm%.fa*}_insp \
-t 30