#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem=90G

####### Set environment variables ###############
module load gcc/10.2.0
PATH=/work/soghigian_lab/apps/conda/envs/TGS-GapCloser:$PATH
racon=/work/soghigian_lab/apps/racon/build/bin
PATH=/work/soghigian_lab/apps/seqtk:$PATH


#activate conda environment before submitting!
#conda activate /work/soghigian_lab/apps/conda/envs/TGS-GapCloser/

cd /scratch/${SLURM_JOB_ID}
echo check /scratch/${SLURM_JOB_ID}

#create soft links to assembly and reads
ln -s /path/to/assembly.fasta

ln -s /path/to/reads1.fastq.gz
ln -s /path/to/reads2.fastq.gz


echo concatenating reads
#concatenate read files together
cat read1.fastq.gz reads2.fastq.gz > reads.combined.fastq.gz

echo converting fastq reads to fasta reads
#convert fastq to fasta
seqtk seq -a reads.combined.fastq.gz > reads.combined.fasta.gz

#define variables
scaf_file=assmebly.fasta
reads=reads.combined.fasta
name=some_name

echo remove / from fasta headers
#convert / in the fasta file to _
sed -i 's!/!_!g' ${reads}


echo running tgs-gapcloser
#--o  is the output directory
#--ne for no error correction
#--racon will use racon to perform error correction (long reads)
#--pilon will use pilon to perform error correction (short reads)
#--pilon will require --ngs where ngs reads can be specified
#--tgstype pb or ont, ont is default

tgsgapcloser --scaff ${scaf_file} \
--reads ${reads} \
--o ${name}  \
--ne \
--minmap_arg '-x map-hifi' \
--tgstype hifi \
--thread 16 \
>pipe.log 2>pipe.err