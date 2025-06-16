#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=12:00:00
#SBATCH --mem=120G

####### Set environment variables ###############
module load lib/zlib/1.2.11
module load gcc/10.2.0
module load perl/5.32.0
PATH=/work/soghigian_lab/apps/minimap2-2.23_x64-linux:$PATH
export PATH=/work/soghigian_lab/apps/R-4.1.2/bin:$PATH
#conda activate /work/soghigian_lab/apps/conda/envs/purge_haplotigs

cd /scratch/${SLURM_JOB_ID}
echo check /scratch/${SLURM_JOB_ID}


cp /work/soghigian_lab/data/pacbio/maryland1/BF05/deepcon/filtered/*.fastq.gz /scratch/${SLURM_JOB_ID} #reads


cp /work/soghigian_lab/gen.morinaga/mosquito_hifi_new2/Bf05/Bf05_canu_out/Bf05.contigs.fasta /scratch/${SLURM_JOB_ID} #assembly

#set variables
asm=Bf05.contigs.fasta

#use minimap2 to assembly to reads
#-a generates CIGAR and outputs to .sam rather than .paf -x is the read type (e.g., hifi, ont, etc) -t is the number of threads
#output from minimap is piped to samtools view and sort functions which generates bam files. 
#samtool view -h outputs headers -F prevents output of alignments with bit value
#samtool sort sorts alignment  -@ is the number of threads  -m is max memory per thread
#-o is output name -T is the temporary alignment prefix
minimap2 -t 20 -ax map-hifi ${asm} ${reads} \
| samtools view -hF 256 - \
|samtools sort -@ 20 -m 4G -o ${asm%.fa*}_aligned.bam -T Nmex03_tmp.ali

#generates read depth histogram -b is the bam you generated in the previous step -g is the assembly  -t is the number of threads
purge_haplotigs readhist -b ${asm%.fa*}_aligned.bam \
-g ${asm} -t 20