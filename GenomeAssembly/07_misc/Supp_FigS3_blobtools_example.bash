#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --time=72:00:00
#SBATCH --mem=160G

####### Set environment variables ###############
#Load conda into your local environment before submitting this job, as below:

#remember to activate the conda environment first!
#conda activate /work/soghigian_lab/apps/conda/envs/sepp

PATH=/work/soghigian_lab/apps/ncbi-blast-2.12.0+/bin:$PATH
PATH=/work/soghigian_lab/apps/diamond:$PATH
PATH=/work/soghigian_lab/apps/seqtk:$PATH
PATH=/work/soghigian_lab/apps/minimap2-2.23_x64-linux:$PATH
PATH=/work/soghigian_lab/apps/bbmap:$PATH


#####prior to running, cat reads together into single file make sure that the extension is fastq.gz
#cat read1.fq.gz reads2.fq.gz > all_reads.fq.gz

#####get the number of scaffolds in the assembly and the scaffold bp length
#stats.sh in=assembly.fasta format=3 > bbmap.txt

#####get the number of bases in the reads file (this just prints to console so pay attention, also kinda slow)
#zcat file.fq.gz | paste - - - - | cut -f2 | wc -c


#####or check the report output from canu

#####copy the assembly to the working directory and gzip it and change the file extension specifically to fasta.gz 
#cp assembly.fa .
#gzip assembly.fa

#set directory and file names
dir=/working/dir/

name=whatever

btk pipeline run --tool blobtoolkit --config config_abs.yaml --threads 40 --workdir ${dir}/${name}

#Very likely the figure generation will fail. Instead, try the below commands:
#Note that if you save as a PNG, fonts may not be visible, depending on your environment. So just choose SVG!

blobtk plot -v snail -d ${dir}/${name}/view/${name} -o ${name}_snail.svg
blobtk plot -v blob -d ${dir}/${name}/view/${name} -o ${name}_blob.svg
blobtk plot -v cumulative -d ${dir}/${name}/view/${name} -o ${name}_cumulative.svg