#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --time=100:00:00
#SBATCH --mem=180G

####### Set environment variables ###############
module load gcc/10.2.0
PATH=/work/soghigian_lab/apps/conda/envs/ragtag/bin:$PATH
PATH=/work/soghigian_lab/apps/minimap2-2.23_x64-linux:$PATH
PATH=/work/soghigian_lab/apps/bioawk:$PATH
PATH=/work/soghigian_lab/apps/samtools-1.15/bin:$PATH
PATH=/work/soghigian_lab/apps/conda/envs/syri/bin:$PATH
PATH=/work/soghigian_lab/apps/seqtk:$PATH
####### Run your script #########################
cd /scratch/${SLURM_JOB_ID}
echo check /scratch/${SLURM_JOB_ID}


#activate syri before submitting the job!
#conda activate /work/soghigian_lab/apps/conda/envs/syri

#copy the L5 Aedes aegypti chromosomal genome to scratch
ln -s /work/soghigian_lab/gen.morinaga/Aeaeg_L5-chr-only.fasta
ln -s /work/soghigian_lab/gen.morinaga/mosquito_hifi_new2/Bf05/assemble_from_dc_reads/Bf05-dc_canu/Bf05-dc_canu_ph/Bf05-dc_canu_ph_insp/Bf05_dc_canu_ph_insp_yahs_redo/Bf05_yahs-redo/tgs3/2024-02-11_Bf05-dc_canu_ph_yahs_jbat_tgs.fa

#create a list of sequence names for the target genome that you want to keeps
#in this case, scaffold_1 scaffold_2 and scaffold_3
cp /work/soghigian_lab/gen.morinaga/mosquito_hifi_new2/dc-assembly_synteny/scafs.txt .

bf_asm=2024-02-11_Bf05-dc_canu_ph_yahs_jbat_tgs.fa

xargs samtools faidx ${bf_asm} < scafs.txt > bf_scafs1.2.3.fa


ref=Aeaeg_L5-chr-only.fasta
bf=bf_scafs1.2.3.fa

#the bf scaffold 3 needs to be reverse complemented
#first split each scaffold into individual fastas
cat ${bf} |\
awk '/^>/ {if(N>0) printf("\n"); printf("%s\n",$0);++N;next;} { printf("%s",$0);} END {printf("\n");}' |\
split -l 2 --additional-suffix=.fasta - seq_

seqtk seq -r seq_ac.fasta > scaf3i.fasta
cat seq_aa.fasta seq_ab.fasta scaf3i.fasta > bf_scafs1.2.3i.fasta

bf=bf_scafs1.2.3i.fasta

#define variables

echo run minimap2
#use minimap to align the renamed reference and reference scaffolded assembly and output as *.sam
#-ax will need to be changed depending on how distantly related the reference species and the species of interest are. if it's close (or the same) just use asm5. if more distant, try asm20
minimap2 -ax asm5 --eqx -t 25 ${ref} ${bf} > ref_Bf05_syri.bam

#create synteny map from the *.sam
#-c is output from minimap 
#-r is the reference chromosome 
#-q is the is your assembly
#-k keeps intermediate files default is false
#-F is input file type. T=table, S=sam, B=bam defaults to T
#--prefix is the prefix for file names
#--nc is the number of cores to use for analysis

echo run syri

syri -c ref_Bf05_syri.bam -r ${ref} -q ${bf} -F S -f --prefix ref_bf_ --nc 25

#do this in an interactive session with the syri environment activated
#--sr is the output from syri. has extension *.out
#you can specify --sr multiple times if you have multiple *.out files you want to compare but they have to be in a specific order when you do the comparison and they have to listed in exactly the same order in genome.txt file
#-W is the width of the plot in inches
#-H is height of the plot in inchaes
#genomes.txt is a tab-separated-file with columns file and name. file is the filename of the assembly, name is whatever you want to call that assembly file in the plot
#make sure that the first line is commented out (i.e., #file)
#optionally, genome.txt can have an additional column called tags. for use see https://github.com/schneebergerlab/plotsr

#plotsr --sr syri.out --genomes genomes.txt -o output_plot.png -W 11 -H 8