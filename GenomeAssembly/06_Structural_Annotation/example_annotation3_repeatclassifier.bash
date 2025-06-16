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
seqkit=/work/soghigian_lab/apps/seqkit
PATH=/work/soghigian_lab/apps/cdhit:$PATH

#activate the environment first, then submit the job!
#conda activate /work/soghigian_lab/apps/conda/envs/RepeatModeler2/


cd rep_db
mkdir rep_classifier #create a rep_classifier directory
cd rep_classifier
ln -s /work/soghigian_lab/gen.morinaga/mosquito_hifi_new2/TEfam1091plusAAEGY.consensi #link to TEfam repeat consensus file (it's just a fasta) to new directory
ln -s  ../*.fa #link to de novo repeat consensus from step 2 to new directory
asm=assembly.fasta #this is the assembly used just to get file name prefix
tefam=TEfam1091plusAAEGY.consensi #variable for TEfam repeat consensus


cat ${asm%.fa*}-families.fa $tefam > combined_consensi.fa
consensi=combined_consensi.fa
cd-hit-est -c 0.8 -n 5 -i $consensi -o combined_consensi_nonred.fa -T 25 -M 75000

cat combined_consensi_nonred.fa | $seqkit fx2tab | grep -v "Unknown" | $seqkit tab2fx > tefam_denovo_known.fa
cat combined_consensi_nonred.fa | $seqkit fx2tab | grep "Unknown" | $seqkit tab2fx > tefam_denovo_unknown.fa

unknown=tefam_denovo_unknown.fa

echo -e begin RepeatClassifier '\n'

#classify just the unknowns using RepeatClassifier
RepeatClassifier -consensi $unknown

#extract the newly identified set from the unknown set and concatenate them into the known set
cat ${unknown}.classified | $seqkit fx2tab | grep -v "Unknown" | $seqkit tab2fx >> tefam_denovo_known.fa

#separate the unknown set
cat ${unknown}.classified | $seqkit fx2tab | grep "Unknown" | $seqkit tab2fx > tefam_denovo_unknown_postclassification.fa

#combined unknown and known
cat tefam_denovo_known.fa tefam_denovo_unknown_postclassification.fa > classified_repeats_for_repeatmasker.fa