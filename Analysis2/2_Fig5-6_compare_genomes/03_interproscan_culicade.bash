#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --time=12:00:00
#SBATCH --mem=120G

module load java/openjdk-23.0.1
PATH=/work/soghigian_lab/apps/interproscan/interproscan-5.67-99.0:$PATH

# interactively create these directories
# mkdir tmp
# mkdir scripts
# mkdir slurm_out
# mkdir gff
# mkdir tsv


#-t to scan nucleotide sequences
#-dp disables match precalculation
#-appl <appl name> different applications to run if multiple are desired, serparate using comma, e.g., CDD,COILS,PANTHER
#different applications are: CDD COILS Gene3D HAMAP MOBIDB PANTHER Pfam IRSF PRINTS PROSITESFLD SMART SUPERFAMILY NCBIFAM
#-i input sequences
#-goterms output GO mappings
#-b basename for output
#-o output name. must provide output format output formats: TSV XML JSON GFF3
#-pa map pathway info
#-dra disable residue annotation. Disable if you dont need this info as it can improve performance
#
# ln -s /work/soghigian_lab/gen.morinaga/interproscan/orthofinder/faa/OrthoFinder/Results_Oct23/Orthogroup_Sequences/OG0000000.fa
# og_seq_dir=/work/soghigian_lab/gen.morinaga/compare_genomes5/comparative_analysis/compare_genomes3/culicidae6/Results_Apr14/Orthogroup_Sequences/

# input=TEMPLATE



#the gff3 output contains amino acid sequences in fasta format after the last row
# interproscan.sh -i ${og_seq_dir}${input}.fa* -goterms -dp -b $input -T tmp -cpu 12 -appl PANTHER -f TSV, GFF3
# mv slurm-${SLURM_JOB_ID}.out ./slurm_out
# mv ${input}.gff3 ./gff
# mv ${input}.tsv ./tsv

mkdir og_sequences_renamed
echo add taxon name to sequence headers
for file in ../Orthogroup_Sequences/*.fa
	do
	name=$(basename $file .fa)
	sed "s|^>\(.*\)|>${name}_\1|I" ${file} >./og_sequences_renamed/${name}.fa
done


echo concatenate sequences into single amino acid fasta
for file in ./og_sequences_renamed/*.fa
do
cat $file >> ./all_ogs.faa
done

echo start interproscan
interproscan.sh -i all_ogs.faa -goterms -dp -b culicidae -T tmp -cpu 25 -appl PANTHER -f TSV, GFF3
