#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --time=7-00:00:00
#SBATCH --mem=90G

####### Set environment variables ###############
module load lib/zlib/1.2.11
module load gcc/10.2.0
module load perl/5.32.0
PATH=/work/soghigian_lab/apps/R-4.1.2/bin:$PATH
PATH=/work/soghigian_lab/apps/conda/envs/orthofinder/bin:$PATH
PATH=/work/soghigian_lab/apps/iqtree-latest/bin:$PATH
#conda activate /work/soghigian_lab/apps/conda/envs/orthofinder


#-M either msa or dendroblast defaults to dendroblast
#-y splits paralog clades below the root into separate HOGS
#-T if -M msa, choose between iqtree fasttree raxml and raxml-ng Note that this will add time to the analysis
#-t number of threads for sequence search
#-a number of parallel threads...? not sure that this helps with anything.
# dir=/work/soghigian_lab/gen.morinaga/compare_genomes5/new_preprocess/5_simple_faa_headers/
#flags like f* and o* are different analysis starts and ends points. 
#-s to start analysis with a tree inferred outside of orthofinder
# orthofinder -f $dir -t 30 -a 5 -o culicidae3 -y -M msa -T iqtree
# orthofinder -fg $fg -t 30 -a 6 -y -M msa -T iqtree
# orthofinder -f $dir -t 30 -a 5 -o culicidae -y -M msa -A mafft -oa
# dir=/work/soghigian_lab/gen.morinaga/compare_genomes5/comparative_analysis/compare_genomes2/culicidae/Results_Feb21
#orthofinder -ft $dir  -t 30 -a 5 -y -M msa -A mafft -T iqtree -s $tree
dir=/work/soghigian_lab/gen.morinaga/compare_genomes5/comparative_analysis/compare_genomes3/culicidae9/faa
orthofinder -f $dir -M msa -A mafft -t 30 -a 15 -o culicidae9
