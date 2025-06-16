#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --time=144:00:00
#SBATCH --mem=120G

####### Set Environment Variables #############

module load cmake
module load perl/5.32.0

#activate conda environment first!
PATH=/work/soghigian_lab/apps/iqtree-latest/bin:$PATH


####set -nt for multicore analysis!
dates=/work/soghigian_lab/gen.morinaga/compare_genomes5/comparative_analysis/compare_genomes3/culicidae9/dates.txt
iqtree2 -s /work/soghigian_lab/gen.morinaga/compare_genomes5/comparative_analysis/compare_genomes3/culicidae9/culicidae9/Results_May28/MultipleSequenceAlignments/SpeciesTreeAlignment.fa -m MFP --alrt 1000 -B 1000 -o Phlebotomus_papatasi -nt 25 --date $dates --date-tip 0