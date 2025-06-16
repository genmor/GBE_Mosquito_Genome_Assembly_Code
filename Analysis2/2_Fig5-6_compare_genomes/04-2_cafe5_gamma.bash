#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --time=5:00:00
#SBATCH --mem=80G

####### Set environment variables ###############
module load gcc/10.2.0
PATH=/work/soghigian_lab/apps/CAFE5/bin:$PATH

tree=../SpeciesTreeAlignment.fa.timetree.nwk
dat=../OrthogroupCounts_Interpro_var_filt.tsv
ltree=lambda_tree.tre

#i is input data table MAKE SURE THAT IT HAS UNIX LINE ENDINGS AND THE FIRST COLUMN IS A GENE DESCRIPTION (CAN BE ALL NA)
#-t is ultrametric newick tree. ensure that node lables have been stripped out of the tree
#-y is a newick tree wherein branches are grouped to specify different lambdas
#-o is output directory
#-s is the number of simulations (needs to have no space) or --simulate = 1000
#-p uses a Poisson distribution for the root frequency distribution; no flag uses uniform distribution
#a value can be given for -p (-p10) otherwise estimated from gene families


####estimate an error model first
output_dir=cafe5_gamma_out

mkdir $output_dir
cd $output_dir
echo start cafe5 error run
cafe5 -i $dat -t $tree -o error -e -p -P 0.01  -c 40 &> error.log
mv error.log ./error


####barebones, base model cafe5 analysis can be done by

# cafe5 -i $dat -t $tree -o out -P 0.01


####if you want to estimate different gamma categories, documentation suggests running several values of k and each one multiple times
####to ensure convergence. Likelihoods cannot be analytically compared between k. A good rule of thumb is to look for likelihoods that
####appear to stabilize across runs for each k. Note that the base model will very likely have the highest likelihood and be VERY stable.

####estimate an optimal k (from 1[base] to 4) across 5 runs
echo start cafe5 gamma analyses
for i in {1..4}
	do
	mkdir ./glam_k_${i};
		for j in {1..5}
		do
			(cafe5 -i $dat -t $tree -p -e./error/Base_error_model.txt -o ./glam_k_${i}/run_${j} -k $i -P 0.01 -c 8 &> glam_k_${i}_run_${j}.log)&
		done
		wait
	done
	