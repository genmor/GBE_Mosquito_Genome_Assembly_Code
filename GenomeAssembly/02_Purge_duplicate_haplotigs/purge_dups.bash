#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=100:00:00
#SBATCH --mem=90G

####### Set environment variables ###############
module load gcc/10.2.0
PATH=/work/soghigian_lab/apps/purge_dups/bin:$PATH
PATH=/work/soghigian_lab/apps/minimap2-2.23_x64-linux:$PATH
####### Run your script #########################

#conda activate /work/soghigian_lab/apps/conda/envs/purge_dups

cd /scratch/${SLURM_JOB_ID}
echo check /scratch/${SLURM_JOB_ID}

#copy reads, assembly, and a copy of this script to scratch
#I do this because for some reason, purge_dups creates *.paf files in the directory with the reads, which we don't need to keep after this
cp /work/soghigian_lab/data/pacbio/maryland1/BF05/deepcon/filtered/*.fastq.gz /scratch/${SLURM_JOB_ID} #reads


cp /work/soghigian_lab/gen.morinaga/mosquito_hifi_new2/Bf05/Bf05_canu_out/Bf05.contigs.fasta /scratch/${SLURM_JOB_ID} #assembly

#set variable names
pb_list=*.fastq.gz #list of reads
pri_asm=Bf05.contigs.fasta #assembly to be purged


###########STEP 1###########
#map assembly to reads using loop
for i in $pb_list;
do
	minimap2 -I6G -t 20 -x map-hifi $pri_asm $i | gzip -c - > $i.paf.gz #I is the number of bases loaded into ram for indexing (def=4G). I've never messed with this. t is number of threads -x is the mapping mode (i.e., hifi, ont, etc...)
done

#create a read depth histogram
pbcstat *.paf.gz

#programatically figure out where hetero/homozygous peaks are in the histogram
calcuts PB.stat > cutoffs 2>calcults.log

#plot histogram. take a look at this at the png and see if you agree with the cutoffs that purge_dups is picking. 
/work/soghigian_lab/apps/purge_dups/scripts/hist_plot.py -c cutoffs PB.stat PB.cov.png

###########STEP 2###########
#you need to run this after STEP 1 finishes either of the following chunks depending on what your calcults.log file says. If it says something like the following, use Chunk 1
# [M::calcuts] Find 2 peaks
# [M::calcuts] Merge local peaks and valleys: 2 peaks remain
# [M::calcuts] Remove peaks and valleys less than 5: 2 peaks remain
# [M::calcuts] Use top 3 frequent read depth
# [M::calcuts] Found a valley in the middle of the peaks, use two-peak mode

#If calcults.log says the following, use Chunk 2
#this means that purge_dups wasn't able to determine where the cut-offs should be for low, medium, or high coverage
#take a look at the histogram and try to work it out
# [M::calcuts] Find 1 peaks
# [M::calcuts] mean: 22, peak: 21, mean larger than peak, treat as diploid assembly
# [W::calcuts] mean is not significantly different with peak, please recheck the cutoffs
#


#####Chunk 1#####
# split_fa $pri_asm > $pri_asm.split
# minimap2 -xasm5 -DP -t 10 $pri_asm.split $pri_asm.split | gzip -c - > $pri_asm.split.self.paf.gz
# purge_dups -2 -T cutoffs -c PB.base.cov $pri_asm.split.self.paf.gz > dups.bed 2> purge_dups.log
# get_seqs -e dups.bed $pri_asm -p Nmex03_canu_pd-man

#####Chunk 2#####
# calcuts -l -m -u PB.stat > cutoffs2 2>calcults2.log #set cutoffs manually -l is low cutoff, -m is middle cutoff, -u is upper cutoff.
# split_fa $pri_asm > $pri_asm.split
# minimap2 -xasm5 -DP -t 10 $pri_asm.split $pri_asm.split | gzip -c - > $pri_asm.split.self.paf.gz
# purge_dups -2 -T cutoffs2 -c PB.base.cov $pri_asm.split.self.paf.gz > dups.bed 2> purge_dups.log
# get_seqs -e dups.bed $pri_asm -p Nmex03_canu_pd-man