#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --time=100:00:00
#SBATCH --mem=75G


#with 25 cpus and 75G memory, the job took about 1.25 hours for Bf05
#cpu efficiency was low (3.98%) and memory efficiency with moderate (48%)

#request fewere cpus (like 4) next time. 

module load gcc/10.2.0
PATH=/work/soghigian_lab/apps/conda/envs/RepeatModeler2/bin:$PATH
PATH=/work/soghigian_lab/apps/bedtools2/bin:$PATH
PATH=/work/soghigian_lab/apps/seqtk:$PATH
PATH=/home/gen.morinaga/software/datamash-1.3/bin:$PATH

##activate conda environment first!!!
#conda activate /work/soghigian_lab/apps/conda/envs/RepeatModeler2/

#make variables of each repeat masker run
# ln -s /work/soghigian_lab/gen.morinaga/mosquito_hifi_new2/MadSpNov/MadSpNov_canu/ph/MadSpNov_canu_ph.fasta

rm1=rep_mask_r1
rm2=rep_mask_r2
rm3=rep_mask_r3

#set names for output directory, file prefixes and the scaffolded draft assembly that will be masked
mkdir rep_mask_combined
out=rep_mask_combined
name=assembly
asm=assembly.fasta

#paths to several utility scripts from RepeatMasker for combining outputs
rmgff=/work/soghigian_lab/apps/conda/envs/RepeatModeler2/share/RepeatMasker/util/rmOutToGFF3.pl
combine=/work/soghigian_lab/apps/conda/envs/RepeatModeler2/share/RepeatMasker/util/combineRMFiles.pl
build=/work/soghigian_lab/apps/conda/envs/RepeatModeler2/share/RepeatMasker/util/buildSummary.pl

echo combining $rm1 $rm2 $rm3
echo -e output directory is $out '\n'

file1=`ls ./repeat_masker/${rm1}/*.align`
file2=`ls ./repeat_masker/${rm2}/*.align`

mkdir ./${out}/tmp1
perl $combine ${file1%.align} ${file2%.align} ./${out}/tmp1/tmp1

echo -e '\n'
file1=`ls ./${out}/tmp1/*.align`
file2=`ls ./repeat_masker/${rm3}/*.align`

mkdir ./${out}/tmp2
perl $combine ${file1%.align} ${file2%.align} ./${out}/tmp2/tmp2

# echo -e '\n'
# file1=`ls ./${out}/tmp2/*.align`
# file2=`ls ./${rm4}/*.align`
# perl $combine ${file1%.align} ${file2%.align} ./${out}/${name}

echo -e building new table '\n'
cat $asm | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > ./${out}/test.length.tsv

$build -species Diptera -genome ./${out}/*.tsv -useAbsoluteGenomeSize ./${out}/tmp2/*.out > ./${out}/${name}.tabulate

echo -e write combined output to *.gff and *.masked.fasta '\n'
$rmgff ./${out}/tmp2/*.out > ./${out}/${name}.gff


#mask scaffolded draft assembly acording to *.gff produced above
bedtools maskfasta -soft -fi $asm -bed ./${out}/*.gff -fo ./${out}/${name}_masked.fasta

gzip ./${out}/tmp1/*.align
gzip ./${out}/tmp2/*.align

echo use tmp2.align.gz for repeat landscape 