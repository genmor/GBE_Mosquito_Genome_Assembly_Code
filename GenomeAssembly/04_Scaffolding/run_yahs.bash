#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --time=24:00:00
#SBATCH --mem=45G

####### Set environment variables ###############
module load lib/zlib/1.2.11
module load gcc/10.2.0
module load perl/5.32.0
module load java/openjdk-23.0.1
PATH=/work/soghigian_lab/apps/bwa:$PATH
PATH=/work/soghigian_lab/apps/bedtools2/bin:$PATH
PATH=/work/soghigian_lab/apps/samtools-1.15/bin:$PATH
PATH=/work/soghigian_lab/apps/yahs:$PATH
PATH=/work/soghigian_lab/apps/fastp:$PATH
juicer_tools="java -Xmx45G -jar /work/soghigian_lab/apps/juicer/juicer_tools_1.22.01.jar pre"
pretext_map=/work/soghigian_lab/apps/conda/envs/pretextmap/bin/PretextMap
pretext_snapshot=/work/soghigian_lab/apps/conda/envs/pretextsnapshot/bin/PretextSnapshot

cd /scratch/${SLURM_JOB_ID}
echo check /scratch/${SLURM_JOB_ID}

#copy assembly and hic reads to scratch re-name reads to something more sensible
ln -s path/to/hic/reads1
ln -s path/to/hic/reads2
ln -s path/to/assembly

#define variables for file processing
asm=Bf05_canu_ph_insp.fa
og_reads1=reads1
og_reads2=reads2
fastp_out=fastp_out

#misc. variables
qual=10

#paths to programs and scripts
filter=/work/soghigian_lab/apps/arima/mapping_pipeline/filter_five_end.pl
combiner=/work/soghigian_lab/apps/arima/mapping_pipeline/two_read_bam_combiner.pl
picard=/work/soghigian_lab/apps/picard.jar
statsQ=/work/soghigian_lab/apps/arima/mapping_pipeline/get_stats.pl

#create directories for files
mkdir ./1_bam
mkdir ./2_filtered_bam
mkdir ./3_paired_bam
mkdir ./tmp
mkdir ./4_de_duped


####trim the first 5 bases from each read file using fastp####
echo -e '\n\ntrim the 1st 5 bases from each read file\n\n'
fastp -i $og_reads1 -I $og_reads2 --out1 ${fastp_out}_r1.fq.gz --out2 ${fastp_out}_r2.fq.gz --html ${fastp_out} --detect_adapter_for_pe -f 5 -w 12 -w 12

####most of this is from the ARIMA alignment/postprocessing pipeline
echo -e '\n\ncreate index\n\n'
bwa index -a bwtsw ${asm}
samtools faidx ${asm}


####comment out depending on bz2 or gz
echo -e '\n\nalign reads\n\n'
r1=${fastp_out}_r1.fq.gz
r2=${fastp_out}_r2.fq.gz
# bwa mem -t25 ${asm} <(bzip -dc ${r1}) | samtools view -@ 25 -Sb - > ./1_bam/${asm%.fa*}_r1.bam
# bwa mem -t25 ${asm} <(bzip -dc ${r2}) | samtools view -@ 25 -Sb - > ./1_bam/${asm%.fa*}_r2.bam

bwa mem -t25 ${asm} ${r1} | samtools view -@ 25 -Sb - > ./1_bam/${asm%.fa*}_r1.bam
bwa mem -t25 ${asm} ${r2} | samtools view -@ 25 -Sb - > ./1_bam/${asm%.fa*}_r2.bam

echo -e '\n\nfilter 5 end\n\n'
samtools view -h ./1_bam/${asm%.fa*}_r1.bam | perl ${filter} | samtools view -Sb - > ./2_filtered_bam/${asm%.fa*}_r1.bam
samtools view -h ./1_bam/${asm%.fa*}_r2.bam | perl ${filter} | samtools view -Sb - > ./2_filtered_bam/${asm%.fa*}_r2.bam

echo -e '\n\ncreate paired end reads and filter\n\n'
perl ${combiner} ./2_filtered_bam/${asm%.fa*}_r1.bam ./2_filtered_bam/${asm%.fa*}_r2.bam samtools ${qual} | samtools view -bS -t ${asm}.fai - | samtools sort -@ 25 -o ./tmp/${asm%.fa*}_tmp.bam -

echo -e '\n\nadd read group\n\n'
java -Xmx4G -Djava.io.tmpdir=temp/ -jar ${picard} AddOrReplaceReadGroups INPUT=./tmp/${asm%.fa*}_tmp.bam OUTPUT=./3_paired_bam/${asm%.fa*}_pe.bam ID=${asm%.fa*} LB=${asm%.fa*} SM=${asm%.fa*} PL=ILLUMINA PU=none

echo -e '\n\n\mark duplicates\n\n'
java -Xmx30G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ -jar ${picard} MarkDuplicates INPUT=./3_paired_bam/${asm%.fa*}_pe.bam OUTPUT=./4_de_duped/${asm%.fa*}_dedupe.bam METRICS_FILE=./4_de_duped/metrics_${asm%.fa*}_dedupe.txt TMP_DIR=./tmp ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

samtools index ./4_de_duped/${asm%.fa*}_dedupe.bam

perl ${statsQ} ./4_de_duped/${asm%.fa*}_dedupe.bam > ./4_de_duped/${asm%.fa*}_dedupe.stats

####run yahs
#if enzymes need to be specified use -e GATC,GANTC,CTNAG,TTAA
echo -e '\n\nrun yahs\n\n'
mkdir ${asm%.fa*}_yahs
yahs ${asm} ./4_de_duped/${asm%.fa*}_dedupe.bam -o ./${asm%.fa*}_yahs/${asm%.fa*}_yahs

####run juicer tools to output hic files for viz using juicebox
echo -e '\n\ngenerate .assembly and .hic files for editing in juicebox\n\n'
mkdir hic
cd hic
cp ../${asm%.fa*}_yahs/${asm%.fa*}_yahs.bin .
cp ../${asm%.fa*}_yahs/${asm%.fa*}_yahs_scaffolds_final.agp .
cp ../${asm%.fa*}_yahs/${asm%.fa*}_yahs_scaffolds_final.fa .
cp ../${asm}.fai .
out=${asm%.fa*}_yahs
bin=${asm%.fa*}_yahs.bin
agp=${asm%.fa*}_yahs_scaffolds_final.agp
scafs=${asm%.fa*}_yahs_scaffolds_final.fa

echo -e '\n\nGenerate input file for juicer_tools non-assembly mode and pretext map (1/4)\n\n'
(juicer pre ${bin} ${agp} ${asm}.fai 2> tmp_juicer_pre.log | LC_ALL=C sort -k2,2d -k6,6d -T . --parallel=25 -S45G | awk 'NF' > alignments_sorted.txt.part) && (mv alignments_sorted.txt.part alignments_sorted.txt)

echo -e '\n\nCreate .chrom.sizes file (2/4)\n\n'
cat tmp_juicer_pre.log | grep "PRE_C_SIZE" | cut -d' ' -f2- >${out}_scaffolds_final.chrom.sizes

echo -e '\n\nCreate .hic file (3/4)\n\n'
(${juicer_tools} alignments_sorted.txt ${out}.hic.part ${out}_scaffolds_final.chrom.sizes) && (mv ${out}.hic.part ${out}.hic)

echo -e '\n\nUse pretextmap and pretext snapshot to create contact maps (4/4)\n\n'
(awk 'BEGIN{print "## pairs format v1.0"} {print "#chromsize:\t"$1"\t"$2} END {print "#columns:\treadID\tchr1\tpos1\tchr2\tpos2\tstrand1\tstrand2"}' ${out}_scaffolds_final.chrom.sizes; awk '{print ".\t"$2"\t"$3"\t"$6"\t"$7"\t.\t."}' alignments_sorted.txt) | ${pretext_map} -o ${out}.pretext

${pretext_snapshot} -m ${out}.pretext --sequences "=full" -o .

echo -e '\n\nCreate input files for juicer_tools assembly mode (1/3)\n\n'
juicer pre -a -o ${out}_JBAT ${bin} ${agp} ${asm}.fai 2>tmp_juicer_pre_JBAT.log

echo -e '\n\nCreate .chrom.sizes files (2/3)\n\n'
cat tmp_juicer_pre_JBAT.log | grep "PRE_C_SIZE" | cut -d' ' -f2- >${out}_JBAT.chrom.sizes

echo -e '\n\nCreate .hic file (3/3)\n\n'
(${juicer_tools} ${out}_JBAT.txt ${out}_JBAT.hic.part ${out}_JBAT.chrom.sizes) && (mv ${out}_JBAT.hic.part ${out}_JBAT.hic)
#after manual editing in Juicebox output final assembly using the following
#contigs.fa is ${asm}; liftover.agp is created in the previous steps JBAT.review.assembly is created by Juicebox; -o is the output name for the final assembly
#juicer post -o out_JBAT out_JBAT.review.assembly out_JBAT.liftover.agp contigs.fa