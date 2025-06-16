#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=12:00:00
#SBATCH --mem=80G

####### Set environment variables ###############
module load gcc/10.2.0
module load bioconda/conda3
module load lib/boost/1.70.0-mpi
module load lib/zlib/1.2.11
export PATH=/work/soghigian_lab/apps/metaeuk/bin/:$PATH
export PATH=/work/soghigian_lab/apps/R-4.1.2/bin:$PATH
export PATH="/work/soghigian_lab/apps/Augustus/bin:$PATH"
export PATH="/work/soghigian_lab/apps/Augustus/scripts:$PATH"
export AUGUSTUS_CONFIG_PATH="/work/soghigian_lab/apps/Augustus/config/"
export PATH="/work/soghigian_lab/apps/ncbi-blast-2.12.0+/bin:$PATH"
export PATH="/work/soghigian_lab/apps/hmmer-3.3.2/bin:$PATH"


##activate conda environment first!
#conda activate /work/soghigian_lab/apps/conda/envs/busco



# cd /scratch/${SLURM_JOB_ID}
# echo check /scratch/${SLURM_JOB_ID}


###for just one assembly
# busco -m genome -i  /work/soghigian_lab/gen.morinaga/mosquito_hifi_new2/Nmex03/Nmex03_synteny/Nmex03_samba_ragtag_test/ragtag_output/Nmex03_insp_samba_ragtag.fasta -o canu_scaf -l /work/soghigian_lab/apps/busco_lineages/diptera_odb10 -c 20 --offline


###for mulitple assemblies
mkdir Bf05_flye_asm
asm_dir=./Bf05_flye_asm
ln -s /work/soghigian_lab/gen.morinaga/mosquito_hifi_new2/Bf05/flye/Bf05_flye.fasta $asm_dir
ln -s /work/soghigian_lab/gen.morinaga/mosquito_hifi_new2/Bf05/flye/pd/Bf05_flye_pd.fasta $asm_dir
ln -s /work/soghigian_lab/gen.morinaga/mosquito_hifi_new2/Bf05/flye/ph/Bf05_flye_ph.fasta $asm_dir

###set output names
spp=Bf05
asm=flye
#run busco in batch mode (autodetects)
busco -m genome -i $asm_dir -o ${spp}_${asm} -l /work/soghigian_lab/apps/busco_lineages/diptera_odb10 -c 20 --offline