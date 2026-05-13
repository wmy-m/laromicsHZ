#!/bin/bash
#SBATCH --job-name=bcf_opt
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=600M
#SBATCH --time=1-00:00:00
#SBATCH --qos=1day
#SBATCH --output=myrun.o%j.%a
#SBATCH --error=myrun.e%j.%a
#SBATCH --array=398,399,410,413-417,419-421,578,663

# Load modules
module load BCFtools/1.19-GCC-13.2.0

cd ./5_bcftools

# Define reference genome and bam file list, and output directory
ref=./GCF_964199755.1_bLarMic1.1_genomic_chr.fna
bamlist=./dataset1.list
outdir=./5_bcftools
bcftoolspath=./99_scripts/bcftools/
# Get the parition to call, convert 0-based to 1-based
part=$(sed ''$SLURM_ARRAY_TASK_ID'q;d' bLarMic1.2Mbpwin.bed | awk '{print $1":"$2+1"-"$3}')
# Get list of bam files to call
bams=$(cat $bamlist | sed ':a;N;$!ba;s/\n/ /g')
# Run calling
${bcftoolspath}bcftools mpileup \
-Ou --max-depth 50 -q 30 -Q 20 -a AD -r ${part} -f ${ref} ${bams} | \
${bcftoolspath}bcftools call \
-m -v --ploidy-file bLarMic1.ploidy.MF.txt -S 20260106_SexG.MF.ped \
-G - -v -a GQ -Oz -o bcfoptpart.${SLURM_ARRAY_TASK_ID}.vcf.gz
