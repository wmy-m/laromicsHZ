#!/bin/bash
#SBATCH --job-name=bcf_conchr
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=200M
#SBATCH --time=00:30:00
#SBATCH --qos=30min
#SBATCH --output=conchr.o%j.%a
#SBATCH --error=conchr.e%j.%a
#SBATCH --array=36

cd ./5_bcftools

# Load module for dependencies
module load BCFtools/1.19-GCC-13.2.0
# Define path to latest bcftools version
bcftoolspath=../99_scripts/bcftools/

# Get file names / partition numbers for respetive chromosome or scaffolds
autosomes=$(printf "chr%02d " {1..33})
specials="chrW chrZ NW"
all_groups="$autosomes $specials"
c=$(echo $all_groups | cut -d " " -f ${SLURM_ARRAY_TASK_ID})
vcffiles=$(grep -n ${c} bLarMic1.2Mbpwin.bed | cut -f 1 -d ":" | awk '{print "bcfoptpart."$1".vcf.gz"}' | sed ':a;N;$!ba;s/\n/ /g')

# Concatenate parts into file
${bcftoolspath}bcftools concat -n -o 20260106_Chapter1.${c}.vcf.gz ${vcffiles}
# Index VCF file
tabix -p vcf 20260106_Chapter1.${c}.vcf.gz

