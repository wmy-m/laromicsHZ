#!/bin/bash
#SBATCH --job-name=conChr
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100M
#SBATCH --time=02:00:00
#SBATCH --qos=6hours
#SBATCH --output=conchr.o%j
#SBATCH --error=conchr.e%j

cd ./5_bcftools

# Load module for dependencies
module load BCFtools/1.19-GCC-13.2.0
# Define path to latest bcftools version
bcftoolspath=../99_scripts/bcftools/

# Make combine file
echo "Running concat for autosomes (excluding Z and W)"
${bcftoolspath}bcftools concat -n -o Chapter1.autosomes.bSNPs.new.vcf.gz $(ls *bSNPs.new.vcf.gz | grep -v -E 'chrW|chrZ|MT')
tabix -p vcf Chapter1.autosomes.bSNPs.new.vcf.gz
