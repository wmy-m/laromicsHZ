#!/bin/bash
#SBATCH --job-name=w_filtSNPs
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=01:00:00
#SBATCH --qos=6hours
#SBATCH --output=myrun.o%j
#SBATCH --error=myrun.e%j

cd ./5_bcftools

# Load module for dependencies
module load BCFtools/1.19-GCC-13.2.0
# Define path to latest bcftools version
bcftoolspath=../99_scripts/bcftools/
# Define directory to bcftools plugins
export BCFTOOLS_PLUGINS=../99_scripts/bcftools/plugins


# Removing males
${bcftoolspath}bcftools view -S ^../males.txt --force-samples \
Chapter1.chrW.bSNPs.new.vcf.gz \
-o  Chapter1.chrW.woMales.bSNPs.new.vcf.gz
tabix -fp vcf Chapter1.chrW.woMales.bSNPs.new.vcf.gz
