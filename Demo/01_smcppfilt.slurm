#!/bin/bash
#SBATCH --job-name=smcppfilt
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=200M
#SBATCH --time=01:00:00
#SBATCH --qos=6hours
#SBATCH --output=myrun.o%j
#SBATCH --error=myrun.e%j

# Load modules
module load BCFtools/1.19-GCC-13.2.0
bcftoolspath=./99_scripts/bcftools/
export BCFTOOLS_PLUGINS=./99_scripts/bcftools/plugins

cd ./7_smcpp

# Run subset on autosomes and filter command
bcftools view -S Demo.indlist -i 'TYPE!~"indel" & N_ALT<=1' ../5_bcftools/datasets/Chapter1.autosomes.bSNPs.new.vcf.gz | \
bcftools view -i 'F_MISSING<0.1 & MAC>0' | \
bgzip -c > Demo.filt.vcf.gz
tabix -fp vcf Demo.filt.vcf.gz
