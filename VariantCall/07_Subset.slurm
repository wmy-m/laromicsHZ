#!/bin/bash
#SBATCH --job-name=bcf_subset
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=06:00:00
#SBATCH --qos=6hours
#SBATCH --output=subset.o%j.%a
#SBATCH --error=subset.e%j.%a
#SBATCH --array=1-2  # number of subset files

# Load module for dependencies
module load BCFtools/1.19-GCC-13.2.0
# Define path to latest bcftools version
bcftoolspath=../99_scripts/bcftools/

cd ./5_bcftools/datasets

# Path to subset files
subset_files=(EU_Cline.txt NA_Cline.txt EU_Pheno.txt NA_Pheno.txt NA_Full.txt EU_full.txt EU_allo.txt NA_allo.txt)
subset_file=${subset_files[$SLURM_ARRAY_TASK_ID-1]}  
subset_name=$(basename $subset_file .txt)

#ZChr
outfile_z=Chapter1.$subset_name.ChrZ.vcf.gz
${bcftoolspath}bcftools view -Oz -S $subset_file ../Chapter1.chrZ.bSNPs.new.vcf.gz -o $outfile_z
tabix -p vcf $outfile_z

#Auto
outfile_auto=Chapter1.$subset_name.auto.new.vcf.gz
${bcftoolspath}bcftools view -Oz -S $subset_file ../Chapter1.autosomes.bSNPs.new.vcf.gz -o $outfile_auto
tabix -p vcf $outfile_auto

#WChr
outfile_w=Chapter1.$subset_name.ChrW.vcf.gz
${bcftoolspath}bcftools view -Oz -S $subset_file --force-samples \
../Chapter1.chrW.woMales.bSNPs.new.vcf.gz -o $outfile_w
tabix -p vcf $outfile_w
