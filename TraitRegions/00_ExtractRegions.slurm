#!/bin/bash

#SBATCH --job-name=subTrait           # Name of your job
#SBATCH --cpus-per-task=1               # Number of cores to reserve
#SBATCH --mem-per-cpu=200M                # Amount of RAM/core to reserve
#SBATCH --time=00:05:00                 # Maximum allocated time
#SBATCH --qos=30min                    # Selected queue to allocate your job
#SBATCH --output=subset.o%j.%a    # Path and name to the file for the STDOUT
#SBATCH --error=subset.e%j.%a     # Path and name to the file for the STDERR
#SBATCH --array=1-7

cd ./8_Cline/TraitRegions

# Load module for dependencies
ml load BCFtools/1.19-GCC-13.2.0


VCF_FILE="/scicore/home/marque0000/GROUP/wu0006/5.2_bcftools/20260106_Chapter1/dataset/Chapter1.NA_Cline.auto.new.bSNPs.2.vcf.gz"
#VCF_FILE="/scicore/home/marque0000/GROUP/wu0006/5.2_bcftools/20260106_Chapter1/dataset/Chapter1.EU_Cline.auto.new.bSNPs.2.vcf.gz"
REGION_FILE="TraitRegions.txt"

# Get region info for this array task
REGION_LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $REGION_FILE)
ID=$(echo $REGION_LINE | awk '{print $1}')
CHR=$(echo $REGION_LINE | awk '{print $2}')
START=$(echo $REGION_LINE | awk '{print $3}')
END=$(echo $REGION_LINE | awk '{print $4}')

OUTPUT="NA.auto.Cline.${ID}.vcf.gz"
#OUTPUT="EU.auto.Cline.${ID}.vcf.gz"

# Extract region

echo $CHR
echo $START
echo $END

bcftools view -r ${CHR}:${START}-${END} -Oz -o $OUTPUT $VCF_FILE
