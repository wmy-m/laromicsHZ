#!/bin/bash

#SBATCH --job-name=subTrait           # Name of your job
#SBATCH --cpus-per-task=1               # Number of cores to reserve
#SBATCH --mem-per-cpu=200M                # Amount of RAM/core to reserve
#SBATCH --time=00:05:00                 # Maximum allocated time
#SBATCH --qos=30min                    # Selected queue to allocate your job
#SBATCH --output=subset.o%j.%a    # Path and name to the file for the STDOUT
#SBATCH --error=subset.e%j.%a     # Path and name to the file for the STDERR
#SBATCH --array=1-1000

cd ./8_Cline/TraitRegions/backgroundLoci

# Load module for dependencies
ml load BCFtools/1.19-GCC-13.2.0

VCF_FILE="../../../5_bcftools/datasets/Chapter1.EU_Cline.auto.new.bSNPs.2.vcf.gz"
REGION_FILE="background.EU.snps.region.txt"

#VCF_FILE="../../../5_bcftools/datasets/Chapter1.NA_Cline.auto.new.bSNPs.2.vcf.gz"
#REGION_FILE="background.NA.snps.region.txt"

# Get region info for this array task
REGION_LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $REGION_FILE)
CHR=$(echo $REGION_LINE | awk '{print $1}')
START=$(echo $REGION_LINE | awk '{print $2}')
END=$(echo $REGION_LINE | awk '{print $3}')
BED=`echo "$CHR:$START-$END"`
OUTPUT="EU.auto.Cline.${BED}.vcf.gz"

# Extract region

echo $CHR
echo $START
echo $END

bcftools view -r ${BED} -Oz -o $OUTPUT $VCF_FILE
