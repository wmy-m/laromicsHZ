#!/bin/bash
#SBATCH --job-name=bSNPs
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=1-00:00:00
#SBATCH --qos=1day
#SBATCH --output=bSNP.o%j.%a
#SBATCH --error=bSNP.e%j.%a
#SBATCH --array=1-36

cd ./5_bcftools

# Load module for dependencies
module load BCFtools/1.19-GCC-13.2.0
# Define path to latest bcftools version
bcftoolspath=./99_scripts/bcftools/
# Define directory to bcftools plugins
export BCFTOOLS_PLUGINS=./99_scripts/bcftools/plugins

# Choose chromosome to work on
autosomes=$(printf "chr%02d " {1..33})
specials="chrW chrZ NW"
all_groups="$autosomes $specials"
c=$(echo $all_groups | cut -d " " -f ${SLURM_ARRAY_TASK_ID})

# Get depth boundary for autosomes (incl. scaffolds), Z chromosome or mtDNA
maxdp=10277   # default autosomes
mindp=5318   # default autosomes

if [ "$c" = "chrZ" ]; then
    maxdp=8067   # Z chromosome
    mindp=3456
elif [ "$c" = "chrW" ]; then
    maxdp=6370   # W chromosome 
    mindp=0   
fi

echo "Running on $c with mindp=$mindp and maxdp=$maxdp"

# Start filtering
# 1) First remove SNPs close to indels
# 2) Remove indels, QUAL<20 and based on INFO stats
# 3) Set genotypes to missing with <4 reads or GQ<20
${bcftoolspath}bcftools filter -Oz -g3 20260106_Chapter1.${c}.vcf.gz | \
${bcftoolspath}bcftools view -Oz -f .,PASS -i 'QUAL>=20'\
' & TYPE!~"indel" & INFO/DP<'${maxdp}' & INFO/DP>'${mindp}\
' & ABS(INFO/RPBZ)<1.96 & ABS(INFO/MQBZ)<1.96 & ABS(INFO/BQBZ)<1.96'\
' & ABS(INFO/MQSBZ)<1.96 & ABS(INFO/SCBZ)<1.96' | \
${bcftoolspath}bcftools plugin setGT -Oz -- -tq -n './.' -i 'SUM(FMT/AD)<4 | FMT/GQ<20' | \
${bcftoolspath}bcftools plugin fill-tags -Oz \
> Chapter1.${c}.bSNPs.new.vcf.gz
# Index filtered VCF file
tabix -fp vcf Chapter1.${c}.bSNPs.new.vcf.gz

