#!/bin/bash
#SBATCH --job-name=bcf_filter
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=1-00:00:00
#SBATCH --qos=1day
#SBATCH --output=bSNPs.o%j.%a
#SBATCH --error=bSNPs.e%j.%a
#SBATCH --array=1-6   

module load BCFtools/1.19-GCC-13.2.0
bcftoolspath=../99_scripts/bcftools/
export BCFTOOLS_PLUGINS=../99_scripts/bcftools/plugins

cd ./5_bcftools/dataset

vcf=$(sed -n "${SLURM_ARRAY_TASK_ID}p" vcflist.txt)
base=$(basename "$vcf" .vcf.gz)
outfile=${base}.bSNPs.2.vcf.gz

echo "Processing $vcf"

#filtering
if [[ "$vcf" == *ChrW* ]]; then
    echo "ChrW filtering"

    ${bcftoolspath}bcftools filter -Oz -g3 "$vcf" | \
    ${bcftoolspath}bcftools view -Oz \
        -i 'N_ALT=1 & MAC>0 & F_MISSING<0.75 & MAF[0]>=0.05' \
        -o "$outfile"

else
    echo "Autosome / ChrZ filtering"

    ${bcftoolspath}bcftools filter -Oz -g3 "$vcf" | \
    ${bcftoolspath}bcftools view -Oz \
        -i 'N_ALT=1 & MAC>0 & F_MISSING<0.5 & MAF[0]>=0.05' \
        -o "$outfile"
fi

tabix -fp vcf "$outfile"
