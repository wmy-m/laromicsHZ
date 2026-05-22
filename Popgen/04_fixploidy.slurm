#!/bin/bash
#SBATCH --job-name=fixploidy
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:30:00
#SBATCH --qos=30min
#SBATCH --output=fixploidy.o%j.%a
#SBATCH --error=fixploidy.e%j.%a
#SBATCH --array=1-32

module load BCFtools/1.19-GCC-13.2.0
bcftoolspath=./99_scripts/bcftools/
export BCFTOOLS_PLUGINS=./99_scripts/bcftools/plugins

cd /5_bcftools/dataset

    vcf=$(sed -n "${SLURM_ARRAY_TASK_ID}p" fixploidy_list.txt)

    base=$(basename "$vcf" .vcf.gz)
    outfile=${base}.diploid.vcf.gz

    echo "Fixing ploidy for: $vcf"

    ${bcftoolspath}bcftools +fixploidy "$vcf" -- -f 2 | \
       bgzip > "$outfile"

    tabix -fp vcf "$outfile"

    echo "Done: $outfile"
