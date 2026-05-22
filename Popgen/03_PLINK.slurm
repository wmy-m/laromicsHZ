#!/bin/bash
#SBATCH --job-name=pca
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:30:00
#SBATCH --qos=30min
#SBATCH --output=pca.o%j.%a
#SBATCH --error=pca.e%j.%a
#SBATCH --array=1-4   

ml load PLINK/2.00a3.7-gfbf-2023a  

cd ./5.2_bcftools/datasets/PCA

# Define your input files
vcf_files=(
Chapter1.EU_full.auto.new.bSNPs.2.vcf.gz
Chapter1.NA_Full.auto.new.bSNPs.2.vcf.gz
)

vcf=${vcf_files[$SLURM_ARRAY_TASK_ID-1]}
prefix=$(basename $vcf .vcf.gz)

echo "Processing $vcf ..."

plink --vcf ./$vcf --allow-extra-chr --chr-set 95 no-xy no-mt --make-bed --out ${prefix}
plink --bfile ${prefix} --missing --allow-extra-chr --chr-set 95 no-xy no-mt --out ${prefix}.miss

# run PCA
plink --bfile ${prefix} --allow-extra-chr --chr-set 95 no-xy no-mt --pca --out ${prefix}.pca

# Filter for Admixture
plink --vcf ./$vcf --allow-extra-chr --chr-set 95 no-xy no-mt --geno 0 --mac 2 --bp-space 10000 \
    --recode vcf-iid --out ${prefix}_adm_filter
plink --vcf ${prefix}_adm_filter.vcf --make-bed --chr-set 95 no-xy no-mt --out ${prefix}_adm \
    --set-missing-var-ids @:# --allow-extra-chr

echo "Finished $vcf"
