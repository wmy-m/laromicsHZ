#!/bin/bash
#SBATCH --job-name=plink           # Name of your job
#SBATCH --cpus-per-task=1               # Number of cores to reserve
#SBATCH --mem-per-cpu=4G                # Amount of RAM/core to reserve
#SBATCH --time=00:30:00                 # Maximum allocated time
#SBATCH --qos=30min                    # Selected queue to allocate your job
#SBATCH --output=plink.o%j    # Path and name to the file for the STDOUT
#SBATCH --error=plink.e%j     # Path and name to the file for the STDERR
#SBATCH --array=1-20

ml load PLINK/2.00a3.7-gfbf-2023a

cd ./6_GEMMA/EU

phenotypes=("NB" "T10" "T9" "T5" "T4" "OrbitalRingHueCat2" "IrisPigment"\
            "Tar_R" "W10_R" "W9_R" "B10_R" "B9_R" \
            "B8_R" "B7_R" "B6_R" "P10C_R" "P9C_R" "BillLength_R" "BillShape_R" "Bill6_R")

pheno="${phenotypes[$SLURM_ARRAY_TASK_ID - 1]}"

# Autosomes
plink --vcf ../../5_bcftools/datasets/Chapter1.EU_Pheno.New.auto.vcf.gz \
      --double-id --allow-extra-chr --chr-set 95 no-xy no-mt \
      --allow-no-sex --set-missing-var-ids @:# \
      --pheno ./EU.PhenoInput_16042026.txt \ #from PhenoInput.R
      --threads 4 --make-bed \
      --out "EU.GEMMA.auto.$pheno" \
      --pheno-name "$pheno"

# ChrZ
plink --vcf ../../5_bcftools/datasets/Chapter1.EU_Pheno.New.ChrZ.vcf.gz \
      --double-id --allow-extra-chr --chr-set 95 no-xy no-mt \
      --allow-no-sex --set-missing-var-ids @:# \
      --pheno ./EU.PhenoInput_16042026.txt \
      --threads 4 --make-bed \
      --out "EU.GEMMA.ChrZ.$pheno" \
      --pheno-name "$pheno"

# ChrW
plink --vcf ../../5_bcftools/datasets/Chapter1.EU_Pheno.New.ChrW.vcf.gz \
      --double-id --allow-extra-chr --chr-set 95 no-xy no-mt \
      --allow-no-sex --set-missing-var-ids @:# \
      --pheno ./EU.PhenoInput_16042026.txt \
      --threads 4 --make-bed \
      --out "EU.GEMMA.ChrW.$pheno" \
      --pheno-name "$pheno"

echo "Finished $pheno"
