#!/bin/bash
#SBATCH --job-name=plink           # Name of your job
#SBATCH --cpus-per-task=1               # Number of cores to reserve
#SBATCH --mem-per-cpu=4G                # Amount of RAM/core to reserve
#SBATCH --time=00:30:00                 # Maximum allocated time
#SBATCH --qos=30min                    # Selected queue to allocate your job
#SBATCH --output=plink.o%j    # Path and name to the file for the STDOUT
#SBATCH --error=plink.e%j     # Path and name to the file for the STDERR
#SBATCH --array=1-29

ml load PLINK/2.00a3.7-gfbf-2023a

cd ./6_GEMMA/NA

phenotypes=("BillHueAngle1" "BillHueAngle2" "MiddleToe_R" "OrbitalRingHueAngle" "IrisHueAngle" \
            "IrisPigHueAngle" "TarsusHueAngle" "TailLength_R" "Tarsus_R" "IrisPig" "Weight_R" "WingLength_R" \
            "P10_R" "P9_R" "P8_R" "P7_R" "P6_R" "P5_R" "BillLength_R" "BillShape_R" "Bill1_R" "Bill2_R" "Bill4_R" "Bill3_R" \
            "Bill5_R" "Mantle" "Tip" "OrbitalRingHueCat" "OrbitalRingHueCat2")

pheno="${phenotypes[$SLURM_ARRAY_TASK_ID - 1]}"

echo "Running phenotype: $pheno"

# Autosomes
plink --vcf ../../5_bcftools/datasets/Chapter1.NA_Pheno.New.auto.vcf.gz \
      --double-id --allow-extra-chr --chr-set 95 no-xy no-mt \
      --allow-no-sex --set-missing-var-ids @:# \
      --pheno ./NA.PhenoInput_16042026.txt \
      --threads 4 --make-bed \
      --out "NA.GEMMA.auto.$pheno" \
      --pheno-name "$pheno"

# ChrZ
plink --vcf ../../5_bcftools/datasets/Chapter1.NA_Pheno.New.ChrZ.vcf.gz \
      --double-id --allow-extra-chr --chr-set 95 no-xy no-mt \
      --allow-no-sex --set-missing-var-ids @:# \
      --pheno ./NA.PhenoInput_16042026.txt \
      --threads 4 --make-bed \
      --out "NA.GEMMA.ChrZ.$pheno" \
      --pheno-name "$pheno"

# ChrW
plink --vcf ../../5_bcftools/datasets/Chapter1.NA_Pheno.New.ChrW.vcf.gz \
      --double-id --allow-extra-chr --chr-set 95 no-xy no-mt \
      --allow-no-sex --set-missing-var-ids @:# \
      --pheno ./NA.PhenoInput_16042026.txt \
      --threads 4 --make-bed \
      --out "NA.GEMMA.ChrW.$pheno" \
      --pheno-name "$pheno"

echo "Finished $pheno"
