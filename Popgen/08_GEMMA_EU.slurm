#!/bin/bash
#SBATCH --job-name=gemma           # Name of your job
#SBATCH --cpus-per-task=1               # Number of cores to reserve
#SBATCH --mem-per-cpu=10G                # Amount of RAM/core to reserve
#SBATCH --time=02:00:00                 # Maximum allocated time
#SBATCH --qos=6hours                    # Selected queue to allocate your job
#SBATCH --output=gemma.o%j.%a    # Path and name to the file for the STDOUT
#SBATCH --error=gemma.e%j.%a     # Path and name to the file for the STDERR
#SBATCH --array=1-20


conda activate gwas
cd .6_GEMMA/EU

phenotypes=("NB" "T10" "T9" "T5" "T4" "OrbitalRingHueCat2" "IrisPigment"\
            "Tar_R" "W10_R" "W9_R" "B10_R" "B9_R" \
            "B8_R" "B7_R" "B6_R" "P10C_R" "P9C_R" "BillLength_R" "BillShape_R" "Bill6_R")

# Get the phenotype name based on the SLURM array task ID
pheno="${phenotypes[$SLURM_ARRAY_TASK_ID - 1]}"
chroms=("auto" "ChrZ" "ChrW")

for chrom in "${chroms[@]}"
do
    echo "Running $pheno on $chrom"

    bfile="EU.GEMMA.${chrom}.${pheno}"

    # Step 1: Relatedness matrix
    gemma -bfile "$bfile" \
          -gk 1 \
          -o "${chrom}.${pheno}.relatedness" \
          -miss 0.2 \
          -outdir ./Univariate

    # Step 2: LMM
    gemma -bfile "$bfile" \
          -lmm 4 \
          -o "${chrom}.${pheno}.lmm" \
          -miss 0.2 \
          -outdir ./Univariate \
          -k "./Univariate/${chrom}.${pheno}.relatedness.cXX.txt"

done

echo "Finished phenotype: $pheno"

outfile="Univariate/${pheno}.ALL.lmm.assoc.txt"

echo "Merging results for $pheno..."

# Take header from autosomes
head -n 1 "Univariate/auto.${pheno}.lmm.assoc.txt" > "$outfile"

# Append all data 
tail -n +2 -q \
    "Univariate/auto.${pheno}.lmm.assoc.txt" \
    "Univariate/ChrZ.${pheno}.lmm.assoc.txt" \
    "Univariate/ChrW.${pheno}.lmm.assoc.txt" \
    >> "$outfile"

echo "Merged file created: $outfile"

awk '{print $1, $3, $13, $14, $15}' "$outfile" | gzip > "${outfile}.reduced.gz"

echo "Reduced file created: ${outfile}.reduced.gz"
