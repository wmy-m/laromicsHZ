#!/bin/bash
#SBATCH --job-name=gemma           # Name of your job
#SBATCH --cpus-per-task=1               # Number of cores to reserve
#SBATCH --mem-per-cpu=16G                # Amount of RAM/core to reserve
#SBATCH --time=12:00:00                 # Maximum allocated time
#SBATCH --qos=1day                    # Selected queue to allocate your job
#SBATCH --output=gemma.o%j.%a    # Path and name to the file for the STDOUT
#SBATCH --error=gemma.e%j.%a     # Path and name to the file for the STDERR
#SBATCH --array=13,14,15,16,17,18,19,20,23,25

conda activate gwas
cd ./6_GEMMA/NA

phenotypes=("BillHueAngle1" "BillHueAngle2" "MiddleToe_R" "OrbitalRingHueAngle" "IrisHueAngle" \
            "IrisPigHueAngle" "TarsusHueAngle" "TailLength_R" "Tarsus_R" "IrisPig" "Weight_R" "WingLength_R" \
            "P10_R" "P9_R" "P8_R" "P7_R" "P6_R" "P5_R" "BillLength_R" "BillShape_R" "Bill1_R" "Bill2_R" "Bill4_R" "Bill3_R" \
            "Bill5_R" "Mantle" "Tip" "OrbitalRingHueCat" "OrbitalRingHueCat2")

# Get the phenotype name based on the SLURM array task ID
pheno="${phenotypes[$SLURM_ARRAY_TASK_ID - 1]}"
chroms=("auto" "ChrZ" "ChrW")

for chrom in "${chroms[@]}"
do
    echo "Running $pheno on $chrom"

    bfile="NA.GEMMA.${chrom}.${pheno}"

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
