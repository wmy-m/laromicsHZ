#!/bin/bash
#SBATCH --job-name=triangular
#SBATCH --cpus-per-task=1
#SBATCH --mem=160G
#SBATCH --time=01:00:00
#SBATCH --qos=6hours
#SBATCH --array=1-6
#SBATCH --output=triangular.o%j.%a
#SBATCH --error=triangular.e%j.%a

ml load R/4.4.1-gfbf-2023b

cd ./5_bcftools/datasets/TriangulaR

JOBS=(
"../Chapter1.NA_Full.auto.new.bSNPs.2.vcf.gz ./NApopmap.txt GLAUC OCCI 1"
"../Chapter1.NA_Full.auto.new.bSNPs.2.vcf.gz ./NApopmap.txt GLAUC OCCI 0.8"
"../Chapter1.NA_Full.auto.new.bSNPs.2.vcf.gz ./NApopmap.txt GLAUC OCCI 0.5"
"../Chapter1.EU_full.auto.new.bSNPs.2.vcf.gz ./EUpopmap.txt CACH ARGE 1"
"../Chapter1.EU_full.auto.new.bSNPs.2.vcf.gz ./EUpopmap.txt CACH ARGE 0.8"
"../Chapter1.EU_full.auto.new.bSNPs.2.vcf.gz ./EUpopmap.txt CACH ARGE 0.5"
)

CMD=${JOBS[$SLURM_ARRAY_TASK_ID-1]}

echo "SLURM ARRAY TASK: $SLURM_ARRAY_TASK_ID"
echo "Rscript TriangulaR.R $CMD"

Rscript 01_TriangulaR.R $CMD
