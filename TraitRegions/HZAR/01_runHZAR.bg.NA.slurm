#!/bin/bash
#SBATCH --job-name=hzar           # Name of your job
#SBATCH --cpus-per-task=4               # Number of cores to reserve
#SBATCH --mem-per-cpu=15G                # Amount of RAM/core to reserve
#SBATCH --time=00:05:00                 # Maximum allocated time
#SBATCH --qos=30min                    # Selected queue to allocate your job
#SBATCH --output=hzarNA.v3.o%j.%a   # Path and name to the file for the STDOUT
#SBATCH --error=hzarNA.v3.e%j.%a     # Path and name to the file for the STDERR
#SBATCH --array=2-1000

# most can run on 15G mem but some need more
ml load R/4.4.1-gfbf-2023b

cd ./8_Cline/TraitRegions/HZAR

echo "Array Task ID: $SLURM_ARRAY_TASK_ID"

DATA_FILE="../backgroundLoci/NA.HI.background.csv"
META_FILE="../NA.meta.txt"

echo "Processing $DATA_FILE"

# Run R script
Rscript 01_runHZAR.bg.NA.R "$DATA_FILE" "$SLURM_ARRAY_TASK_ID" "$META_FILE"
