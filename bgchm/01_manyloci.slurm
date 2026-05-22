#!/bin/bash

#SBATCH --job-name=bgc           # Name of your job
#SBATCH --cpus-per-task=4               # Number of cores to reserve
#SBATCH --mem-per-cpu=2G                # Amount of RAM/core to reserve
#SBATCH --time=1-00:00:00                 # Maximum allocated time
#SBATCH --qos=1day                    # Selected queue to allocate your job
#SBATCH --output=bgc0.o%j.%a    # Path and name to the file for the STDOUT
#SBATCH --error=bgc0.e%j.%a     # Path and name to the file for the STDERR
#SBATCH --array=1-6

ml load R/4.4.1-gfbf-2023b

cd ./8_Cline/

RDATA_FILES=(
  "EU.ChrW.Cline_input.RData"
  "EU.ChrZ.Cline_input.RData"
  "EU.auto.Cline_input.RData"
  "NA.ChrW.Cline_input.RData"
  "NA.ChrZ.Cline_input.RData"
  "NA.auto.Cline_input.RData"
)

RDATA_FILE=${RDATA_FILES[$((SLURM_ARRAY_TASK_ID - 1))]}

echo "Running on file: ${RDATA_FILE}"

Rscript bgchm_Step1_manyloci.R "${RDATA_FILE}"
