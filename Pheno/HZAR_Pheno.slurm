#!/bin/bash

#SBATCH --job-name=hzar           # Name of your job
#SBATCH --cpus-per-task=8               # Number of cores to reserve
#SBATCH --mem-per-cpu=1G                # Amount of RAM/core to reserve
#SBATCH --time=01:00:00                 # Maximum allocated time
#SBATCH --qos=6hours                    # Selected queue to allocate your job
#SBATCH --output=hzarEU.o%j.%a   # Path and name to the file for the STDOUT
#SBATCH --error=hzarEU.e%j.%a     # Path and name to the file for the STDERR
#SBATCH --array=1-4

ml load R/4.4.1-gfbf-2023b

cd ./8_Cline/EU_HZAR/Pheno
Rscript HZAR_Pheno_EU.R $SLURM_ARRAY_TASK_ID

cd ./8_Cline/NA_HZAR/Pheno
Rscript HZAR_Pheno_NA.R $SLURM_ARRAY_TASK_ID
