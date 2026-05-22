#!/bin/bash

#SBATCH --job-name=bgc           # Name of your job
#SBATCH --cpus-per-task=1               # Number of cores to reserve
#SBATCH --mem-per-cpu=12G                # Amount of RAM/core to reserve
#SBATCH --time=1-00:00:00                 # Maximum allocated time
#SBATCH --qos=1day                    # Selected queue to allocate your job
#SBATCH --output=bgc2.auto.o%j.%a   # Path and name to the file for the STDOUT
#SBATCH --error=bgc2.auto.e%j.%a     # Path and name to the file for the STDERR
#SBATCH --array=27

ml load R/4.4.1-gfbf-2023b

cd ./8_Cline

# this file can be split so each chromosome class can be run concurrently 
#Rscript --vanilla bgchm_Step2_fitSnps_EU_auto.R $SLURM_ARRAY_TASK_ID
#Rscript --vanilla bgchm_Step2_fitSnps_EU_ChrZ.R $SLURM_ARRAY_TASK_ID
#Rscript --vanilla bgchm_Step2_fitSnps_EU_ChrW.R 

#Rscript --vanilla bgchm_Step2_fitSnps_NA_auto.R $SLURM_ARRAY_TASK_ID
#Rscript --vanilla bgchm_Step2_fitSnps_NA_ChrZ.R $SLURM_ARRAY_TASK_ID
Rscript --vanilla bgchm_Step2_fitSnps_NA_ChrW.R 

