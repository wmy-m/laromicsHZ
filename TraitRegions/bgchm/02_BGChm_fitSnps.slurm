#!/bin/bash

#SBATCH --job-name=bgc           # Name of your job
#SBATCH --cpus-per-task=1               # Number of cores to reserve
#SBATCH --mem-per-cpu=12G                # Amount of RAM/core to reserve
#SBATCH --time=1-00:00:00                 # Maximum allocated time
#SBATCH --qos=1day                   # Selected queue to allocate your job
#SBATCH --output=bgc2.auto.o%j.%a   # Path and name to the file for the STDOUT
#SBATCH --error=bgc2.auto.e%j.%a     # Path and name to the file for the STDERR
#SBATCH --array=9

ml load R/4.4.1-gfbf-2023b

cd ./8_Cline/TraitRegions/BGChm/manyloci


#Rscript --vanilla 02_BGChm_fitSnps.R $SLURM_ARRAY_TASK_ID

