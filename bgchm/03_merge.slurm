#!/bin/bash

#SBATCH --job-name=bgc           # Name of your job
#SBATCH --cpus-per-task=4               # Number of cores to reserve
#SBATCH --mem-per-cpu=2G                # Amount of RAM/core to reserve
#SBATCH --time=03:00:00                 # Maximum allocated time
#SBATCH --qos=6hours                    # Selected queue to allocate your job
#SBATCH --output=bgc3.o%j    # Path and name to the file for the STDOUT
#SBATCH --error=bgc3.e%j     # Path and name to the file for the STDERR

ml load R/4.4.1-gfbf-2023b

cd ./8_Cline/EU

Rscript bgchm_Step3_Merge_ChrW_EU.R
Rscript bgchm_Step3_Merge_ChrZ_EU.R
Rscript bgchm_Step3_Merge_auto_EU.R

cd ./8_Cline/NA

Rscript bgchm_Step3_Merge_ChrW_NA.R
Rscript bgchm_Step3_Merge_ChrZ_NA.R
Rscript bgchm_Step3_Merge_auto_NA.R
