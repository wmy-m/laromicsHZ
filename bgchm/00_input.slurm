#!/bin/bash

#SBATCH --job-name=bgc           # Name of your job
#SBATCH --cpus-per-task=4               # Number of cores to reserve
#SBATCH --mem-per-cpu=2G                # Amount of RAM/core to reserve
#SBATCH --time=06:00:00                 # Maximum allocated time
#SBATCH --qos=6hours                    # Selected queue to allocate your job
#SBATCH --output=bgc0.o%j.a    # Path and name to the file for the STDOUT
#SBATCH --error=bgc0.e%j.%a     # Path and name to the file for the STDERR

ml load R/4.4.1-gfbf-2023b

cd ./8_Cline

Rscript bgchm_Step0_input.R

