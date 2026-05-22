#!/bin/bash
#SBATCH --job-name=smcppestd5
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=25G
#SBATCH --time=01:00:00
#SBATCH --qos=6hours
#SBATCH --output=est.o%j%a
#SBATCH --error=est.e%j%a
#SBATCH --array=2-4
	
# Load modules
module load SMC++/1.15.4-foss-2022a
cd ./7_smcpp

# Choose comparison
comp=$(echo LarArg LarCac LarGla LarOcc | cut -f $SLURM_ARRAY_TASK_ID -d " ")
# Run estimation
smc++ estimate -o run${comp}d10 7.90e-09 infiles/${comp}.chr[0-9]*.d[1-10].smc.gz
