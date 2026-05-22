#!/bin/bash
#SBATCH --job-name=smcppsplit
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=25G
#SBATCH --time=00:30:00
#SBATCH --qos=30min
#SBATCH --output=myrun.o%j
#SBATCH --error=myrun.e%j
#SBATCH --array=2

module load SMC++/1.15.4-foss-2022a
cd ./7_smcpp

# Choose comparison
comp=$(echo LarArg.LarCac LarGla.LarOcc | cut -f $SLURM_ARRAY_TASK_ID -d " ")
# Run estimation
smc++ split -o run${comp}002 run${comp%.*}d10/model.final.json run${comp#*.}d10/model.final.json infiles/${comp}.[0-9]*.smc.gz infiles/${comp#*.}.${comp%.*}.[0-9]*.smc.gz infiles/${comp%.*}.chr[0-9]*.d[1-10].smc.gz infiles/${comp#*.}.chr[0-9]*.d[1-10].smc.gz
