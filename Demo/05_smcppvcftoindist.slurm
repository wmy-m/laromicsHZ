#!/bin/bash
#SBATCH --job-name=smcppvcftoindist
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=200M
#SBATCH --time=00:30:00
#SBATCH --qos=30min
#SBATCH --output=myrun.o%j
#SBATCH --error=myrun.e%j
#SBATCH --array=1-33

# Load modules
module load SMC++/1.15.4-foss-2022a

cd ./7_smcpp

# Choose chromsome to work on
chroms=({1..33})
chr=${chroms[$((SLURM_ARRAY_TASK_ID-1))]}

# Define dataset and species to work on
# do this for all species
dat=Demo
spe=LarGla #CHANGE THIS

# Get population strings
pop=$(cat ${dat}.indlist | grep ${spe} | sed ':a;N;$!ba;s/\n/,/g')
# Loop over five distinguished individuals
for i in {1..5}; do
	ind=$(cat ${dat}.indlist | grep ${spe} | sed ''${i}'q;d')

	# Convert chromosome
	smc++ vcf2smc -d ${ind} ${ind} ${dat}.filt.recoded.vcf.gz \
	infiles/${spe}.chr${SLURM_ARRAY_TASK_ID}.d${i}.smc.gz \
	${chr} ${spe}:${pop} -m ${dat}.miss.fixed.bed.gz
done
