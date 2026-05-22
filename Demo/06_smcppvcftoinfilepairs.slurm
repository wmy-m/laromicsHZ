#!/bin/bash
#SBATCH --job-name=smcppvcftoinfilepairs
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

# Get population strings
largla=$(cat Demo.indlist | grep Gla | sed ':a;N;$!ba;s/\n/,/g')
larocc=$(cat Demo.indlist | grep Occ | sed ':a;N;$!ba;s/\n/,/g')
larcac=$(cat Demo.indlist | grep Cac | sed ':a;N;$!ba;s/\n/,/g')
lararg=$(cat Demo.indlist | grep Arg | sed ':a;N;$!ba;s/\n/,/g')

# Convert chromosome
smc++ vcf2smc Demo.filt.recoded.vcf.gz infiles/LarGla.LarOcc.${SLURM_ARRAY_TASK_ID}.smc.gz ${chr} LarGla:${largla} LarOcc:${larocc} -m Demo.miss.fixed.bed.gz
smc++ vcf2smc Demo.filt.recoded.vcf.gz infiles/LarOcc.LarGla.${SLURM_ARRAY_TASK_ID}.smc.gz ${chr} LarOcc:${larocc} LarGla:${largla} -m Demo.miss.fixed.bed.gz
smc++ vcf2smc Demo.filt.recoded.vcf.gz infiles/LarCac.LarArg.${SLURM_ARRAY_TASK_ID}.smc.gz ${chr} LarCac:${larcac} LarArg:${lararg} -m Demo.miss.fixed.bed.gz
smc++ vcf2smc Demo.filt.recoded.vcf.gz infiles/LarArg.LarCac.${SLURM_ARRAY_TASK_ID}.smc.gz ${chr} LarArg:${lararg} LarCac:${larcac} -m Demo.miss.fixed.bed.gz
