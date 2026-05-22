#!/bin/bash
#SBATCH --job-name=smcppbamdepth
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=200M
#SBATCH --time=01:00:00
#SBATCH --qos=6hours
#SBATCH --output=depth.o%j
#SBATCH --error=depth.e%j
#SBATCH --array=1-40

module load SAMtools/1.18-GCC-12.3.0 BEDTools/2.31.0-GCC-12.3.0

cd ./7_smcpp

ind=$(sed ''$SLURM_ARRAY_TASK_ID'q;d' Demo.indlist)
bam='../4_MD/'${ind}'.bam'
samtools depth -aa ${bam} | \
awk -v "t=${t}" 'BEGIN{OFS="\t"}{if($3>4){print $1,$2-1,$2}}' | \
bedtools merge > beds/${ind}.DP5.bed
