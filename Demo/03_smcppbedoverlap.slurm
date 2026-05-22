#!/bin/bash
#SBATCH --job-name=smcppbedoverlap
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=25G
#SBATCH --time=00:30:00
#SBATCH --qos=30min
#SBATCH --output=myrun.o%j
#SBATCH --error=myrun.e%j

module load BEDTools/2.31.0-GCC-12.3.0 BCFtools/1.18-GCC-12.3.0
bcftoolspath=./99_scripts/bcftools/
export BCFTOOLS_PLUGINS=./99_scripts/bcftools/plugins

cd ./7_smcpp

indlist=Demo.indlist
totalinds=40
beds=$(for i in $(cat ${indlist}); do ls beds/* | grep $i; done | sed ':a;N;$!ba;s/\n/ /g')
bedtools multiinter -i ${beds} | awk -v "t=${totalinds}" 'BEGIN{OFS="\t"}{if($4/t<0.9){print $1,$2,$3}}' > ${indlist%.indlist}.missing.bed
bedtools sort -i ${indlist%.indlist}.missing.bed > ${indlist%.indlist}.missing.sorted.bed
bedtools merge -i ${indlist%.indlist}.missing.sorted.bed > ${indlist%.indlist}.miss.bed
bgzip ${indlist%.indlist}.miss.bed
tabix -p bed ${indlist%.indlist}.miss.bed.gz
