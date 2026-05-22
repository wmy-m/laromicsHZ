#!/bin/bash
#SBATCH --job-name=winfst
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100M
#SBATCH --time=00:29:00
#SBATCH --qos=30min
#SBATCH --output=myrun.o%j
#SBATCH --error=myrun.e%j

cd ./5_bcftools/datasets/FST

for j in EU.auto EU.ChrW.forced EU.ChrZ.forced NA.auto NA.ChrW.forced NA.ChrZ.forced; do
  for i in 10000 25000 100000 250000; do
    tail -n +2 /scicore/home/marque0000/GROUP/wu0006/5.2_bcftools/20250916_mpile_all/FST/${j}.wcfst.txt | \
    awk -v w=${i} 'BEGIN{c="";s=1;e=w;a=0;b=0;n=0;OFS="\t"}{if(c==""){c=$1};if(c==$1){if($3!="-nan"){if(s<$2 && e>=$2){a+=$3;b+=$4;n+=1}else{if(b>0){print c,s,e,n,a/b}else{print c,s,e,n,"nan"};for(i=s;i<$2-w;i+=w){s+=w;e+=w;a=0;b=0;n=0};a+=$3;b+=$4;n+=1}}}else{if(b>0){print c,s,e,n,a/b}else{print c,s,e,n,"nan"};c=$1;s=1;e=w;a=0;b=0;n=0;if($3!="-nan"){a+=$3;b+=$4;n+=1}}}END{if(b>0){print c,s,e,n,a/b}else{print c,s,e,n,"nan"}}' \
> ${j}.${i}win.fst.txt
  done
done
