#!/bin/bash
#SBATCH --job-name=wcfst
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=200M
#SBATCH --time=00:30:00
#SBATCH --qos=30min
#SBATCH --output=wcfst.o%j
#SBATCH --error=wcfst.e%j

cd ./5_bcftools/dataset/FST

python3 ./99_scripts/vcf2wcfst.py \
    -i ../Chapter1.EU_FST.auto.new.bSNPs.2.vcf.gz \
    -o EU.auto \
    -p Arg.txt -q Cac.txt

python3 ./99_scripts/vcf2wcfst.py \
    -i ../Chapter1.EU_FST.ChrZ.bSNPs.2.diploid.vcf.gz \
    -o EU.ChrZ.forced \
    -p Arg.txt -q Cac.txt

python3 ./99_scripts/vcf2wcfst.py \
    -i ../Chapter1.EU_FST.ChrW.bSNPs.2.diploid.vcf.gz \
    -o EU.ChrW.forced \
    -p ArgW.txt -q CacW.txt

python3 ./99_scripts/vcf2wcfst.py \
    -i ../Chapter1.NA_FST.auto.new.bSNPs.2.vcf.gz \
    -o NA.auto \
    -p Occ.txt -q Gla.txt

python3 ./99_scripts/vcf2wcfst.py \
    -i ../Chapter1.NA_FST.ChrZ.bSNPs.2.diploid.vcf.gz \
    -o NA.ChrZ.forced \
    -p Occ.txt -q Gla.txt

python3 ./99_scripts/vcf2wcfst.py \
    -i ../Chapter1.NA_FST.ChrW.bSNPs.2.diploid.vcf.gz \
    -o NA.ChrW.forced \
    -p OccW.txt -q GlaW.txt
