#!/bin/bash
#SBATCH --job-name=bcf_dp
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=500M
#SBATCH --time=02:00:00
#SBATCH --qos=6hours
#SBATCH --output=bcfDepth.o%j
#SBATCH --error=bcfDepth.e%j

# Load module for dependencies
module load BCFtools/1.19-GCC-13.2.0


cd ./5_bcftools/
# Define path to latest bcftools version
bcftoolspath=./99_scripts/bcftools/
# Extract DP values from autosomes
${bcftoolspath}bcftools query -f '%INFO/DP' 20260106_Chapter1.chr{01..33}.vcf.gz 20260106_Chapter1.NW.vcf.gz > tmp.dp.autosomes
# Extract DP values from Z chromosome
${bcftoolspath}bcftools query -f '%INFO/DP' 20260106_Chapter1.chrZ.vcf.gz > tmp.dp.Zchrom
# Extract DP values from W chromosome
${bcftoolspath}bcftools query -f '%INFO/DP' 20260106_Chapter1.chrW.vcf.gz > tmp.dp.Wchrom


module load R/4.4.2-gfbf-2024a
Rscript -e 'a=read.table("tmp.dp.autosomes");b=summary(a$V1);print(b);'\
'pdf("tmp.dp.autosomes.histogram.pdf",width=10,height=5);'\
'hist(as.numeric(a$V1),breaks=1000,main="autosomes",xlab="DP");'\
'abline(v=c(2*b[4],b[4]+c(1,-1)*1.5*(b[5]-b[2])),col=c(2,4,3));'\
'legend("topright",c("2*mean","mean+1.5*interquartilerange","mean-1.5*interquartilerange"),lwd=1,col=c(2,4,3),bty="n");dev.off();'
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#      1    7735    8615    7798    9388   48885 
# MaxDP = 7798+(1.5(9388-7735))=10277
# MinDP = 7798-(1.5(9388-7735))=5318

Rscript -e 'a=read.table("tmp.dp.Zchrom");b=summary(a$V1);print(b);'\
'pdf("tmp.dp.Zchrom.histogram.pdf",width=10,height=5);'\
'hist(as.numeric(a$V1),breaks=1000,main="Zchrom",xlab="DP");'\
'abline(v=c(2*b[4],b[4]+c(1,-1)*1.5*(b[5]-b[2])),col=c(2,4,3));'\
'legend("topright",c("2*mean","mean+1.5*interquartilerange","mean-1.5*interquartilerange"),lwd=1,col=c(2,4,3),bty="n");dev.off();'
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#      1    5591    6432    5762    7128   47814 
# MaxDP = 5762+(1.5(7128-5591))=8067
# MinDP = 5762-(1.5(7128-5591))=3456      

Rscript -e 'a=read.table("tmp.dp.Wchrom");b=summary(a$V1);print(b);'\
'pdf("tmp.dp.Wchrom.histogram.pdf",width=10,height=5);'\
'hist(as.numeric(a$V1),breaks=1000,main="Wchrom",xlab="DP");'\
'abline(v=c(2*b[4],b[4]+c(1,-1)*1.5*(b[5]-b[2])),col=c(2,4,3));'\
'legend("topright",c("2*mean","mean+1.5*interquartilerange","mean-1.5*interquartilerange"),lwd=1,col=c(2,4,3),bty="n");dev.off();'
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#     1     634    2364    3044    3017   61248 
# MaxDP = 2572+(1.5(2816-284))=6370
# MinDP = 2572-(1.5(2816-284))=0 (negative changed to 0)

Rscript -e 'a=read.table("tmp.dp.mtDNA");b=summary(a$V1);print(b);'\
'pdf("tmp.dp.mtDNA.histogram.pdf",width=10,height=5);'\
'hist(as.numeric(a$V1),breaks=1000,main="mtDNA",xlab="DP");'\
'abline(v=c(2*b[4],b[4]+c(1,-1)*1.5*(b[5]-b[2])),col=c(2,4,3));'\
'legend("topright",c("2*mean","mean+1.5*interquartilerange","mean-1.5*interquartilerange"),lwd=1,col=c(2,4,3),bty="n");dev.off();'
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   2545   53884   54852   52880   55475   61186 
# MaxDP = 52880+(1.5(55475-53884))=55267
# MinDP = 52880-(1.5(55475-53884))=50494
