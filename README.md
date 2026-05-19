# laromicsHZ

This repository contains details on the data processing and analysis steps used to study the <i> Larus argentatus </i> / <i> Larus cachinnans </i> and 
<i> Larus glaucscens </i> / <i> Larus occidentalis </i> hybrid zone through looking at the genomic divergence of involved species and the genetic architecture of
phenotypic species differences. Calculations and bioinformatic analyses were performed at sciCORE (http://scicore.unibas.ch/) scientific computing center at University of Basel. 

## Phenotypic data analyses
Phenotypic data was retrieved from past studies (NA: Bell 1996 The Condor; EU: Gay et al. 2007 Molecular Ecology, Gat et al. 2009 Heredity, Neubauer et al. 2009 The Auk) and merged with our recent fieldwork data. The R scripts are filed under __Pheno__ folder.

Size traits were first standarisded and PCA was conducted on the phenotype dataset using `prcomp` in R. Linear discriminant analysis was carried out, with allopatric population used as test population. This was ran on a local computer. 

```
Pheno.NA.R
Pheno.EU.R
```

## Initial genomic data processing
All scripts for this section are filed under __Bioinfo__ folder
### Adaptor trimming and filtering 
`Cutadapt` was used to trim adaptors and filter read quality

**Run Cutadapt**  
This was ran as an array in slurm. All outputs are stored in ./1_cutadapt.
```
# Make a list of R1 fastq files
find . -name "*_R1_001.fastq.gz" > R1_files_list.txt

# Run Cutadapt
sbatch 01_cutadapt_array.slurm
```
### Read alignment using BWA and SAMtools
`BWA-mem2` was used to aligned the reads to the _Larus michahellis_ reference genome [GCF_964199755.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_964199755.1/). `SAMtools` was used for filtering and sorting the bam files. 

**Index reference genome**  
The refseq assembly does not include the mitochondrion genome. The reads are also aligned to the mitochondrion scaffold in the GenBank assembly. 
```
bwa-mem2 index GCF_964199755.1_bLarMic1.1_genomic_chr.fna
bwa-mem2 index GCA_964199755.1_bLarMic1.1_genomic_chr_MT.fa
```
**Run BWA-mem2 and SAMtools** 
All outputs are stored in ./2_bwa. 
```
# Make a list of R1 file names
# These file names do not include R1 or R2 information.
find ./1_cutadapt -name "*_R1.fastq.gz" | sed 's/_R1.fastq.gz//' \
| xargs -n 1 basename > sample_names.txt

# Run BWA-mem2 and Samtools
sbatch 02_bwa_array.slurm
```
### Merging samples sequenced on multiple lanes
Each sample was sequenced on multiple lanes. `Picard` was used to add read groups and then merge the bam files of the same sample together. 

**Add read group with Picard AddOrRemoveReadGroups**  
For RGPU, I have used flowcell.lane instead of instrument.flowcell for the ease of coding. The logic remains the same. All outputs are stored in ./3_RG.
```
# Make a list of bam files
find ./2_bwa -name "*.bam" > bam_files_list.txt

# Run Picard AddOrRemoveReadGroups
sbatch 03_RG_array.slurm
```
**Merge samples with Picard MarkDuplicates**  
Samples were merged and duplicates were also marked. The script was ran in ./4_MD with the finalised bam and bai files moved to ./01_bam_refMic
```
# Extract unique sample names from BAM file names and save to a text file
find /scicore/home/marque0000/GROUP/wu0006/3_RG -name "*.bam" | awk -F'_' '{print $5}' | sort | uniq > sample_names.txt

# Run Picard MarkDuplicates
sbatch 04_MD_array.slurm

# Index new bam files
sbatch 05_picard_bai_array.slurm
```
### Quality check on bam files using Qualimap and SAMtools
`SAMtools quickcheck` was used to check validity of the files. `Qualimap` was used to check overall bam file quality, with the coverage and mapping quality being extracted. All outputs are stored in ./4_MD.
for a closer inspection. 

1. SAMtools quickcheck 
```
# Checking validity of bam files
cd 01_bam_refMic
samtools quickcheck -v ./*.bam
```
2. Qualimap  
If the bam files passed `quickcheck`, it will be processed by `Qualimap`.
```
# Run Qualimap
sbatch 06_qualimap_array.slurm

# check if all runs finished properly
grep -r "Finished" qualimap_bamqc.o* | wc -l
```
Mapping percentage, mean mapping quality and mean coverage was extracted from the `Qualimap` report of each samples and summarised. 
```
# Make summary file
output_file="qualimap20260421_summary_final.txt"

# Loop through each sample directory
for sample_dir in */; do
    # Extract the sample name 
    sample_name=$(basename "$sample_dir" | sed 's/_results$//')
    
    # Define the path to the file inside the sample directory
    text_file="${sample_dir}/genome_results.txt"  

    # Extract the required lines from the text file
    num_mapped_reads=$(grep -E '^[[:space:]]*number of mapped reads =' "$text_file" | awk -F'= ' '{print $2}' | awk '{print $1}' | tr -d ',')
    percentage_mapped=$(grep -E '^[[:space:]]*number of mapped reads =' "$text_file" | awk -F'[()]' '{print $2}') # Extract the percentage
    mean_mapping_quality=$(grep -E '^[[:space:]]*mean mapping quality =' "$text_file" | awk -F'= ' '{print $2}')
    mean_coverage_data=$(grep -E '^[[:space:]]*mean coverageData =' "$text_file" | awk -F'= ' '{print $2}')

    # Check if any of the extracted values are missing
    if [[ -z "$num_mapped_reads" || -z "$percentage_mapped" || -z "$mean_mapping_quality" || -z "$mean_coverage_data" ]]; then
        echo "Missing values for sample $sample_name" >&2
    fi

    # Append the data to the output file
    echo -e "${sample_name}\t${num_mapped_reads}\t${percentage_mapped}\t${mean_mapping_quality}\t${mean_coverage_data}" >> $output_file
done

# Calculate the average coverage
average=$(awk -F'\t' 'NR>1 {sum+=$5; count++} END {if (count > 0) print sum/count}' $output_file)

# Append the average to the output file
echo -e "\nAverage Coverage Data:\t$average" >> $output_file
```

## Variant calling
All scripts for this section are filed under __VariantCall__ folder
Variants are called using `bcftools mpileup`. There are five datasets for different analyses: 1 – Hybrid zone;  2. – Demography; 3 – Allopatric; 4 – Cline; 5 – Trait mapping. SNPs were first called for all samples in this study and filtered. The five datasets were subset from this original SNP call. Additional filters were applied after subsetting. All outputs are stored in ./5_bcftools and subdirectory ./5_bcftools/datasets for the subsets.

**Determine genetic sex**  
Coverage depth of the first 15 autosomes were compared with Chr Z to determine genetic sex. Females should be halved of the autosome depth because of the ZW system. 

```
for i in ./qualimap; do
for j in $(ls $i | grep -v summary); do
echo ${j%_results} $(grep "mean coverageData = " $i"/"$j"/genome_results.txt" | sed 's/.*\= \(.*\)X/\1/g') $(grep "chrZ" $i"/"$j"/genome_results.txt" | cut -f 5) | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$3/$2}'
done; done > sexing_ZtomeanDPratio.txt

or i in ./qualimap; do 
  for j in $(ls $i | grep -v summary); do
    echo ${j%_results} \
      $(awk '$1=="chr01"{print $4}' $i/$j/genome_results.txt) \ 
      $(awk '$1=="chr02"{print $4}' $i/$j/genome_results.txt) \
      $(awk '$1=="chr03"{print $4}' $i/$j/genome_results.txt) \
      $(awk '$1=="chr04"{print $4}' $i/$j/genome_results.txt) \
      $(awk '$1=="chr05"{print $4}' $i/$j/genome_results.txt) \
      $(awk '$1=="chr06"{print $4}' $i/$j/genome_results.txt) \
      $(awk '$1=="chr07"{print $4}' $i/$j/genome_results.txt) \
      $(awk '$1=="chr08"{print $4}' $i/$j/genome_results.txt) \
      $(awk '$1=="chr09"{print $4}' $i/$j/genome_results.txt) \
      $(awk '$1=="chr10"{print $4}' $i/$j/genome_results.txt) \
      $(awk '$1=="chr11"{print $4}' $i/$j/genome_results.txt) \
      $(awk '$1=="chr12"{print $4}' $i/$j/genome_results.txt) \
      $(awk '$1=="chr13"{print $4}' $i/$j/genome_results.txt) \
      $(awk '$1=="chr14"{print $4}' $i/$j/genome_results.txt) \
      $(awk '$1=="chr15"{print $4}' $i/$j/genome_results.txt) \
      $(awk '$1=="chrZ"{print $4}' $i/$j/genome_results.txt) | \
awk 'BEGIN{OFS="\t"}{for(i=2;i<=16;i++){a+=$i};a=a/15;print $1,a,$17,$17/a}'
done; done > sexing_Ztomean15chrratio.txt
```
**Creating paritions to parallelise variant calling**  
Create a bed file for windows of 2Mbp using the index. This is to speed up variant calling by running multiple windows simultaneously. 
```
module load BEDTools/2.31.0-GCC-12.3.0
bedtools makewindows -g GCF_964199755.1_bLarMic1.1_genomic_chr.fna.fai -w 2000000 > bLarMic1.2Mbpwin.bed
```

**Running bcftools mpileup**  
1. Call SNPs for each of the paritions
```
sbatch 01_mpileup.slurm
```
2. Concatenate the paritions to chromosome levels
```
sbatch 02_concatChr.slurm
```

**Variant filtering**
`BCFtools` were used for filtering. There are two rounds of filtering. The first filtering contain quality filters, the second filtering was tailored to the anaylses. 

1. Determine depth filters. Information from '%INFO/DP' of the vcf files were extarcted from each chromosome class and depth limits were calculated in R. Minimium and maximum depth were defined as mean ± (1.5(interquartile range).
```
sbatch 03_bcfDepth.slurm
```
2. Initial filtering for quality
Indels, SNPs closed to indels, SNPs where QUAL <20, multiallelic sites, invariant sites and sites with >50% missing data were removed. Genotypes with <4 reads or GQ <20 were set to missing. Maximum and minimum depth filters were also applied.
```
sbatch 04_filtSNP.slurm
```
4. Concatenate all autosomes to make downstream analyses easier
```
sbatch 05_concatAuto.slurm
```
5. Removing the male individuals from ChrW so the males genotype will not be treated as missing data in ChrW.
```
sbatch 06_RemoveM.slurm
```
6. Subsetting to specific dataset
The samples in each dataset are given in supplementary table S4. The subsetting for Dataset (2) demography was done separately later (see below). 
```
sbatch 07_Subset.slurm
```
8. Second filter specific to the analyses
Multiallelic sites, non-polymorphic sites, sites with >50% missing data and MAF < 0.05 were removed for each dataset. For ChrW, the missing data filter was relaxed to >75%.  
```
sbatch 08_filtSNP_R2.slurm
```

## Basic population genomics analyses
Population genomic structure was explored using PCA with `PLINK`, `Admixture` and `TriangulaR`. All scripts in this section are filed under __Popgen__ folder. Only the autosomes were used for these analyses.

`TriangulaR` was used to identify ancestry-informative markers (AFD=0.5 to 1) and using those markers to calculate hybrid indices, interclass heterozygosity and building triangle plots. The allopatric populations are used as the parental populations.  

```
sbatch 01_TriangulaR.R
```

`PLINK` was used to run PCA and prepare the input file for `Admixture`. The dataset is first splited into the two hybrid zones and before `PLINK`. 
```
sbatch 02_Subset.slurm
```
In `PLINK`, the missingness was calculated and pca was ran using `--missing` and `--pca` respectively. The vcf was further filtered by removing missing genotypes, a minor allele count is observed in at least 2 individuals and keeping only one SNP every 10,000 bp for `ADMIXTURE`.
```
sbatch 03_PLINK.slurm
```
`ADMIXTURE` was ran using the `PLINK` output. The chromosome names were first removed to keep them all as intergers so `ADMIXTURE` can read it. The best _K_ was then identified.
```
awk '{$1="0";print $0}' Chapter1.NA_Full.auto.new.bSNPs.2_adm.bim > Chapter1.NA_Full.auto.new.bSNPs.2_adm.bim.tmp
mv Chapter1.NA_Full.auto.new.bSNPs.2_adm.tmp Chapter1.NA_Full.auto.new.bSNPs.2_adm.bim

awk '{$1="0";print $0}' Chapter1.EU_Full.auto.new.bSNPs.2_adm.bim > Chapter1.EU_Full.auto.new.bSNPs.2_adm.bim.tmp
mv Chapter1.EU_Full.auto.new.bSNPs.2_adm.tmp Chapter1.EU_Full.auto.new.bSNPs.2_adm.bim

for K in 2 3 4 5 6 7; \
do admixture --cv Chapter1.EU_Full.auto.new.bSNPs.2_adm.bed $K | tee logEU${K}.out; done
grep -h CV logEU*.out

for K in 2 3 4 5 6 7; \
do admixture --cv Chapter1.NA_Full.auto.new.bSNPs.2_adm.bed $K | tee logNA${K}.out; done
grep -h CV logNA*.out
```

## Demographic history

## Genomic scans 
To investigate genomic architecture, genomic scans were carried out using FST and genomic cline parameters. 
FST was calculated separately for each chromosome class

