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
**Clinal analyses**  
Phenotype-geographic clines were fitted using `HZAR` for shortlisted traits. The `HZAR` pipeline was turnt into a R function to run selected traits as an array in slurm. The phenotype dataset was subsetted to only containing inidivduals from the genomic clinal dataset.
```
sbatch HZAR_pheno.slurm
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
Population genomic structure was explored using PCA with `PLINK`, `Admixture` and `TriangulaR`. All scripts in this section are filed under __Popgen__ folder. Only the autosomes from the hybrid zone dataset (1) were used for these analyses.

`TriangulaR` was used to identify ancestry-informative markers (AFD=0.5 to 1) and using those markers to calculate hybrid indices, interclass heterozygosity and building triangle plots. A population map was created with all the samples being tested, with the allopatric populations set to the parental populations.  

```
sbatch 01_TriangulaR.R
```

`PLINK` was used to run PCA and prepare the input file for `Admixture`. The dataset is first splited into the two hybrid zones and before using `PLINK`. 
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
`SMC++` was used to estimate effective population size of the four hybrid zone species. We further estimated split times for each hybrid zone. This analysis was done using the Demographic dataset (2), with 10 individuals per species. The scripts are filed under __Demo__ folder.

We first readied the dataset by filtering for biallelic SNPs, max. 10% missing genotypes and autosomes only. We created a missingness mask to omit regions where more than 10% of individuals have missing data. 

```
# Filter VCF file for SNPs: biallelic SNPs, max. 10% missing genotypes, autosomes only
sbatch 01_smcppfilt.slurm

# Create missingness mask 
# b) Mask for covered sites: 90% of individuals covered, minus mappability mask
# (1b) Generate BED file for each individual with at least 5x coverage
sbatch 02_smcppbamdepth.slurm

# Overlap bed files to get a mask where <90% of individuals are covered, combine with missingness mask
sbatch 03_smcppbedoverlap.slurm

# Recode vcf file into vcf-iid
# This will cause the chromosome names to go a bit haywire
# Chr name in BED file needs to renamed

plink --vcf Demo.filt.vcf.gz --keep-allele-order --allow-extra-chr \
    --chr-set 95 no-xy no-mt --recode vcf-iid --out Demo.filt.recoded
    bgzip -c Demo.filt.recoded.vcf > Demo.filt.recoded.vcf.gz
    tabix Demo.filt.recoded.vcf.gz

zcat Demo.miss.bed.gz | sed 's/^chr0*//' > Demo.miss.fixed.bed
bgzip Demo.miss.fixed.bed
tabix Demo.miss.fixed.bed.gz

# Convert VCF files to smcpp input files and run smcpp
sbatch 04_smcppvcftoinfile.slurm
sbatch 05_smcppvcftoindist.slurm #This is ran four times for each species
sbatch 06_smcppvcftoinfilepairs.slurm #for 2 species 

```

We then ran smc++ estimation on 10 distinugished individiual
```
sbatch 07_smcppestd10.slurm

# Plotting the results
for i in Arg Cac Gla Occ; do
    smc++ plot -g 10 -x 1e4 1e7 runLar${i}d10/runLar${i}d10.plot.pdf runLar${i}d10/model.final.json
    done
```
We followed up by running smc++ split time estimation using the paired files. 
```
sbatch 08_smcppsplit.slurm

# plot split 
    for i in LarArg.LarCac LarGla.LarOcc; do
    smc++ plot -g 10 -x 1e4 1e7 run${i}002/run${i}002.plot.pdf run${i}002/model.final.json
    done
```
## Genomic scans 
To investigate genomic architecture, genomic scans were carried out using FST and genomic cline parameters. 
FST was calculated separately for each chromosome class using a custom script for Weir and Cockerham’s FST. The Allopatric dataset was used for the FST analyses. The scripts for FST are filed under __Popgen__ folder. 

First, the ploidy for chrW and chrZ was forced to diploid notation using `BCFtools`. 
```
ls *ChrW* *ChrZ* > fixploidy_list.txt
sbatch 04_fixploidy.slurm
```
The custom script calculates locus-based FST with the numerator and denominator values separately. The locus-based FST was then binned into different window sizes for better visualisation. 
```
sbatch 05_wcfst.slurm

#Calculating the average FST values from the locus-based FST results
##global FST
cat EU.auto.wcfst.txt EU.ChrW.forced.wcfst.txt EU.ChrZ.forced.wcfst.txt| \
awk '{if($3!="-nan"){a+=$3;b+=$4}}END{print a/b}'
cat NA.auto.wcfst.txt NA.ChrW.forced.wcfst.txt NA.ChrZ.forced.wcfst.txt| \
awk '{if($3!="-nan"){a+=$3;b+=$4}}END{print a/b}'

##autosome FST
awk '{if($1 != "chrZ" && $1 !~ /^chrW/ && $3 != "-nan") {a+=$3; b+=$4}} END {print a/b}' < NA.auto.wcfst.txt
awk '{if($1 != "chrZ" && $1 !~ /^chrW/ && $3 != "-nan") {a+=$3; b+=$4}} END {print a/b}' < EU.auto.wcfst.txt

##Z FST
awk '{if($1 == "chrZ" && $3 != "-nan") {a+=$3; b+=$4}} END {print a/b}' < NA.ChrZ.forced.wcfst.txt
awk '{if($1 == "chrZ" && $3 != "-nan") {a+=$3; b+=$4}} END {print a/b}' < EU.ChrZ.forced.wcfst.txt

## Chr W
awk '{if($1 == "chrW" && $3 != "-nan") {a+=$3; b+=$4}} END {print a/b}' < NA.ChrW.forced.wcfst.txt
awk '{if($1 == "chrW" && $3 != "-nan") {a+=$3; b+=$4}} END {print a/b}' < EU.ChrW.forced.wcfst.txt

# Making windowed FST outputs
sbatch 06_winfst.slurm
```

Genomic clines were fitted across the genome using `bgchm` to observe changes in the cline slope and center. We used the clinal dataset (4) for this analysis. The dataset was first filtered to only retain ancestry informative SNPs (AFD threshold = 0.5) before running `bgchm`. The script for this analysis was filed under __bgchm__. 

VCFs that were filtered to AFD threshold = 0.5 was exported from TriangulaR script (vcf.diff), and their positions were retrieved. The clinal datasets were subsetted to these postions.
```
# extracting SNPs with AFD > 0.5
sbatch 01_extractPos.slurm

# subset clinal datasets to shortlisted SNPs
sbatch 02_extractSNPs.slurm
```
Input files for `bgchm` were generated from the vcf files. The vcf was first converted to 012 format using `VCFtools` and then splited in P0, P1 and hybrids population. 

```
# Making 012 input
sbatch 03_vcf012.slurm

# creating P0 and P1 population file from triangulaR population map
# these two files are under the popgen folder
awk '$2 == "ARGE" {print $1}' ./EUpopmap.txt > EU_P0.txt
awk '$2 == "CACH" {print $1}' ./EUpopmap.txt > EU_P1.txt

awk '$2 == "GLAUC" {print $1}' ./NApopmap.txt > NA_P0.txt
awk '$2 == "OCCI" {print $1}' ./NApopmap.txt > NA_P1.txt

# creating input files for EU
# this is the same for NA, just replace the file_prefizes, P0 and P1 files
file_prefixes=("EU.ChrZ.Cline" "EU.auto.Cline" "EU.ChrW.Cline")
    for prefix in "${file_prefixes[@]}"; do
        # Create dynamic header from .012.pos file and store it in a variable
        header=$(awk 'BEGIN {printf "Sample "} {printf "%s_%s ", $1, $2} END {print ""}' ${prefix}.012.pos)

        # Create the main .txt file with the header
        echo "$header" > ${prefix}.txt
        awk 'NR==FNR {names[NR]=$1; next} { $1 = names[FNR]; print }' ${prefix}.012.indv ${prefix}.012 >> ${prefix}.txt

        # Add header and filter out entries using P1.txt
        echo "$header" > ${prefix}_P1.txt
        grep -f EU_P1.txt ${prefix}.txt >> ${prefix}_P1.txt

        # Add header and filter out entries using P0.txt
        echo "$header" > ${prefix}_P0.txt
        grep -f EU_P0.txt ${prefix}.txt >> ${prefix}_P0.txt

        # Add header and filter out hybrids
        grep -vFf EU_P1.txt -vFf EU_P0.txt ${prefix}.txt >> ${prefix}_hybrids.txt

        # Print row and column counts for the main file
        awk 'END {print "Rows: " NR; print "Columns: " NF}' ${prefix}.txt

        # Print row and column counts for the P0 file
        awk 'END {print "Rows: " NR; print "Columns: " NF}' ${prefix}_P0.txt

        # Print row and column counts for the P1 file
        awk 'END {print "Rows: " NR; print "Columns: " NF}' ${prefix}_P1.txt

        # Print row and column counts for the hybrids file
        awk 'END {print "Rows: " NR; print "Columns: " NF}' ${prefix}_hybrids.txt
    done
```
All loci was fitted using `bgchm` in a scalable, parallelizable manner as [suggested](https://github.com/zgompert/bgc-hm). This was ran separately for the chromosome classes.

```
# Creating a dataframe by merging all 3 files (P0, P1, hybrids) into one matrix
sbatch 00_input.slurm
# Estimation of cline SDs and hybrid indices
sbatch 01_manyloci.slurm
# Estimate clines for all of the loci in parallel
sbatch 02_fitSnps.slurm
# Combine estimates from each batch (for autosomes and chrZ)
# the outputs were first moved to a new folder for each hybrid zone
# chrW can be run without parallelisation for EU 
sbatch 03_merge.slurm
```

## Genome-wide association studies
GWAS was carried out using `GEMMA` to locate regions of interest that could be involved in phenotype. We used the standardised trait values from the phenotypic analyses. See ./Pheno/PhenoInput.R for input file generation. The Trait mapping dataset (5) was used for this section. This is all individuals except for the ones from the allopatric population. The scripts are filed under __popgen__ folder

We first made the input using `PLINK` before running the lmm model in `GEMMA` with correction for relatedness. 

```
# Make input file using PLINK
sbatch 07_GEMMAInput_EU.slurm
sbatch 07_GEMMAInput_NA.slurm

# Run GEMMA for lmm model
sbatch 08_GEMMA_EU.slurm
sbatch 08_GEMMA_NA.slurm

```
The results were plotted in R to quickly locate traits with regions of interest

```
# Plotting loop
files <- list.files(path = "./", pattern = "*.gemma.assoc.txt.reduced.gz", full.names = TRUE)
# Loop through each file
for (file in files) {
  # Read in the data
  b <- read.table(file, h = TRUE)
  
  # Calculate Bonferroni and FDR adjusted p-values
  b <- cbind(b,
             p_wald_bonferroni = p.adjust(b$p_wald, "bonferroni"),
             p_wald_fdr = p.adjust(b$p_wald, "fdr"))
  # Remove non-finite or zero p-values
  valid <- is.finite(b$p_wald) & b$p_wald > 0
  b <- b[valid, ]
  
  # Skip file if no valid p-values
  if (nrow(b) == 0) next
  
  significant_threshold <- -log10(0.05 / nrow(b))
  plot_filename <- paste0("plot_", gsub(".*/|\\.gemma\\.assoc\\.txt\\.reduced\\.gz", "", file), ".png")
    png(plot_filename, width = 10, height = 5, res = 300, units = "in")
  # Create the basic plot
  plot(1:nrow(b), -log10(b$p_wald),
       type = "n", xaxt = "n", xlab = "Chromosomes", ylab = expression(log[10](P)),
       ylim = c(0, max(-log10(b$p_wald), na.rm = TRUE)), main = paste(basename(file)))
  # Add points with color based on significance
  points(1:nrow(b), -log10(b$p_wald),
         col = c("#00000044", "#FF880044")[1 + 1 * (-log10(b$p_wald) > significant_threshold)],
         pch = 16)
    axis(1, labels = FALSE, tck = 0)
  # Close the plot device
  dev.off()
    cat("Completed plot for", file, "\n")
}
```

## Clinal analyses for regions of interest
Geographic clines based on hybrid indices derived from trait loci versus the genome background were used to investigate regions that gave a peak in GWAS analyses. These regions were identified from the GWAS analyses and listed in `TraitRegions.txt`. All scripts from this section are filed under __TraitRegions__ folder.

These regions were extracted from the clinal datasets using `BCFtools` and their hybrid indices was calculated using `TriangulaR`. 

```
# Extract SNPs from region of interest
# this was ran for separately for each hybrid zone
sbatch 00_ExtractRegions.slurm

# Calculating hybrid indices
01_HI.R

# Running HZAR for 6 models
# 3 chain, 3 runs, chain length = 100000 
./HZAR/02_runHZAR_EU.R
./HZAR/02_runHZAR_NA.R
```
Clines from the trait regions will be compared to 1000 randomly selected background loci following [Schield et al. 2024](https://github.com/drewschield/hirundo_speciation_genomics/) pipeline. 

```
ml load BEDTools/2.31.0-GCC-12.3.0
ml load BCFtools/1.18-GCC-12.3.0 #a different version

bcftools query -f '%CHROM\t%POS\n' ./5_bcftools/datasets/Chapter1.NA_Cline.auto.new.bSNPs.2.vcf.gz | \
grep -v '^NW_' |awk '{OFS="\t"}{print $1,$2-1,$2}' | bedtools intersect -v -wb -a - -b candidate.bed | \
shuf -n 1000 | awk '{OFS="\t"}{print $1,$3}' > background.NA.snps.txt

bcftools query -f '%CHROM\t%POS\n' ./5_bcftools/datasets/Chapter1.EU_Cline.auto.new.bSNPs.2.vcf.gz | \
grep -v '^NW_' | awk '{OFS="\t"}{print $1,$2-1,$2}' | bedtools intersect -v -wb -a - -b candidate.bed | \
shuf -n 1000 | awk '{OFS="\t"}{print $1,$3}' > background.EU.snps.txt

bcftools query -f '%CHROM\t%POS\n' ./5_bcftools/datasets/Chapter1.NA_Cline.auto.new.bSNPs.2.vcf.gz  > all.NA.snps.txt
bcftools query -f '%CHROM\t%POS\n' ./5_bcftools/datasets/Chapter1.EU_Cline.auto.new.bSNPs.2.vcf.gz  > all.EU.snps.txt

pop="NA"
outfile=background.$pop.snps.region.txt
rm -f $outfile
while read region; do
	chrom=`echo "$region" | cut -f 1`;
	snp=`echo "$region" | cut -f 2`;
	tmp=`grep $chrom all.$pop.snps.txt | grep -w -B49 -A50 $snp`
	start=`echo $tmp | cut -d' ' -f2`
	end=`echo $tmp | rev | cut -d' ' -f 1 | rev`
	echo -e "$chrom\t$start\t$end" >> $outfile
done < background.$pop.snps.txt

pop="EU"
outfile=background.$pop.snps.region.txt
rm -f $outfile
while read region; do
	chrom=`echo "$region" | cut -f 1`;
	snp=`echo "$region" | cut -f 2`;
	tmp=`grep $chrom all.$pop.snps.txt | grep -w -B49 -A50 $snp`
	start=`echo $tmp | cut -d' ' -f2`
	end=`echo $tmp | rev | cut -d' ' -f 1 | rev`
	echo -e "$chrom\t$start\t$end" >> $outfile
done < background.$pop.snps.txt

# Make vcf files for each background region
sbatch 03_ExtractRegions_bg.slurm

# Get HI for each background region
04_HI.R

# Run HZAR
sbatch ./HZAR/01_runHZAR.bg.EU.slurm
sbatch ./HZAR/01_runHZAR.bg.NA.slurm
```





