# laromicsHZ

This repository contains details on the data processing and analysis steps used to study the <i> Larus argentatus </i> / <i> Larus cachinnans </i> and 
<i> Larus glaucscens </i> / <i> Larus occidentalis </i> hybrid zone through looking at the genomic divergence of involved species and the genetic architecture of
phenotypic species differences. Calculations and bioinformatic analyses were performed at sciCORE (http://scicore.unibas.ch/) scientific computing center at University of Basel. 

## Initial data processing
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
`BWA-mem2` was used to aligned the reads to the _Larus michahellis_ reference genome [GCF_964199755.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_964199755.1/).
`SAMtools` was used for filtering and sorting the bam files.
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
For RGPU, I have used flowcell.lane instead of instrument.flowcell for the ease of coding. The logic remains the same. 
```
# Make a list of bam files
find ./2_bwa -name "*.bam" > bam_files_list.txt

# Run Picard AddOrRemoveReadGroups
sbatch 03_RG_array.slurm
```
**Merge samples with Picard MarkDuplicates**  
Samples were merged and duplicates were also marked. These finalised bam files were stored in ./01_bam_refMic
```
# Extract unique sample names from BAM file names and save to a text file
find /scicore/home/marque0000/GROUP/wu0006/3_RG -name "*.bam" | awk -F'_' '{print $5}' | sort | uniq > sample_names.txt

# Run Picard MarkDuplicates
sbatch 04_MD_array.slurm

# Index new bam files
sbatch 05_picard_bai_array.slurm
```
### Quality check on bam files using Qualimap and SAMtools
`SAMtools quickcheck` was used to check validity of the files. `Qualimap` was used to check overall bam file quality, with the coverage and mapping quality being extracted
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





