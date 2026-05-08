#!/bin/bash

# The previous line is mandatory

#SBATCH --job-name=picard_rg    # Name of your job
#SBATCH --cpus-per-task=1       # Number of cores to reserve
#SBATCH --mem-per-cpu=1G        # Amount of RAM/core to reserve
#SBATCH --time=01:00:00         # Maximum allocated time
#SBATCH --qos=6hours             # Selected queue to allocate your job
#SBATCH --output=picard_rg.o%j.%a    # Path and name to the file for the STDOUT
#SBATCH --error=picard_rg.e%j.%a     # Path and name to the file for the STDERR
#SBATCH --array=299           # wc -l bam_files_list.txt

# Picard was installed in cutadapt conda environment
conda activate cutadapt

cd ./3_RG

# Define the output directory
output_dir="./3_RG/"

# Select bam file in the array
bam_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" bam_files_list.txt)
echo "Current sample name: $bam_file"

# Extract full sample name from the BAM file
full_sample_name=$(basename "$bam_file" .bam)

# Extract specific components from the full sample name
sample=$(echo "$full_sample_name" | cut -d '_' -f 1)
instrument=$(echo "$full_sample_name" | cut -d '_' -f 3)
flowcell=$(echo "$full_sample_name" | cut -d '_' -f 4 | sed 's/^L//')

# Construct the read group parameters
rgpu="${instrument}.${flowcell}"
rglb="Lib_${sample}"
if [[ "$full_sample_name" == *"-NEW"* ]]; then
    rglb="${rglb}_2"
fi
rgsample="$sample"

# Print rgpu, rgsample, and rglb for checking
echo "RGPU: $rgpu"
echo "RGSAMPLE: $rgsample"
echo "RGLB: $rglb"

# Construct the output BAM file path
output_bam="${output_dir}${full_sample_name}.bam"

# Run the picard AddOrReplaceReadGroups command
picard AddOrReplaceReadGroups \
    I=$bam_file \
    RGPL=illumina \
    RGPU=$rgpu \
    RGLB=$rglb \
    RGSM=$rgsample \
    RGID=$rgsample \
    O=$output_bam

conda deactivate

