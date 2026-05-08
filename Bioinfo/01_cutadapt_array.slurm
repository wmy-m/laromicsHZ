#!/bin/bash
#The previous line is mandatory

#SBATCH --job-name=cutadapt    #Name of your job
#SBATCH --cpus-per-task=20    #Number of cores to reserve
#SBATCH --mem-per-cpu=1G     #Amount of RAM/core to reserve
#SBATCH --time=00:30:00      #Maximum allocated time
#SBATCH --qos=30min         #Selected queue to allocate your job
#SBATCH --output=20250903_cutadapt_array.o%j.%a   #Path and name to the file for the STDOUT
#SBATCH --error=20250903_cutadapt_array.e%j.%a    #Path and name to the file for the STDERR
#SBATCH --array=1-365   # Array range (Upper limit is wc -l R1_files_list.txt)

# cutadapt is installed in a conda environment 
conda activate cutadapt 

# Output directory
output_dir="./1_cutadapt/"
cd "$output_dir" || exit

# Get the R1 file for the current array job
R1_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ./R1_files_list.txt)

# My sample name looks like 
# GFB-17416_H5LG7DSXC_4_LarCacYP6DY_S32_L004_R1_001.fastq.gz

# Extract the sample name from R1 file
full_sample_name=$(basename "$R1_file" _R1_001.fastq.gz)

# Construct the corresponding R2 file path
R2_file="${R1_file/_R1/_R2}"

echo "Current sample name: $full_sample_name"

# Check if both R1 and R2 files exist
if [ -f "$R1_file" ] && [ -f "$R2_file" ]; then
    # Check if output files already exist
    if [ -f "${full_sample_name}_R1.fastq.gz" ] && [ -f "${full_sample_name}_R2.fastq.gz" ]; then
        echo "Output files already exist for sample $full_sample_name. Skipping..."
    else
        # Run cutadapt command for this pair of R1 and R2 files with the extracted adapter sequences
        cutadapt -q 20 -m 30 --cores=0 --poly-a \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -o "${full_sample_name}_R1.fastq.gz" \
        -p "${full_sample_name}_R2.fastq.gz" \
        "$R1_file" \
        "$R2_file"
    fi
else
    echo "R1 or R2 file is missing for sample $full_sample_name"
fi

conda deactivate
