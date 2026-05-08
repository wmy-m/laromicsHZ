#!/bin/bash

# The previous line is mandatory

#SBATCH --job-name=bwa    # Name of your job
#SBATCH --cpus-per-task=20    # Number of cores to reserve
#SBATCH --mem-per-cpu=5G     # Amount of RAM/core to reserve
#SBATCH --time=02:00:00      # Maximum allocated time
#SBATCH --qos=6hours         # Selected queue to allocate your job
#SBATCH --output=bwa2_array.o%j.%a   # Path and name to the file for the STDOUT
#SBATCH --error=bwa2_array.e%j.%a    # Path and name to the file for the STDERR
#SBATCH --array=1-375 

# BWA-mem2 is installed in the cutadapt conda environment
conda activate cutadapt

cd ./2_bwa
# Read the sample name corresponding to the current array task ID 2150%50
sample_name=$(sed -n "$((SLURM_ARRAY_TASK_ID))p" sample_names.txt)

echo "Current sample name: $sample_name"

# Directory containing input fastq files
input_dir="./1_cutadapt"

# Directory for output BAM files
output_dir="./2_bwa"

# Find the corresponding R1 and R2 files
R1_file="${input_dir}/${sample_name}_R1.fastq.gz"
R2_file="${input_dir}/${sample_name}_R2.fastq.gz"

echo "Processing sample: $sample_name"
echo "R1 file: $R1_file"
echo "R2 file: $R2_file"

# Check if both R1 and R2 files exist
if [ -f "$R1_file" ] && [ -f "$R2_file" ]; then

    # Run bwa-mem2 and samtools commands
    bwa-mem2 mem -t 20 -M ./GCF_964199755.1_bLarMic1.1_genomic_chr.fna \
    "$R1_file" "$R2_file" | \
    samtools view -bhq 20 -@ 20 - | \
    samtools sort -m 5G -@ 20 -o ./2_bwa/${sample_name}.bam
    
else
    echo "R1 or R2 file is missing for sample $sample_name :("
fi

conda deactivate
