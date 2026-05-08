#!/bin/bash

# The previous line is mandatory

#SBATCH --job-name=picard_bai    # Name of your job
#SBATCH --cpus-per-task=1       # Number of cores to reserve
#SBATCH --mem-per-cpu=1G        # Amount of RAM/core to reserve
#SBATCH --time=00:10:00         # Maximum allocated time
#SBATCH --qos=30min             # Selected queue to allocate your job
#SBATCH --output=picard_bai.o%j.%a    # Path and name to the file for the STDOUT
#SBATCH --error=picard_bai.e%j.%a     # Path and name to the file for the STDERR
#SBATCH --array=1-113            # Adjust this based on the number of BAM files

# Picard was installed in the cutadapt conda environment
conda activate cutadapt

cd ./4_MD/

# Directory paths
input_dir="./4_MD/"
output_dir="./4_MD/"

# Get the sample name 
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sample_names.txt)
echo "Current sample name: $sample"

#  Construct the input BAM file path
input_bam="${input_dir}${sample}.bam"

# Construct the output BAM file path
output_bam="${output_dir}${sample}.bam"

# If the input BAM file does not exist, exit
if [ ! -f "$input_bam" ]; then
    echo "BAM file missing for sample $sample"
    exit 1
fi

# Run picard BuildBamIndex
export _JAVA_OPTIONS="-Xmx1g"

picard BuildBamIndex \
    USE_JDK_DEFLATER=true \
    USE_JDK_INFLATER=true \
    INPUT="$input_bam"

conda deactivate
