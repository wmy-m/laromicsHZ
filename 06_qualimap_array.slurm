#!/bin/bash

# The previous line is mandatory

#SBATCH --job-name=qualimap_bamqc    # Name of your job
#SBATCH --cpus-per-task=1            # Number of cores to reserve
#SBATCH --mem-per-cpu=16G            # Amount of RAM/core to reserve
#SBATCH --time=02:00:00              # Maximum allocated time
#SBATCH --qos=6hours                  # Selected queue to allocate your job
#SBATCH --output=qualimap_bamqc.o%j  # Path and name to the file for the STDOUT
#SBATCH --error=qualimap_bamqc.e%j   # Path and name to the file for the STDERR
#SBATCH --array=1-113                # Adjust this based on the number of samples

# Qualimap is installed via conda
conda activate qualimap

# Directory paths
bam_dir="./4_MD"
output_dir="./qualimap"

# Get the sample name 
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sample_names.txt)
echo "Current sample name: $sample"

# Construct the input BAM file path
input_bam="${bam_dir}/${sample}.bam"

# Check if the input BAM file exists
if [ ! -f "$input_bam" ]; then
    echo "BAM file missing for sample $sample"
    exit 1
fi

# Run qualimap bamqc
qualimap bamqc -bam "$input_bam" --collect-overlap-pairs --java-mem-size=16G \
-outdir "${output_dir}/${sample}_results" 

conda deactivate
