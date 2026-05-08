#!/bin/bash

# The previous line is mandatory

#SBATCH --job-name=picard_md    # Name of your job
#SBATCH --cpus-per-task=1       # Number of cores to reserve
#SBATCH --mem-per-cpu=1G        # Amount of RAM/core to reserve
#SBATCH --time=02:00:00         # Maximum allocated time
#SBATCH --qos=6hours            # Selected queue to allocate your job
#SBATCH --output=picard_md.o%j.%a    # Path and name to the file for the STDOUT
#SBATCH --error=picard_md.e%j.%a     # Path and name to the file for the STDERR
#SBATCH --array=55 #wc -l sample_names.txt

# Picard was installed in cutadapt conda environment
conda activate cutadapt

cd ./4_MD

# Directory paths
bam_dir="./3_RG"
output_dir="./4_MD"
tmp_dir="$TMPDIR"

# Get the sample name 
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sample_names.txt)
echo $sample 

# Find BAM files for the current sample
input_files=$(find "$bam_dir" -regextype posix-extended -regex ".*/.*_${sample}(-NEW)?_.*\.bam")
# If no input files are found, exit
if [ -z "$input_files" ]; then
    echo "BAM files missing for sample $sample"
    exit 1
fi

# Generate INPUT parameters for Picard
input_params=""
for bam_file in $input_files; do
    input_params+=" INPUT=$bam_file"
done

# Define output files
output_bam="${output_dir}/${sample}.bam"
metrics_file="${output_dir}/${sample}.txt"

# Skip if output BAM already exists
if [ -f "$output_bam" ]; then
    echo "Output BAM already exists for $sample, skipping..."
    exit 0
fi
# Construct the Picard command
picard_cmd="picard MarkDuplicates $input_params \
    OUTPUT=\"$output_bam\" \
    METRICS_FILE=\"$metrics_file\" \
    TAGGING_POLICY=All \
    REMOVE_DUPLICATES=true \
    MAX_RECORDS_IN_RAM=300000 \
    USE_JDK_DEFLATER=true \
    USE_JDK_INFLATER=true \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
    TMP_DIR=\"$tmp_dir\""

# Print the Picard command for checking
echo "Constructed Picard command:"
echo "$picard_cmd"

# A minority samples require larger memory to process
# Using JAVA_OPTIONS to allocate more memory for only these samples
#export _JAVA_OPTIONS="-Xmx35g" 

eval "$picard_cmd"

conda deactivate
