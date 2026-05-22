    #!/bin/bash
    #SBATCH --job-name=bcf_subset
    #SBATCH --cpus-per-task=1
    #SBATCH --mem-per-cpu=1G
    #SBATCH --time=06:00:00
    #SBATCH --qos=6hours
    #SBATCH --output=subset.o%j.%a
    #SBATCH --error=subset.e%j.%a
    #SBATCH --array=1-2 

    # Load module for dependencies
    module load BCFtools/1.19-GCC-13.2.0
    # Define path to latest bcftools version
    bcftoolspath=./99_scripts/bcftools/

    cd ./5.2_bcftools/datasets/PCA

    # Path to subset files
    subset_files=(EU.txt NA.txt)
    subset_file=${subset_files[$SLURM_ARRAY_TASK_ID-1]}  
    subset_name=$(basename $subset_file .txt)

    #Auto
    outfile_auto=Chapter1.$subset_name.auto.new.vcf.gz
    ${bcftoolspath}bcftools view -Oz -S $subset_file ./Chapter1.autosomes.bSNPs.new.2.vcf.gz -o $outfile_auto
    tabix -p vcf $outfile_auto
