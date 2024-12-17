#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=02:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=ribo_trim
#SBATCH --output=trim_%j.out
#SBATCH --error=trim_%j.err

# Set up directories
INPUT_DIR="/global/scratch/users/enricocalvane/riboseq/Xu2017"
OUTPUT_DIR="/global/scratch/users/enricocalvane/riboseq/Xu2017/trimmed"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Process each fastq file
for file in ${INPUT_DIR}/*.fastq; do
    # Extract base filename
    base=$(basename ${file} .fastq)
    
    # Run Trim Galore
    trim_galore \
        --fastqc \
        --length 25 \
        --trim-n \
        --output_dir ${OUTPUT_DIR} \
        ${file}
done

# Generate MultiQC report
multiqc ${OUTPUT_DIR} -o ${OUTPUT_DIR}/multiqc_report

echo "Trimming and quality control completed."
