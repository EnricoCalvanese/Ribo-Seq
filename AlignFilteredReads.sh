#!/bin/bash
#SBATCH --job-name=bowtie2_mapping
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=bowtie2_mapping_%A_%a.out
#SBATCH --error=bowtie2_mapping_%A_%a.err
#SBATCH --array=0-5

# Set up input and output directories
REFERENCE_DIR="/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference"
INPUT_DIR="/global/scratch/users/enricocalvane/riboseq/Xu2017/rRNA_filtered"
OUTPUT_DIR="/global/scratch/users/enricocalvane/riboseq/Xu2017/bowtie2_aligned"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Array of input files
INPUT_FILES=(SRR4203374_non_rRNA_tRNA.fastq.gz SRR4203375_non_rRNA_tRNA.fastq.gz \
             SRR4203376_non_rRNA_tRNA.fastq.gz SRR4203377_non_rRNA_tRNA.fastq.gz \
             SRR4203378_non_rRNA_tRNA.fastq.gz SRR4203379_non_rRNA_tRNA.fastq.gz \
             SRR4203380_non_rRNA_tRNA.fastq.gz SRR4203381_non_rRNA_tRNA.fastq.gz)

# Select current file based on array task ID
INPUT_FILE=${INPUT_FILES[${SLURM_ARRAY_TASK_ID}]}
SAMPLE_NAME=$(basename ${INPUT_FILE} _non_rRNA_tRNA.fastq.gz)

# Load Bowtie2 module (adjust based on your system's module system)
module load bowtie2
module load samtools

# Build Bowtie2 index (if not already done)
if [ ! -f "${REFERENCE_DIR}/TAIR10_index.1.bt2" ]; then
    bowtie2-build ${REFERENCE_DIR}/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa ${REFERENCE_DIR}/TAIR10_index
fi

# Perform Bowtie2 alignment
bowtie2 -x ${REFERENCE_DIR}/TAIR10_index \
        -U ${INPUT_DIR}/${INPUT_FILE} \
        -S ${OUTPUT_DIR}/${SAMPLE_NAME}_aligned.sam \
        --very-sensitive \
        -N 1 \
        -L 20 \
        --np 0 \
        --no-unal \
        -p ${SLURM_CPUS_PER_TASK}

# Convert SAM to BAM and sort
samtools view -bS ${OUTPUT_DIR}/${SAMPLE_NAME}_aligned.sam | \
samtools sort -o ${OUTPUT_DIR}/${SAMPLE_NAME}_sorted.bam

# Index the BAM file
samtools index ${OUTPUT_DIR}/${SAMPLE_NAME}_sorted.bam

# Optional: Remove intermediate SAM file to save space
rm ${OUTPUT_DIR}/${SAMPLE_NAME}_aligned.sam

# Print completion message
echo "Alignment completed for ${SAMPLE_NAME}"
