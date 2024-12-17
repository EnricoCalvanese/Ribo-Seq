#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=06:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL

# Load required modules
module load bio/bowtie2/2.5.1-gcc-11.4.0
module load bio/samtools/1.17-gcc-11.4.0

# Set directories
INPUT_DIR="/global/scratch/users/enricocalvane/riboseq/Xu2017/trimmed"
OUTPUT_DIR="/global/scratch/users/enricocalvane/riboseq/Xu2017/rRNA_filtered"
REFERENCE="/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/rRNA_tRNA_index"

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Process each trimmed fastq file
for file in ${INPUT_DIR}/*_trimmed.fq; do
    # Extract base filename
    base=$(basename ${file} _trimmed.fq)
    
    # Bowtie2 alignment to rRNA/tRNA index
    # --very-fast for quick alignment
    # --un-gz will output reads that do NOT align to the reference
    bowtie2 -x ${REFERENCE} \
            -U ${file} \
            --very-fast \
            --un-gz ${OUTPUT_DIR}/${base}_non_rRNA_tRNA.fastq.gz \
            -S ${OUTPUT_DIR}/${base}_rRNA_tRNA_alignments.sam

    # Optional: compress the SAM file to save space
    gzip ${OUTPUT_DIR}/${base}_rRNA_tRNA_alignments.sam
done

echo "rRNA/tRNA filtering completed."
