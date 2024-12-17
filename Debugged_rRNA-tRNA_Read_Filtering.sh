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
    
    # Verbose debugging output
    echo "Processing file: ${file}"
    echo "Reference index: ${REFERENCE}"
    
    # Check if input file exists and is readable
    if [ ! -f "${file}" ]; then
        echo "Error: Input file ${file} does not exist!"
        continue
    fi
    
    # Try alignment with verbose output
    bowtie2 -x ${REFERENCE} \
            -U ${file} \
            -v 3 \
            --un-gz ${OUTPUT_DIR}/${base}_non_rRNA_tRNA.fastq.gz \
            -S ${OUTPUT_DIR}/${base}_rRNA_tRNA_alignments.sam \
            --verbose
    
    # Check the exit status
    if [ $? -ne 0 ]; then
        echo "Alignment failed for ${file}"
    fi
done

echo "rRNA/tRNA filtering completed."
