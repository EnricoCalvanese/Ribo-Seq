#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=03:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=indexing_alignment_hisat.log

# Set variables
REFERENCE_DIR="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/reads"
REFERENCE_FA="/global/scratch/users/enricocalvane/riboseq/Athaliana_447_TAIR10.fa"
INDEX_PREFIX="${REFERENCE_DIR}/Athaliana_447_TAIR10_index"
OUTPUT_DIR="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/reads/hisat_alignment"

cd /global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/reads

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

module load bio/hisat2/2.2.1-gcc-11.4.0

echo "Step 1: Building HISAT2 index"
hisat2-build "$REFERENCE_FA" "$INDEX_PREFIX" 2> "${REFERENCE_DIR}/log2_hisat2_build_log.txt"
if [ $? -ne 0 ]; then
    echo "Error: Failed to build HISAT2 index"
    exit 1
fi

echo "Step 2: Aligning sequencing data"

# Align each pair of files explicitly
hisat2 -p 24 -x "$INDEX_PREFIX" -1 /global/scratch/users/enricocalvane/riboseq/imb2/filtered_reads/LZT101-1.filtered.1.fq \
    -2 /global/scratch/users/enricocalvane/riboseq/imb2/filtered_reads/LZT101-1.filtered.2.fq \
    > "${OUTPUT_DIR}/LZT101-1.sam" \
    2> "${OUTPUT_DIR}/log_LZT101-1.txt"

hisat2 -p 24 -x "$INDEX_PREFIX" -1 /global/scratch/users/enricocalvane/riboseq/imb2/filtered_reads/LZT101-2.filtered.1.fq \
    -2 /global/scratch/users/enricocalvane/riboseq/imb2/filtered_reads/LZT101-2.filtered.2.fq \
    > "${OUTPUT_DIR}/LZT101-2.sam" \
    2> "${OUTPUT_DIR}/log_LZT101-2.txt"

hisat2 -p 24 -x "$INDEX_PREFIX" -1 /global/scratch/users/enricocalvane/riboseq/imb2/filtered_reads/LZT102-1.filtered.1.fq \
    -2 /global/scratch/users/enricocalvane/riboseq/imb2/filtered_reads/LZT102-1.filtered.2.fq \
    > "${OUTPUT_DIR}/LZT102-1.sam" \
    2> "${OUTPUT_DIR}/log_LZT102-1.txt"

hisat2 -p 24 -x "$INDEX_PREFIX" -1 /global/scratch/users/enricocalvane/riboseq/imb2/filtered_reads/LZT102-2.filtered.1.fq \
    -2 /global/scratch/users/enricocalvane/riboseq/imb2/filtered_reads/LZT102-2.filtered.2.fq \
    > "${OUTPUT_DIR}/LZT102-2.sam" \
    2> "${OUTPUT_DIR}/log_LZT102-2.txt"

hisat2 -p 24 -x "$INDEX_PREFIX" -1 /global/scratch/users/enricocalvane/riboseq/imb2/filtered_reads/LZT103-1.filtered.1.fq \
    -2 /global/scratch/users/enricocalvane/riboseq/imb2/filtered_reads/LZT103-1.filtered.2.fq \
    > "${OUTPUT_DIR}/LZT103-1.sam" \
    2> "${OUTPUT_DIR}/log_LZT103-1.txt"

hisat2 -p 24 -x "$INDEX_PREFIX" -1 /global/scratch/users/enricocalvane/riboseq/imb2/filtered_reads/LZT103-2.filtered.1.fq \
    -2 /global/scratch/users/enricocalvane/riboseq/imb2/filtered_reads/LZT103-2.filtered.2.fq \
    > "${OUTPUT_DIR}/LZT103-2.sam" \
    2> "${OUTPUT_DIR}/log_LZT103-2.txt"

hisat2 -p 24 -x "$INDEX_PREFIX" -1 /global/scratch/users/enricocalvane/riboseq/imb2/filtered_reads/LZT104-1.filtered.1.fq \
    -2 /global/scratch/users/enricocalvane/riboseq/imb2/filtered_reads/LZT104-1.filtered.2.fq \
    > "${OUTPUT_DIR}/LZT104-1.sam" \
    2> "${OUTPUT_DIR}/log_LZT104-1.txt"

hisat2 -p 24 -x "$INDEX_PREFIX" -1 /global/scratch/users/enricocalvane/riboseq/imb2/filtered_reads/LZT104-2.filtered.1.fq \
    -2 /global/scratch/users/enricocalvane/riboseq/imb2/filtered_reads/LZT104-2.filtered.2.fq \
    > "${OUTPUT_DIR}/LZT104-2.sam" \
    2> "${OUTPUT_DIR}/log_LZT104-2.txt"

echo "All samples processed successfully. SAM files are in $OUTPUT_DIR."
