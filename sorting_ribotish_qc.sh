#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=06:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL

# Directory variables for cleaner code
UNIQUE_READS="/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads"
REF_DIR="/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference"

cd /global/scratch/users/enricocalvane/riboseq/imb2

cd unique_reads

module load bio/samtools/1.17-gcc-11.4.0

# Sort only Ribo-seq BAM files
echo "Sorting Ribo-seq BAM files..."
samtools sort ${UNIQUE_READS}/LZT103-1_uniq.bam -o ${UNIQUE_READS}/LZT103-1_uniq_sort.bam
samtools sort ${UNIQUE_READS}/LZT103-2_uniq.bam -o ${UNIQUE_READS}/LZT103-2_uniq_sort.bam
samtools sort ${UNIQUE_READS}/LZT104-1_uniq.bam -o ${UNIQUE_READS}/LZT104-1_uniq_sort.bam
samtools sort ${UNIQUE_READS}/LZT104-2_uniq.bam -o ${UNIQUE_READS}/LZT104-2_uniq_sort.bam

# Index sorted Ribo-seq BAM files
echo "Indexing sorted Ribo-seq BAM files..."
samtools index ${UNIQUE_READS}/LZT103-1_uniq_sort.bam
samtools index ${UNIQUE_READS}/LZT103-2_uniq_sort.bam
samtools index ${UNIQUE_READS}/LZT104-1_uniq_sort.bam
samtools index ${UNIQUE_READS}/LZT104-2_uniq_sort.bam

# Convert GFF3 to GTF
echo "Converting GFF3 to GTF..."
gffread ${REF_DIR}/Arabidopsis_thaliana.TAIR10.60.gff3 -T -o ${REF_DIR}/Arabidopsis_thaliana.TAIR10.60.gtf

# Run Ribotish quality for Ribo-seq samples only
echo "Running Ribotish quality control on Ribo-seq samples..."
echo "Processing WT Ribo-seq Rep1 (LZT103-1)..."
ribotish quality -b ${UNIQUE_READS}/LZT103-1_uniq_sort.bam -g ${REF_DIR}/Arabidopsis_thaliana.TAIR10.60.gtf > RF_WT_Rep1_quality.txt

echo "Processing WT Ribo-seq Rep2 (LZT103-2)..."
ribotish quality -b ${UNIQUE_READS}/LZT103-2_uniq_sort.bam -g ${REF_DIR}/Arabidopsis_thaliana.TAIR10.60.gtf > RF_WT_Rep2_quality.txt

echo "Processing imb2 Ribo-seq Rep1 (LZT104-1)..."
ribotish quality -b ${UNIQUE_READS}/LZT104-1_uniq_sort.bam -g ${REF_DIR}/Arabidopsis_thaliana.TAIR10.60.gtf > RF_imb2_Rep1_quality.txt

echo "Processing imb2 Ribo-seq Rep2 (LZT104-2)..."
ribotish quality -b ${UNIQUE_READS}/LZT104-2_uniq_sort.bam -g ${REF_DIR}/Arabidopsis_thaliana.TAIR10.60.gtf > RF_imb2_Rep2_quality.txt

echo "Quality control analysis complete!"
