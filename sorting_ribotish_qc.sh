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

# First, sort all BAM files
echo "Sorting BAM files..."
samtools sort ${UNIQUE_READS}/LZT101-1_uniq.bam -o ${UNIQUE_READS}/LZT101-1_uniq_sort.bam
samtools sort ${UNIQUE_READS}/LZT101-2_uniq.bam -o ${UNIQUE_READS}/LZT101-2_uniq_sort.bam
samtools sort ${UNIQUE_READS}/LZT102-1_uniq.bam -o ${UNIQUE_READS}/LZT102-1_uniq_sort.bam
samtools sort ${UNIQUE_READS}/LZT102-2_uniq.bam -o ${UNIQUE_READS}/LZT102-2_uniq_sort.bam
samtools sort ${UNIQUE_READS}/LZT103-1_uniq.bam -o ${UNIQUE_READS}/LZT103-1_uniq_sort.bam
samtools sort ${UNIQUE_READS}/LZT103-2_uniq.bam -o ${UNIQUE_READS}/LZT103-2_uniq_sort.bam
samtools sort ${UNIQUE_READS}/LZT104-1_uniq.bam -o ${UNIQUE_READS}/LZT104-1_uniq_sort.bam
samtools sort ${UNIQUE_READS}/LZT104-2_uniq.bam -o ${UNIQUE_READS}/LZT104-2_uniq_sort.bam

# Index all sorted BAM files
echo "Indexing sorted BAM files..."
samtools index ${UNIQUE_READS}/LZT101-1_uniq_sort.bam
samtools index ${UNIQUE_READS}/LZT101-2_uniq_sort.bam
samtools index ${UNIQUE_READS}/LZT102-1_uniq_sort.bam
samtools index ${UNIQUE_READS}/LZT102-2_uniq_sort.bam
samtools index ${UNIQUE_READS}/LZT103-1_uniq_sort.bam
samtools index ${UNIQUE_READS}/LZT103-2_uniq_sort.bam
samtools index ${UNIQUE_READS}/LZT104-1_uniq_sort.bam
samtools index ${UNIQUE_READS}/LZT104-2_uniq_sort.bam

# Convert GFF3 to GTF
echo "Converting GFF3 to GTF..."
gffread ${REF_DIR}/Arabidopsis_thaliana.TAIR10.60.gff3 -T -o ${REF_DIR}/Arabidopsis_thaliana.TAIR10.60.gtf

# Run Ribotish quality for each sorted BAM file
echo "Running Ribotish quality control..."
ribotish quality -b ${UNIQUE_READS}/LZT101-1_uniq_sort.bam -g ${REF_DIR}/Arabidopsis_thaliana.TAIR10.60.gtf > LZT101-1_quality.txt
ribotish quality -b ${UNIQUE_READS}/LZT101-2_uniq_sort.bam -g ${REF_DIR}/Arabidopsis_thaliana.TAIR10.60.gtf > LZT101-2_quality.txt
ribotish quality -b ${UNIQUE_READS}/LZT102-1_uniq_sort.bam -g ${REF_DIR}/Arabidopsis_thaliana.TAIR10.60.gtf > LZT102-1_quality.txt
ribotish quality -b ${UNIQUE_READS}/LZT102-2_uniq_sort.bam -g ${REF_DIR}/Arabidopsis_thaliana.TAIR10.60.gtf > LZT102-2_quality.txt
ribotish quality -b ${UNIQUE_READS}/LZT103-1_uniq_sort.bam -g ${REF_DIR}/Arabidopsis_thaliana.TAIR10.60.gtf > LZT103-1_quality.txt
ribotish quality -b ${UNIQUE_READS}/LZT103-2_uniq_sort.bam -g ${REF_DIR}/Arabidopsis_thaliana.TAIR10.60.gtf > LZT103-2_quality.txt
ribotish quality -b ${UNIQUE_READS}/LZT104-1_uniq_sort.bam -g ${REF_DIR}/Arabidopsis_thaliana.TAIR10.60.gtf > LZT104-1_quality.txt
ribotish quality -b ${UNIQUE_READS}/LZT104-2_uniq_sort.bam -g ${REF_DIR}/Arabidopsis_thaliana.TAIR10.60.gtf > LZT104-2_quality.txt
