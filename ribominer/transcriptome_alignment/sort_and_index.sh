#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=sort_and_index.log

# Directory variables for cleaner code
UNIQUE_READS="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/transcriptome_aligned_reads/unique_reads"

cd ${UNIQUE_READS}

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
