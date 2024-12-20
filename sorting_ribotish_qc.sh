#!/bin/bash

# First, sort all BAM files
samtools sort /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT101-1_uniq.bam -o /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT101-1_uniq_sort.bam
samtools sort /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT101-2_uniq.bam -o /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT101-2_uniq_sort.bam
samtools sort /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT102-1_uniq.bam -o /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT102-1_uniq_sort.bam
samtools sort /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT102-2_uniq.bam -o /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT102-2_uniq_sort.bam
samtools sort /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT103-1_uniq.bam -o /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT103-1_uniq_sort.bam
samtools sort /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT103-2_uniq.bam -o /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT103-2_uniq_sort.bam
samtools sort /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT104-1_uniq.bam -o /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT104-1_uniq_sort.bam
samtools sort /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT104-2_uniq.bam -o /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT104-2_uniq_sort.bam

# Convert GFF3 to GTF
gffread /global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gff3 -T -o /global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf

# Run Ribotish quality for each sorted BAM file
ribotish quality -b /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT101-1_uniq_sort.bam -g /global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf > LZT101-1_quality.txt
ribotish quality -b /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT101-2_uniq_sort.bam -g /global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf > LZT101-2_quality.txt
ribotish quality -b /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT102-1_uniq_sort.bam -g /global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf > LZT102-1_quality.txt
ribotish quality -b /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT102-2_uniq_sort.bam -g /global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf > LZT102-2_quality.txt
ribotish quality -b /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT103-1_uniq_sort.bam -g /global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf > LZT103-1_quality.txt
ribotish quality -b /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT103-2_uniq_sort.bam -g /global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf > LZT103-2_quality.txt
ribotish quality -b /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT104-1_uniq_sort.bam -g /global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf > LZT104-1_quality.txt
ribotish quality -b /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT104-2_uniq_sort.bam -g /global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf > LZT104-2_quality.txt
