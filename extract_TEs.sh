#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=06:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL

module load bio/bedtools2/2.31.0-gcc-11.4.0
module load bio/samtools/1.17-gcc-11.4.0

# Create directories for filtered files
mkdir -p filtered_bams/TEup filtered_bams/TEdown

# For TEup genes
for bam in LZT103-1_uniq_sort.bam LZT103-2_uniq_sort.bam LZT104-1_uniq_sort.bam LZT104-2_uniq_sort.bam; do
    echo "Processing $bam for TEup genes..."
    output="filtered_bams/TEup/${bam/.bam/_TEup.bam}"
    samtools view -b -L TEup_genes.bed $bam > $output
    samtools index $output
done

# For TEdown genes
for bam in LZT103-1_uniq_sort.bam LZT103-2_uniq_sort.bam LZT104-1_uniq_sort.bam LZT104-2_uniq_sort.bam; do
    echo "Processing $bam for TEdown genes..."
    output="filtered_bams/TEdown/${bam/.bam/_TEdown.bam}"
    samtools view -b -L TEdown_genes.bed $bam > $output
    samtools index $output
done

# Verify the files
echo -e "\nVerifying read counts in original and filtered files:"
echo -e "\nOriginal files:"
for bam in LZT103-1_uniq_sort.bam LZT103-2_uniq_sort.bam LZT104-1_uniq_sort.bam LZT104-2_uniq_sort.bam; do
    count=$(samtools view -c $bam)
    echo "$bam: $count reads"
done

echo -e "\nTEup filtered files:"
for bam in filtered_bams/TEup/*_TEup.bam; do
    count=$(samtools view -c $bam)
    basename=$(basename $bam)
    echo "$basename: $count reads"
done

echo -e "\nTEdown filtered files:"
for bam in filtered_bams/TEdown/*_TEdown.bam; do
    count=$(samtools view -c $bam)
    basename=$(basename $bam)
    echo "$basename: $count reads"
done
