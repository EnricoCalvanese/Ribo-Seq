#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2_bigmem
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=8:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=utr_filter
#SBATCH --output=utr_filter_%j.out
#SBATCH --error=utr_filter_%j.err

# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/systemPipeR/attempt2

# Create necessary directories
mkdir -p {logs,results,temp}

# Load required modules
module load bio/bedtools2/2.31.0-gcc-11.4.0

# Define input files
GTF="/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf"
UTR_BED="/global/scratch/users/enricocalvane/riboseq/imb2/ribotish/reference/tair10_5utr.sorted.bed"

# Step 1: Count total number of unique protein-coding genes
echo "Counting total protein-coding genes..."
total_genes=$(awk '$3=="transcript" {gsub("\"", ""); split($10,a,";"); print a[1]}' "$GTF" | \
              sed 's/gene_id gene://' | sort -u | wc -l)
echo "Total protein-coding genes: $total_genes"

# Step 2: Process 5' UTR BED file to get UTR lengths
echo "Processing 5' UTR lengths..."
awk -v OFS="\t" '{
    len = $3 - $2;
    split($4, id, ":");  # Split on colon to get transcript ID
    print id[2], len;    # Print transcript ID and length
}' "$UTR_BED" > temp/utr_lengths.txt

# Step 3: Filter for UTRs ≥ 15 nt and get unique genes
echo "Filtering UTRs >= 15 nt..."
awk '$2 >= 15 {print $1}' temp/utr_lengths.txt | \
    sed 's/\.[0-9]*$//' | sort -u > temp/genes_with_long_utrs.txt

# Count genes with sufficient UTR length
filtered_genes=$(wc -l < temp/genes_with_long_utrs.txt)

# Create a summary report
{
    echo "UTR Filtering Analysis Summary"
    echo "============================="
    echo
    echo "Initial number of protein-coding genes: $total_genes"
    echo "Genes with 5' UTR >= 15 nt: $filtered_genes"
    echo "Filtered out genes: $((total_genes - filtered_genes))"
    echo
    echo "Percentage retained: $(awk "BEGIN {printf \"%.2f\", ($filtered_genes/$total_genes)*100}")%"
} > results/utr_filtering_summary.txt

echo "Analysis complete. Check results/utr_filtering_summary.txt for summary."

# Clean up temporary files
rm -rf temp/utr_lengths.txt

# Save filtered gene list for next steps
mv temp/genes_with_long_utrs.txt results/filtered_genes.txt
