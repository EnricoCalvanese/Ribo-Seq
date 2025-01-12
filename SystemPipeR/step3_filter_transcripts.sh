#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=1:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=filter_ref_utrs
#SBATCH --output=filter_ref_utrs_%j.out
#SBATCH --error=filter_ref_utrs_%j.err

# Stop script if any command fails
set -e

# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/systemPipeR/attempt2

# Create output directories
mkdir -p results/reference_analysis/logs

# Process the reference gene list to get unique genes
echo "Processing reference gene list..."
sort -u results/xu2017sequences.txt > results/reference_analysis/unique_reference_genes.txt
total_ref_genes=$(wc -l < results/reference_analysis/unique_reference_genes.txt)
echo "Found ${total_ref_genes} unique reference genes"

# Create a temporary directory for processing
temp_dir=$(mktemp -d)
trap 'rm -rf "$temp_dir"' EXIT

# Extract gene IDs from our UTR sequences
echo "Processing UTR sequences..."
awk '/^>/ {
    # Extract gene ID from sequence header
    gsub(">", "", $0)
    split($0, parts, "\\.")
    gsub("transcript:", "", parts[1])
    print parts[1]
}' results/sequences/utr_sequences.fa | sort -u > "$temp_dir/our_genes.txt"

total_our_genes=$(wc -l < "$temp_dir/our_genes.txt")
echo "Found ${total_our_genes} unique genes in our UTR set"

# Find intersection of our genes with reference genes
comm -12 "$temp_dir/our_genes.txt" results/reference_analysis/unique_reference_genes.txt > \
    results/reference_analysis/genes_to_keep.txt

genes_to_keep=$(wc -l < results/reference_analysis/genes_to_keep.txt)
echo "Found ${genes_to_keep} genes present in both sets"

# Create a filtered FASTA file containing only UTRs from reference genes
echo "Creating filtered FASTA file..."
awk -v genes_file="results/reference_analysis/genes_to_keep.txt" '
    BEGIN {
        while ((getline gene < genes_file) > 0) {
            genes[gene] = 1
        }
        keep = 0
    }
    /^>/ {
        header = $0
        gsub(">", "", header)
        split(header, parts, "\\.")
        gsub("transcript:", "", parts[1])
        keep = (genes[parts[1]] == 1)
        if (keep) print $0
    }
    !/^>/ {
        if (keep) print $0
    }
' results/sequences/utr_sequences.fa > results/reference_analysis/filtered_utrs.fa

# Calculate sequence statistics
echo "Calculating sequence statistics..."
total_sequences=$(grep -c "^>" results/reference_analysis/filtered_utrs.fa)
echo "Total filtered sequences: ${total_sequences}"

# Calculate length statistics
awk '/^>/ {if (seq) print length(seq); seq=""; next} {seq=seq $0} END {if (seq) print length(seq)}' \
    results/reference_analysis/filtered_utrs.fa > "$temp_dir/lengths.txt"

# Calculate min, max, mean, and median lengths
min_len=$(sort -n "$temp_dir/lengths.txt" | head -n1)
max_len=$(sort -n "$temp_dir/lengths.txt" | tail -n1)
mean_len=$(awk '{ sum += $1 } END { print int(sum/NR) }' "$temp_dir/lengths.txt")
median_len=$(sort -n "$temp_dir/lengths.txt" | awk -v n=$(wc -l < "$temp_dir/lengths.txt") \
    'NR == int((n+1)/2)')

# Create summary report
cat > results/reference_analysis/filtering_summary.txt << EOF
UTR Filtering Summary for Reference Genes
=======================================

Reference Gene Statistics:
- Total unique reference genes: ${total_ref_genes}
- Total genes in our UTR set: ${total_our_genes}
- Genes found in both sets: ${genes_to_keep}
- Coverage of reference set: $(echo "scale=2; 100 * ${genes_to_keep}/${total_ref_genes}" | bc)%

Filtered Sequence Statistics:
- Total UTR sequences retained: ${total_sequences}
- Average sequences per gene: $(echo "scale=2; ${total_sequences}/${genes_to_keep}" | bc)

Length Distribution of Filtered UTRs:
- Minimum length: ${min_len} nt
- Maximum length: ${max_len} nt
- Mean length: ${mean_len} nt
- Median length: ${median_len} nt

Notes:
- Original reference set contains ${total_ref_genes} unique genes out of 9887 total entries
- Each gene may have multiple UTR sequences due to alternative transcripts
- All sequences are from genes identified in the Xu 2017 study
EOF

echo "Filtering complete. Check results/reference_analysis/filtering_summary.txt for details."
