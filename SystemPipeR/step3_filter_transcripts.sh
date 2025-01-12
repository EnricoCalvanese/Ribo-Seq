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

# First, identify all available transcripts and their genes
echo "Cataloging all available transcripts..."
awk '/^>/ {
    header = $0
    gsub(">", "", header)
    # Split into gene and transcript version
    split(header, parts, "\\.")
    gsub("transcript:", "", parts[1])
    # Store full header and extract gene ID
    print parts[1] "\t" header
}' results/sequences/utr_sequences.fa | sort -k1,1 > "$temp_dir/all_transcripts.txt"

# Find primary (.1) transcripts where available
echo "Identifying primary transcripts..."
awk '{
    gene = $1
    transcript = $2
    if (transcript ~ /\.1$/) {
        primary[gene] = transcript
    } else if (!(gene in primary)) {
        # Store first transcript as backup if no .1 exists
        if (!(gene in backup)) {
            backup[gene] = transcript
        }
    }
}
END {
    # Output primary or backup transcript for each gene
    for (gene in primary) {
        print gene "\t" primary[gene] "\t" "primary"
    }
    for (gene in backup) {
        if (!(gene in primary)) {
            print gene "\t" backup[gene] "\t" "backup"
        }
    }
}' "$temp_dir/all_transcripts.txt" > "$temp_dir/representative_transcripts.txt"

# Identify missing reference genes
echo "Identifying missing reference genes..."
comm -23 results/reference_analysis/unique_reference_genes.txt \
    <(cut -f1 "$temp_dir/representative_transcripts.txt" | sort) > "$temp_dir/missing_genes.txt"

missing_count=$(wc -l < "$temp_dir/missing_genes.txt")
echo "Found ${missing_count} missing reference genes"

# Create a list of transcripts to keep
echo "Creating filtered transcript list..."
while IFS=$'\t' read -r gene transcript type; do
    if grep -q "^${gene}$" results/reference_analysis/unique_reference_genes.txt; then
        echo "${transcript}"
    fi
done < "$temp_dir/representative_transcripts.txt" > "$temp_dir/transcripts_to_keep.txt"

# Create filtered FASTA file
echo "Creating filtered FASTA file..."
awk -v transcripts_file="$temp_dir/transcripts_to_keep.txt" '
    BEGIN {
        while ((getline transcript < transcripts_file) > 0) {
            keep_transcript[">" transcript] = 1
        }
        keep = 0
    }
    /^>/ {
        keep = ($0 in keep_transcript)
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
awk '/^>/ {
    if (seq) {
        print length(seq)
        seq=""
    }
    next
}
{
    seq = seq $0
}
END {
    if (seq) print length(seq)
}' results/reference_analysis/filtered_utrs.fa > "$temp_dir/lengths.txt"

min_len=$(sort -n "$temp_dir/lengths.txt" | head -n1)
max_len=$(sort -n "$temp_dir/lengths.txt" | tail -n1)
mean_len=$(awk '{ sum += $1 } END { print int(sum/NR) }' "$temp_dir/lengths.txt")
median_len=$(sort -n "$temp_dir/lengths.txt" | awk -v n=$(wc -l < "$temp_dir/lengths.txt") \
    'NR == int((n+1)/2)')

# Save list of missing genes with details
echo "Gene ID" > results/reference_analysis/missing_genes_details.txt
cat "$temp_dir/missing_genes.txt" >> results/reference_analysis/missing_genes_details.txt

# Create comprehensive summary report
cat > results/reference_analysis/filtering_summary.txt << EOF
UTR Filtering Summary for Reference Genes
=======================================

Reference Gene Statistics:
- Total unique reference genes: ${total_ref_genes}
- Genes with representative transcripts: $(( total_ref_genes - missing_count ))
- Missing genes: ${missing_count}

Transcript Selection:
- Primary (.1) transcripts used where available
- Backup transcripts used when primary not available
- Total sequences selected: ${total_sequences}

Length Distribution of Selected UTRs:
- Minimum length: ${min_len} nt
- Maximum length: ${max_len} nt
- Mean length: ${mean_len} nt
- Median length: ${median_len} nt

Missing Genes:
- List saved in 'missing_genes_details.txt'
- These genes require manual investigation
- Possible reasons for missing genes:
  * Gene not present in original annotation
  * No 5' UTR annotation available
  * UTR filtered out in previous steps

Next Steps:
- Investigate missing genes in original annotation
- Consider retrieving UTRs for missing genes from genome
- Proceed with uORF prediction on available sequences
EOF

echo "Filtering complete. Check results/reference_analysis/ for detailed reports."
