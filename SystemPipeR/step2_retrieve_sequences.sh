#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=8:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=utr_sequences
#SBATCH --output=utr_sequences_%j.out
#SBATCH --error=utr_sequences_%j.err

# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/systemPipeR/attempt2

# Load required modules
module load bio/bedtools2/2.31.0-gcc-11.4.0

# Define input files
GENOME="/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
UTR_BED="/global/scratch/users/enricocalvane/riboseq/imb2/ribotish/reference/tair10_5utr.sorted.bed"
FILTERED_GENES="results/filtered_genes.txt"

# Create necessary directories
mkdir -p {results/sequences,temp}

# Step 1: Filter the BED file to only include UTRs from our filtered genes
echo "Filtering BED file for selected genes..."
while read gene; do
    grep -P "transcript:${gene}\\." "$UTR_BED"
done < "$FILTERED_GENES" > temp/filtered_utrs.bed

# Sort the filtered BED file
sort -k1,1 -k2,2n temp/filtered_utrs.bed > temp/filtered_utrs.sorted.bed

# Step 2: Extract sequences using bedtools
echo "Extracting sequences..."
bedtools getfasta -fi "$GENOME" \
    -bed temp/filtered_utrs.sorted.bed \
    -s \
    -name \
    -fo results/sequences/utr_sequences.fa

# Step 3: Generate sequence statistics
echo "Generating sequence statistics..."
{
    echo "5' UTR Sequence Retrieval Summary"
    echo "================================"
    echo
    echo "Input Statistics:"
    echo "----------------"
    echo "Number of filtered genes: $(wc -l < $FILTERED_GENES)"
    echo "Number of UTR regions extracted: $(wc -l < temp/filtered_utrs.sorted.bed)"
    echo
    echo "Sequence Statistics:"
    echo "-------------------"
    echo "Number of sequences retrieved: $(grep -c "^>" results/sequences/utr_sequences.fa)"
    echo
    echo "Length Distribution:"
    echo "-------------------"
    awk '/^>/ {next} {print length}' results/sequences/utr_sequences.fa | \
        sort -n | \
        awk '
            BEGIN {
                min=999999; max=0; sum=0; count=0
            }
            {
                sum+=$1; count++
                if($1<min) min=$1
                if($1>max) max=$1
                lengths[count]=$1
            }
            END {
                mean=sum/count
                if(count%2==0) median=(lengths[count/2]+lengths[count/2+1])/2
                else median=lengths[int(count/2)+1]
                print "Minimum length:", min, "nt"
                print "Maximum length:", max, "nt"
                print "Mean length:", int(mean), "nt"
                print "Median length:", int(median), "nt"
            }'
} > results/sequences/sequence_stats.txt

# Create FASTA index for future use
samtools faidx results/sequences/utr_sequences.fa

echo "Analysis complete. Check results/sequences/ for output files."

# Clean up temporary files
rm -rf temp/filtered_utrs.bed temp/filtered_utrs.sorted.bed
