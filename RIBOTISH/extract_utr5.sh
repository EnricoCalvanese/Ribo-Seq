#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=2:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=extract_utr5
#SBATCH --output=extract_utr5_%j.out
#SBATCH --error=extract_utr5_%j.err

# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/ribotish

# Load required module
module load bio/bedtools2/2.31.0-gcc-11.4.0

# Create directory for reference files if it doesn't exist
mkdir -p reference

# Extract 5' UTR regions from GTF
# We'll process each field carefully and create a proper BED format
awk '
BEGIN {OFS="\t"}
$3=="five_prime_utr" {
    # Extract gene_id from the 9th field
    match($9, /gene_id "([^"]+)"/, gene)
    # Extract transcript_id from the 9th field
    match($9, /transcript_id "([^"]+)"/, transcript)
    
    # Print in BED format:
    # chr start end name score strand
    print $1,                    # chromosome
          $4-1,                  # start (0-based)
          $5,                    # end
          transcript[1],         # transcript ID as name
          ".",                   # score (using . as placeholder)
          $7                     # strand
}' /global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf > reference/tair10_5utr.bed

# Sort the BED file
sort -k1,1 -k2,2n reference/tair10_5utr.bed > reference/tair10_5utr.sorted.bed

# Let's verify the output
echo "First few lines of the sorted BED file:"
head -n 3 reference/tair10_5utr.sorted.bed

# Count total entries
echo "Total number of 5' UTR regions extracted:"
wc -l reference/tair10_5utr.sorted.bed

echo "5' UTR extraction complete"
