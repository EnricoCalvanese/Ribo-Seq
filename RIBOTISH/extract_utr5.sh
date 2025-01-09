#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=2:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=derive_utr5
#SBATCH --output=derive_utr5_%j.out
#SBATCH --error=derive_utr5_%j.err

# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/ribotish

# Create directory for reference files if it doesn't exist
mkdir -p reference

# Run the Python script to derive 5' UTRs
python3 derive_utr5.py

# Sort the BED file
module load bio/bedtools2/2.31.0-gcc-11.4.0
sort -k1,1 -k2,2n reference/tair10_5utr.bed > reference/tair10_5utr.sorted.bed

# Show some statistics
echo "First few lines of the sorted BED file:"
head -n 3 reference/tair10_5utr.sorted.bed

echo "Total number of 5' UTR regions derived:"
wc -l reference/tair10_5utr.sorted.bed
