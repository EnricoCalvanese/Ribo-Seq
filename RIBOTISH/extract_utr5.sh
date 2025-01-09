#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2_bigmem
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
awk '$3=="five_prime_utr" {print $1"\t"$4-1"\t"$5"\t"$9"\t.\t"$7}' \
    /global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf \
    | sed 's/gene_id "//;s/"; transcript_id "/\t/;s/";//' \
    > reference/tair10_5utr.bed

# Sort the BED file
sort -k1,1 -k2,2n reference/tair10_5utr.bed > reference/tair10_5utr.sorted.bed

echo "5' UTR extraction complete"
