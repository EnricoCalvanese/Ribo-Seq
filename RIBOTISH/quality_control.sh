#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=1:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=ribotish_test
#SBATCH --output=ribotish_test_%j.out
#SBATCH --error=ribotish_test_%j.err

# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/ribotish

# Create test directory
mkdir -p test_output

# Run ribotish with only the required parameters first
ribotish quality \
    -b /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT103-1_uniq_sort.bam \
    -o test_output/basic_test.txt

echo "Basic test complete"

# If the basic test works, try adding one parameter at a time
if [ $? -eq 0 ]; then
    echo "Testing with gene annotation..."
    ribotish quality \
        -b /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT103-1_uniq_sort.bam \
        -g /global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf \
        -o test_output/gene_test.txt
fi

# Print the version of ribotish we're using
echo "Ribotish version:"
ribotish --version
