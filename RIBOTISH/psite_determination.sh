#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2_bigmem
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=psite_determ
#SBATCH --output=psite_determ_%j.out
#SBATCH --error=psite_determ_%j.err

# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/ribotish

# Create output directory
mkdir -p psite_results

# Array of sample names
samples=(LZT103-1 LZT103-2 LZT104-1 LZT104-2)

# Function to process one sample
process_sample() {
    sample=$1
    echo "Processing P-site determination for ${sample}"
    
    ribo-TISH -b /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/${sample}_uniq_sort.bam \
        -g /global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
        -s /global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf \
        --align-type forward \
        -o psite_results/${sample} \
        --aTIS
}

# Process all samples in parallel
for sample in "${samples[@]}"; do
    process_sample "$sample" &
done

# Wait for all background processes to complete
wait

echo "P-site determination complete for all samples"
