#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
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

# Create output directory if it doesn't exist
mkdir -p psite_results

# Define paths to required files
GENOME="/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
GTF="/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf"
UTR_BED="reference/tair10_5utr.sorted.bed"

# Array of sample information
declare -A samples=(
    ["LZT103-1"]="WT_Rep1"
    ["LZT103-2"]="WT_Rep2"
    ["LZT104-1"]="imb2_Rep1"
    ["LZT104-2"]="imb2_Rep2"
)

# Function to process one sample
process_sample() {
    local sample_id=$1
    local sample_name=${samples[$sample_id]}
    
    echo "Processing quality control for ${sample_id} (${sample_name})"
    
    # Create sample-specific output directory
    mkdir -p psite_results/${sample_name}
    
    # Run ribotish quality control
    ribotish quality \
        -b /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/${sample_id}_uniq_sort.bam \
        -g ${GENOME} \
        -s ${GTF} \
        --utr5 ${UTR_BED} \
        --forward \
        -o psite_results/${sample_name} \
        -p 24
        
    echo "Completed processing ${sample_id}"
}

# Process all samples in parallel
for sample_id in "${!samples[@]}"; do
    process_sample "$sample_id" &
done

# Wait for all background processes to complete
wait

echo "Quality control complete for all samples"

# Create a summary report
echo "Creating summary report..."
echo "Quality Control Summary" > psite_results/summary.txt
echo "======================" >> psite_results/summary.txt
echo "" >> psite_results/summary.txt

for sample_id in "${!samples[@]}"; do
    sample_name=${samples[$sample_id]}
    echo "Sample: ${sample_id} (${sample_name})" >> psite_results/summary.txt
    echo "-------------------------" >> psite_results/summary.txt
    
    if [ -d "psite_results/${sample_name}" ]; then
        echo "Quality control results available" >> psite_results/summary.txt
        # Add any available metrics from the output directory
        if [ -f "psite_results/${sample_name}/quality.txt" ]; then
            echo "Metrics from quality control:" >> psite_results/summary.txt
            cat "psite_results/${sample_name}/quality.txt" >> psite_results/summary.txt
        fi
    else
        echo "Quality control failed or incomplete" >> psite_results/summary.txt
    fi
    echo "" >> psite_results/summary.txt
done
