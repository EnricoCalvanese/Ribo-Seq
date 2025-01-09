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
mkdir -p quality_results

# Set the distance range as an environment variable
# This helps avoid shell parsing issues with negative numbers
export DIST_RANGE="-40,20"

# Function to process one sample
process_sample() {
    local sample_id=$1
    local bam_path="/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/${sample_id}_uniq_sort.bam"
    local output_prefix="quality_results/${sample_id}"
    
    echo "Processing quality control for ${sample_id}"
    
    # Use the environment variable for the distance parameter
    ribotish quality \
        -b "${bam_path}" \
        -g "/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf" \
        -o "${output_prefix}_qual.txt" \
        -f "${output_prefix}_qual.pdf" \
        -r "${output_prefix}.para.py" \
        -l "25,35" \
        -d "${DIST_RANGE}" \
        -p 6 \
        -v
        
    echo "Completed processing ${sample_id}"
}

# Array of sample IDs with descriptions
samples=(
    "LZT103-1"  # WT Rep1
    "LZT103-2"  # WT Rep2
    "LZT104-1"  # imb2 Rep1
    "LZT104-2"  # imb2 Rep2
)

# Process each sample
for sample_id in "${samples[@]}"; do
    process_sample "$sample_id"
done

echo "Quality control complete for all samples"

# Create a summary report
{
    echo "Quality Control Summary"
    echo "======================"
    echo ""
    
    for sample_id in "${samples[@]}"; do
        echo "Sample: ${sample_id}"
        echo "-------------------------"
        
        qual_file="quality_results/${sample_id}_qual.txt"
        if [ -f "$qual_file" ]; then
            echo "Quality control completed successfully"
            echo "Summary metrics:"
            head -n 5 "$qual_file"
        else
            echo "Quality control failed or incomplete"
        fi
        echo ""
    done
} > quality_results/summary.txt

echo "All processing complete. Check quality_results/summary.txt for results."
