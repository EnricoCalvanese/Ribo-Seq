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

# Function to process one sample
process_sample() {
    local sample_id=$1
    local bam_path="/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/${sample_id}_uniq_sort.bam"
    local output_prefix="quality_results/${sample_id}"
    
    echo "Processing quality control for ${sample_id}"
    
    # Run ribotish quality with correct parameter formatting
    # Note: -d -40,20 is passed as a single argument with comma separation
    ribotish quality \
        -b "${bam_path}" \
        -g "/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf" \
        -o "${output_prefix}_qual.txt" \
        -f "${output_prefix}_qual.pdf" \
        -r "${output_prefix}.para.py" \
        -l "25,35" \
        -d "-40,20" \
        -p 6 \
        -v
        
    echo "Completed processing ${sample_id}"
}

# Array of sample IDs with descriptions
declare -A samples=(
    ["LZT103-1"]="WT_Rep1"
    ["LZT103-2"]="WT_Rep2"
    ["LZT104-1"]="imb2_Rep1"
    ["LZT104-2"]="imb2_Rep2"
)

# Process samples
for sample_id in "${!samples[@]}"; do
    process_sample "$sample_id" &
done

# Wait for all background processes to complete
wait

echo "Quality control complete for all samples"

# Create a summary report
{
    echo "Quality Control Summary"
    echo "======================"
    echo ""
    
    for sample_id in "${!samples[@]}"; do
        echo "Sample: ${sample_id} (${samples[$sample_id]})"
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
