#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=quality_control
#SBATCH --output=quality_control_%j.out
#SBATCH --error=quality_control_%j.err

# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/ribotish

# Create output directory if it doesn't exist
mkdir -p quality_results

# Define paths to required files
GENOME="/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
GTF="/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf"

# Function to process one sample
process_sample() {
    local sample_id=$1
    local bam_path="/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/${sample_id}_uniq_sort.bam"
    local output_prefix="quality_results/${sample_id}"
    
    echo "Processing quality control for ${sample_id}"
    
    # Run ribotish quality with correct parameters
    ribotish quality \
        -b "${bam_path}" \
        -g "${GTF}" \
        -o "${output_prefix}_qual.txt" \
        -f "${output_prefix}_qual.pdf" \
        -r "${output_prefix}.para.py" \
        -l 25,35 \
        -d -40,20 \
        --bins 20 \
        -p 6 \
        -v
        
    echo "Completed processing ${sample_id}"
}

# Array of sample IDs
samples=(
    "LZT103-1"  # WT Rep1
    "LZT103-2"  # WT Rep2
    "LZT104-1"  # imb2 Rep1
    "LZT104-2"  # imb2 Rep2
)

# Process samples in parallel, but not all at once to avoid overloading
# We'll process 4 samples using 6 cores each (24 cores total)
for sample_id in "${samples[@]}"; do
    process_sample "$sample_id" &
done

# Wait for all background processes to complete
wait

echo "Quality control complete for all samples"

# Create a summary report
echo "Creating summary report..."
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
