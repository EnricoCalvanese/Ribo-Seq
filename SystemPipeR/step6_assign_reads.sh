#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=8:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=count_reads
#SBATCH --output=count_reads_%j.out
#SBATCH --error=count_reads_%j.err

# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/systemPipeR/attempt2

# Load required module
module load python

# Function to count reads and verify output
count_reads() {
    local bam=$1
    local gtf=$2
    local output=$3
    local feature_type=$4
    local sample_name=$(basename "$bam" _uniq_sort.bam)
    
    echo "Processing ${feature_type} counts for ${sample_name}..."
    
    # Run htseq-count
    htseq-count \
        -f bam \
        -r name \
        -s yes \
        -t ${feature_type} \
        -i gene_id \
        ${bam} ${gtf} > ${output}
    
    # Verify the output
    if [ ! -s "${output}" ]; then
        echo "Error: No output generated for ${sample_name} ${feature_type}"
        exit 1
    fi
    
    # Count number of features with non-zero counts
    local nonzero=$(awk '$2 > 0' ${output} | wc -l)
    echo "${sample_name} ${feature_type}: Found ${nonzero} features with non-zero counts"
}

# Create a directory for count summaries
mkdir -p results/counts/summaries

# Process each sample
BAM_DIR="/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads"
SAMPLES="LZT103-1 LZT103-2 LZT104-1 LZT104-2"

echo "Starting read counting process..."
echo "Using properly positioned uORF and mORF coordinates"

for sample in ${SAMPLES}; do
    echo "Processing sample ${sample}..."
    
    # Count reads for uORFs
    count_reads \
        ${BAM_DIR}/${sample}_uniq_sort.bam \
        results/counts/predicted_uorfs.gtf \
        results/counts/${sample}_uorf_counts.txt \
        uORF
    
    # Count reads for mORFs
    count_reads \
        ${BAM_DIR}/${sample}_uniq_sort.bam \
        results/counts/morfs.gtf \
        results/counts/${sample}_morf_counts.txt \
        mORF
    
    # Create a summary of the counts
    echo "Creating summary for ${sample}..."
    {
        echo "Count Summary for ${sample}"
        echo "========================="
        echo "uORF Counts:"
        echo "-----------"
        echo "Total features: $(wc -l < results/counts/${sample}_uorf_counts.txt)"
        echo "Features with counts > 0: $(awk '$2 > 0' results/counts/${sample}_uorf_counts.txt | wc -l)"
        echo "Maximum count: $(sort -k2,2n results/counts/${sample}_uorf_counts.txt | tail -n 1)"
        echo
        echo "mORF Counts:"
        echo "-----------"
        echo "Total features: $(wc -l < results/counts/${sample}_morf_counts.txt)"
        echo "Features with counts > 0: $(awk '$2 > 0' results/counts/${sample}_morf_counts.txt | wc -l)"
        echo "Maximum count: $(sort -k2,2n results/counts/${sample}_morf_counts.txt | tail -n 1)"
    } > results/counts/summaries/${sample}_count_summary.txt
done

echo "Counting complete. Check results/counts/summaries/ for detailed statistics."
