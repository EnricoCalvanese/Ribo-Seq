#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=06:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL

# Directory variables
UNIQUE_READS="/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads"
REF_DIR="/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference"
GTF_FILE="${REF_DIR}/Arabidopsis_thaliana.TAIR10.60.gtf"

cd /global/scratch/users/enricocalvane/riboseq/imb2

# Create output directory for counts
mkdir final_read_count

cd final_read_count

# Function to process each sample
process_sample() {
    local input_bam=$1
    local sample_name=$2
    local description=$3
    
    echo "Processing ${description}..."
    
    # Sort by name (required for htseq-count)
    echo "Sorting ${sample_name} by name..."
    samtools sort -n ${UNIQUE_READS}/${input_bam} -o ${UNIQUE_READS}/${sample_name}_namesort.bam
    
    # Run htseq-count
    echo "Counting ${sample_name}..."
    htseq-count \
        -f bam \
        -r name \
        -s yes \
        -t CDS \
        -i gene_id \
        -m union \
        --nonunique none \
        ${UNIQUE_READS}/${sample_name}_namesort.bam \
        ${GTF_FILE} \
        > $/global/scratch/users/enricocalvane/riboseq/imb2/final_read_count/${sample_name}_htseq_count.txt \
        2> $/global/scratch/users/enricocalvane/riboseq/imb2/final_read_count/${sample_name}_htseq.log
}

# Process RNA-seq samples
echo "Processing RNA-seq samples..."
process_sample "LZT101-1_uniq_sort.bam" "LZT101-1" "WT-RNA-Seq Rep1"
process_sample "LZT101-2_uniq_sort.bam" "LZT101-2" "WT-RNA-Seq Rep2"
process_sample "LZT102-1_uniq_sort.bam" "LZT102-1" "imb2-RNA-Seq Rep1"
process_sample "LZT102-2_uniq_sort.bam" "LZT102-2" "imb2-RNA-Seq Rep2"

# Process Ribo-seq samples
echo "Processing Ribo-seq samples..."
process_sample "LZT103-1_uniq_sort.bam" "LZT103-1" "WT-Ribo-Seq Rep1"
process_sample "LZT103-2_uniq_sort.bam" "LZT103-2" "WT-Ribo-Seq Rep2"
process_sample "LZT104-1_uniq_sort.bam" "LZT104-1" "imb2-Ribo-Seq Rep1"
process_sample "LZT104-2_uniq_sort.bam" "LZT104-2" "imb2-Ribo-Seq Rep2"

echo "All samples processed! Count files are in ${UNIQUE_READS}/counts/"
