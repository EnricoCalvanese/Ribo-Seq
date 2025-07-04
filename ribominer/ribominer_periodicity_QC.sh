#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=ribominer_periodicity.log

# RiboMiner Periodicity Analysis Script
# This script performs periodicity analysis on ribosome profiling data

# Define variables
WORK_DIR="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer"
BAM_DIR="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/reads/unique_reads"
RIBOCODE_ANNOT="${WORK_DIR}/prepared_transcripts"
LONGEST_TRANSCRIPTS_INFO="${WORK_DIR}/longest.transcripts.info.txt"
OUTPUT_DIR="${WORK_DIR}/periodicity_analysis"

# Create output directory
mkdir -p "$OUTPUT_DIR"
cd "$WORK_DIR"

# Define BAM files and their information
declare -A BAM_FILES=(
    ["LZT103-1"]="LZT103-1_uniq_sort.bam"
    ["LZT103-2"]="LZT103-2_uniq_sort.bam"
    ["LZT104-1"]="LZT104-1_uniq_sort.bam"
    ["LZT104-2"]="LZT104-2_uniq_sort.bam"
)

declare -A BAM_LEGENDS=(
    ["LZT103-1"]="WT-Ribo-Seq-Rep1"
    ["LZT103-2"]="WT-Ribo-Seq-Rep2"
    ["LZT104-1"]="imb2-Ribo-Seq-Rep1"
    ["LZT104-2"]="imb2-Ribo-Seq-Rep2"
)

echo "Starting RiboMiner Periodicity Analysis"
echo "Work directory: $WORK_DIR"
echo "BAM directory: $BAM_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "========================================"

# Periodicity Analysis for each sample
echo ""
echo "Analyzing periodicity for each sample..."
echo "======================================="

for sample in "${!BAM_FILES[@]}"; do
    bam_file="${BAM_DIR}/${BAM_FILES[$sample]}"
    output_prefix="${OUTPUT_DIR}/${sample}_periodicity"
    
    echo "Processing $sample (${BAM_LEGENDS[$sample]})..."
    
    Periodicity \
        -i "$bam_file" \
        -a "$RIBOCODE_ANNOT" \
        -o "$output_prefix" \
        -c "$LONGEST_TRANSCRIPTS_INFO" \
        -L 25 \
        -R 35
    
    if [[ $? -eq 0 ]]; then
        echo "✓ Periodicity analysis completed for $sample"
    else
        echo "✗ Periodicity analysis failed for $sample"
        exit 1
    fi
done

# Create summary report
echo ""
echo "Creating summary report..."
echo "========================="

SUMMARY_FILE="${OUTPUT_DIR}/periodicity_summary.txt"
cat > "$SUMMARY_FILE" << EOF
RiboMiner Periodicity Analysis Summary
Generated on: $(date)

Analysis Directory: $WORK_DIR
BAM Files Directory: $BAM_DIR
Output Directory: $OUTPUT_DIR

Samples Analyzed:
EOF

for sample in "${!BAM_FILES[@]}"; do
    echo "- $sample: ${BAM_FILES[$sample]} (${BAM_LEGENDS[$sample]})" >> "$SUMMARY_FILE"
done

cat >> "$SUMMARY_FILE" << EOF

Files Generated:
- Periodicity Analysis: *_periodicity.pdf

Next Steps:
1. Review periodicity plots (*_periodicity.pdf) to determine optimal read lengths and offsets
2. Use the periodicity results to create an attributes.txt file for downstream analysis
3. Good ribosome profiling data should show strong 3-nt periodicity

Output files are located in: $OUTPUT_DIR
EOF

echo ""
echo "Periodicity analysis completed!"
echo "Summary report: $SUMMARY_FILE"
echo "Output files are located in: $OUTPUT_DIR"
