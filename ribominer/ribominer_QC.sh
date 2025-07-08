#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=ribominer_QC.log

# RiboMiner Periodicity and Metaplots Analysis Script
# This script performs periodicity analysis and generates metaplots for ribosome profiling data

# Define variables
WORK_DIR="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer"
BAM_DIR="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/transcriptome_aligned_reads/unique_reads"
RIBOCODE_ANNOT="${WORK_DIR}/prepared_transcripts"
LONGEST_TRANSCRIPTS_INFO="${WORK_DIR}/longest.transcripts.info.txt"
OUTPUT_DIR="${WORK_DIR}/periodicity_analysis"
METAPLOTS_DIR="${WORK_DIR}/metaplots"

# Create output directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$METAPLOTS_DIR"
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



echo "Starting RiboMiner Analysis"
echo "Work directory: $WORK_DIR"
echo "BAM directory: $BAM_DIR"
echo "Periodicity output directory: $OUTPUT_DIR"
echo "Metaplots output directory: $METAPLOTS_DIR"
echo "========================================"

# STEP 2: Generate Metaplots for each sample
echo ""
echo "STEP 2: Generating metaplots for each sample..."
echo "==============================================="

for sample in "${!BAM_FILES[@]}"; do
    bam_file="${BAM_DIR}/${BAM_FILES[$sample]}"
    output_prefix="${METAPLOTS_DIR}/${sample}_metaplots"
    
    echo "Generating metaplots for $sample (${BAM_LEGENDS[$sample]})..."
    
    # Generate metaplots using the correct RiboMiner syntax
    metaplots \
        -a "$RIBOCODE_ANNOT" \
        -r "$bam_file" \
        -o "$output_prefix"
    
    if [[ $? -eq 0 ]]; then
        echo "✓ Metaplots generated for $sample"
    else
        echo "✗ Metaplots generation failed for $sample"
        exit 1
    fi
done

# STEP 3: Create summary report
echo ""
echo "STEP 3: Creating summary report..."
echo "================================="

SUMMARY_FILE="${WORK_DIR}/analysis_summary.txt"
cat > "$SUMMARY_FILE" << EOF
RiboMiner Analysis Summary
Generated on: $(date)
Analysis Directory: $WORK_DIR
BAM Files Directory: $BAM_DIR
Periodicity Output Directory: $OUTPUT_DIR
Metaplots Output Directory: $METAPLOTS_DIR

Samples Analyzed:
EOF

for sample in "${!BAM_FILES[@]}"; do
    echo "- $sample: ${BAM_FILES[$sample]} (${BAM_LEGENDS[$sample]})" >> "$SUMMARY_FILE"
done

cat >> "$SUMMARY_FILE" << EOF

Files Generated:

Periodicity Analysis:
- Individual periodicity plots: ${OUTPUT_DIR}/*_periodicity.pdf
- Use these to determine optimal read lengths and offsets for your data

Metaplots:
- Individual sample metaplots: ${METAPLOTS_DIR}/*_metaplots.pdf
- These plots show periodicity patterns for ribosome footprints of different lengths

Analysis Parameters:
- Read length range for periodicity: 25-35 nt
- Metaplots: Generated using RiboMiner's quality control pipeline

Next Steps:
1. Review periodicity plots to confirm 3-nt periodicity and identify optimal read lengths
2. Review metaplots to assess periodicity patterns across different read lengths
3. Record read lengths and offsets with good periodicity for downstream analysis
4. Use results to filter and process data for further RiboMiner analysis

Quality Check:
- Good ribosome profiling data should show:
  * Strong 3-nt periodicity in periodicity plots
  * Clear ribosome occupancy peak at start codons in metaplots
  * Gradual decrease in occupancy at stop codons in metaplots

Output directories:
- Periodicity: $OUTPUT_DIR
- Metaplots: $METAPLOTS_DIR
EOF

echo ""
echo "Analysis completed successfully!"
echo "Summary report: $SUMMARY_FILE"
echo ""
echo "Output locations:"
echo "- Periodicity analysis: $OUTPUT_DIR"
echo "- Metaplots: $METAPLOTS_DIR"
echo ""
echo "Review the generated plots to assess data quality and ribosome profiling characteristics."
