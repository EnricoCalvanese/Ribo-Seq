#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=02:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=ribominer_quality_control.log

# RiboMiner Quality Control Script
# This script performs comprehensive quality control on ribosome profiling data
# including periodicity checking, frame distribution, and length distribution

# Define variables
SIF="ribocode_ribominer_latest.sif"
WORK_DIR="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer"
BAM_DIR="/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads"
RIBOCODE_ANNOT="${WORK_DIR}/prepared_transcripts"
LONGEST_TRANSCRIPTS_INFO="${WORK_DIR}/longest.transcripts.info.txt"
QC_OUTPUT_DIR="${WORK_DIR}/quality_control"

# Create output directory
mkdir -p "$QC_OUTPUT_DIR"
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

echo "Starting RiboMiner Quality Control Analysis"
echo "Work directory: $WORK_DIR"
echo "BAM directory: $BAM_DIR"
echo "Output directory: $QC_OUTPUT_DIR"
echo "=========================================="

# Check if BAM files exist and are indexed
echo "Checking BAM files..."
for sample in "${!BAM_FILES[@]}"; do
    bam_file="${BAM_DIR}/${BAM_FILES[$sample]}"
    if [[ ! -f "$bam_file" ]]; then
        echo "ERROR: BAM file not found: $bam_file"
        exit 1
    fi
    
    # Check if index exists
    if [[ ! -f "${bam_file}.bai" ]]; then
        echo "WARNING: Index not found for $bam_file"
        echo "Creating index..."
        singularity exec "$SIF" samtools index "$bam_file"
    fi
    echo "✓ $bam_file (indexed)"
done

# Step 1: Periodicity Analysis for each sample
echo ""
echo "Step 1: Analyzing periodicity for each sample..."
echo "================================================"

echo "Applying RiboMiner bug fix..."

# Run periodicity analysis with the fix applied inline
for sample in "${!BAM_FILES[@]}"; do
    bam_file="${BAM_DIR}/${BAM_FILES[$sample]}"
    output_prefix="${QC_OUTPUT_DIR}/${sample}_periodicity"
    
    echo "Processing $sample..."
    
    # Use sed to fix the bug on-the-fly and run the analysis
    singularity exec --writable-tmpfs "$SIF" /bin/bash -c "
        # Fix the Python script
        sed -i 's/transcript_dict\.key()/transcript_dict.keys()/g' /root/miniconda3/lib/python3.7/site-packages/RiboMiner/Periodicity.py
        
        # Run the analysis
        /root/miniconda3/bin/Periodicity \
            -i '$bam_file' \
            -a '$RIBOCODE_ANNOT' \
            -o '$output_prefix' \
            -c '$LONGEST_TRANSCRIPTS_INFO' \
            -L 25 \
            -R 35
    "
    
    if [[ $? -eq 0 ]]; then
        echo "✓ Periodicity analysis completed for $sample"
    else
        echo "✗ Periodicity analysis failed for $sample"
    fi
done

# Step 2: Frame Distribution Analysis
echo ""
echo "Step 2: Analyzing frame distribution..."
echo "======================================"

# Note: This step requires attributes.txt file to be created first
# We'll create a temporary attributes.txt with estimated parameters
# Users will need to refine this based on periodicity results

echo "Creating temporary attributes.txt file..."
echo "NOTE: You will need to update this file based on periodicity results"

TEMP_ATTRIBUTES="${QC_OUTPUT_DIR}/temp_attributes.txt"
echo -e "bamFiles\treadLengths\tOffsets\tbamLegends" > "$TEMP_ATTRIBUTES"

for sample in "${!BAM_FILES[@]}"; do
    bam_file="${BAM_DIR}/${BAM_FILES[$sample]}"
    # Using common ribosome profiling parameters - adjust based on periodicity results
    read_lengths="27,28,29,30"
    offsets="12,12,13,14"  # Adjust these based on your periodicity analysis
    legend="${BAM_LEGENDS[$sample]}"
    
    echo -e "${bam_file}\t${read_lengths}\t${offsets}\t${legend}" >> "$TEMP_ATTRIBUTES"
done

echo "Temporary attributes.txt created at: $TEMP_ATTRIBUTES"
echo "Please review periodicity results and update read lengths and offsets accordingly"

# Run frame distribution analysis with temporary attributes
echo "Running frame distribution analysis..."
singularity exec "$SIF" \
    /root/miniconda3/bin/RiboDensityOfDiffFrames \
    -f "$TEMP_ATTRIBUTES" \
    -c "$LONGEST_TRANSCRIPTS_INFO" \
    -o "${QC_OUTPUT_DIR}/frame_distribution" \
    --plot yes

if [[ $? -eq 0 ]]; then
    echo "✓ Frame distribution analysis completed"
else
    echo "✗ Frame distribution analysis failed"
fi

# Step 3: Length Distribution Analysis
echo ""
echo "Step 3: Analyzing length distribution..."
echo "======================================="

for sample in "${!BAM_FILES[@]}"; do
    bam_file="${BAM_DIR}/${BAM_FILES[$sample]}"
    output_prefix="${QC_OUTPUT_DIR}/${sample}_length_dist"
    
    echo "Processing length distribution for $sample..."
    
    singularity exec "$SIF" \
        /root/miniconda3/bin/LengthDistribution \
        -i "$bam_file" \
        -o "$output_prefix" \
        -f bam
    
    if [[ $? -eq 0 ]]; then
        echo "✓ Length distribution analysis completed for $sample"
    else
        echo "✗ Length distribution analysis failed for $sample"
    fi
done

# Step 4: Create summary report
echo ""
echo "Step 4: Creating summary report..."
echo "================================="

SUMMARY_FILE="${QC_OUTPUT_DIR}/quality_control_summary.txt"
cat > "$SUMMARY_FILE" << EOF
RiboMiner Quality Control Summary
Generated on: $(date)

Analysis Directory: $WORK_DIR
BAM Files Directory: $BAM_DIR
Output Directory: $QC_OUTPUT_DIR

Samples Analyzed:
EOF

for sample in "${!BAM_FILES[@]}"; do
    echo "- $sample: ${BAM_FILES[$sample]} (${BAM_LEGENDS[$sample]})" >> "$SUMMARY_FILE"
done

cat >> "$SUMMARY_FILE" << EOF

Files Generated:
1. Periodicity Analysis: *_periodicity.pdf
2. Frame Distribution: frame_distribution_*.pdf and frame_distribution_*.txt
3. Length Distribution: *_length_dist.pdf and *_length_dist.txt
4. Temporary Attributes File: temp_attributes.txt

Next Steps:
1. Review periodicity plots (*_periodicity.pdf) to determine optimal read lengths and offsets
2. Update the attributes.txt file with the correct parameters
3. Re-run frame distribution analysis if needed with updated parameters
4. Proceed with metagene analysis using the finalized attributes.txt

Important Notes:
- The temp_attributes.txt file contains estimated parameters
- Review periodicity results to optimize read lengths and offsets
- Good ribosome profiling data should show strong 3-nt periodicity
- Length distribution should peak around 28-30 nt
- Frame distribution should show enrichment in one reading frame
EOF

echo "Quality control analysis completed!"
echo "Summary report: $SUMMARY_FILE"
echo ""
echo "IMPORTANT: Please review the periodicity plots and update the attributes.txt file"
echo "with the optimal read lengths and offsets before proceeding to metagene analysis."
echo ""
echo "Output files are located in: $QC_OUTPUT_DIR"
