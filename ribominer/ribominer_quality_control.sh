#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=ribominer_quality_control.log

# RiboMiner Quality Control Script
# This script performs comprehensive quality control on ribosome profiling data
# including periodicity checking, frame distribution, and length distribution

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

# Step 1: Periodicity Analysis for each sample with targeted fix
echo ""
echo "Step 1: Analyzing periodicity for each sample..."
echo "================================================"

for sample in "${!BAM_FILES[@]}"; do
    bam_file="${BAM_DIR}/${BAM_FILES[$sample]}"
    output_prefix="${QC_OUTPUT_DIR}/${sample}_periodicity"
    
    echo "Processing $sample (${BAM_LEGENDS[$sample]})..."
    echo "Applying targeted fix to line 39 of Periodicity.py..."
    
    # Use --writable-tmpfs to create a writable overlay and fix the specific line
    singularity exec --writable-tmpfs "$SIF" /bin/bash -c "
        # Show the problematic line before fixing
        echo 'Before fix - Line 39:'
        sed -n '39p' /root/miniconda3/lib/python3.7/site-packages/RiboMiner/Periodicity.py
        
        # Fix specifically line 39: change .key() to .keys()
        # Using a more precise sed command to target only the problematic line
        sed -i '39s/transcript_dict\.key()/transcript_dict.keys()/' /root/miniconda3/lib/python3.7/site-packages/RiboMiner/Periodicity.py
        
        # Verify the fix
        echo 'After fix - Line 39:'
        sed -n '39p' /root/miniconda3/lib/python3.7/site-packages/RiboMiner/Periodicity.py
        
        # Also check line 52 to make sure it's already correct
        echo 'Line 52 (should already be correct):'
        sed -n '52p' /root/miniconda3/lib/python3.7/site-packages/RiboMiner/Periodicity.py
        
        # Now run the analysis
        echo 'Running Periodicity analysis with fixed code...'
        /root/miniconda3/bin/Periodicity \\
            -i '$bam_file' \\
            -a '$RIBOCODE_ANNOT' \\
            -o '$output_prefix' \\
            -c '$LONGEST_TRANSCRIPTS_INFO' \\
            -L 25 \\
            -R 35
    "
    
    if [[ $? -eq 0 ]]; then
        echo "✓ Periodicity analysis completed for $sample"
    else
        echo "✗ Periodicity analysis failed for $sample"
        echo "Debugging: Let's check what went wrong..."
        
        # If it still fails, let's debug further
        singularity exec --writable-tmpfs "$SIF" /bin/bash -c "
            # Apply the fix again
            sed -i '39s/transcript_dict\.key()/transcript_dict.keys()/' /root/miniconda3/lib/python3.7/site-packages/RiboMiner/Periodicity.py
            
            # Check the transcript matching
            echo 'Debugging transcript matching...'
            python3 -c \"
import pysam
import sys
sys.path.insert(0, '/root/miniconda3/lib/python3.7/site-packages')
from RiboMiner.RiboMiner import load_transcripts_for_analysis

# Load transcript data
select_trans, transcript_dict = load_transcripts_for_analysis('$LONGEST_TRANSCRIPTS_INFO')
print(f'Loaded {len(select_trans)} transcripts from longest transcripts file')
print(f'Transcript dict has {len(transcript_dict)} entries')

# Check BAM file references
bamfile = pysam.AlignmentFile('$bam_file', 'rb')
bam_refs = set(bamfile.references)
print(f'BAM file has {len(bam_refs)} references')

# Check intersection
intersection = set(bam_refs).intersection(set(select_trans)).intersection(set(transcript_dict.keys()))
print(f'Intersection has {len(intersection)} transcripts')

if len(intersection) == 0:
    print('No matching transcripts found!')
    print('First 5 BAM references:', list(bam_refs)[:5])
    print('First 5 select_trans:', list(select_trans)[:5])
    print('First 5 transcript_dict keys:', list(transcript_dict.keys())[:5])
else:
    print('Found matching transcripts:', len(intersection))

bamfile.close()
\"
        "
    fi
done

