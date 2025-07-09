#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=ribominer_comprehensive_qc.log

# RiboMiner Comprehensive Quality Control Script
# This script performs comprehensive quality control analysis on ribosome profiling data

# Define variables
WORK_DIR="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer"
BAM_DIR="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/transcriptome_aligned_reads/unique_reads"
GENOME_BAM_DIR="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/genome_aligned_reads/unique_reads"
RIBOCODE_ANNOT="${WORK_DIR}/prepared_transcripts"
LONGEST_TRANSCRIPTS_INFO="${WORK_DIR}/longest.transcripts.info.txt"
ATTRIBUTES_FILE="${WORK_DIR}/attributes_genome.txt"
MODIFIED_GTF="/global/scratch/users/enricocalvane/riboseq/Araport11_GTF_genes_ribominer.gtf"
OUTPUT_DIR="${WORK_DIR}/comprehensive_qc_analysis"

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

echo "Starting RiboMiner Comprehensive Quality Control Analysis"
echo "Work directory: $WORK_DIR"
echo "BAM directory: $BAM_DIR"
echo "Genome BAM directory: $GENOME_BAM_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Attributes file: $ATTRIBUTES_FILE"
echo "Modified GTF: $MODIFIED_GTF"
echo "========================================"

# Check if required files exist
if [[ ! -f "$LONGEST_TRANSCRIPTS_INFO" ]]; then
    echo "Error: longest.transcripts.info.txt not found at $LONGEST_TRANSCRIPTS_INFO"
    exit 1
fi

if [[ ! -f "$ATTRIBUTES_FILE" ]]; then
    echo "Error: attributes_genome.txt not found at $ATTRIBUTES_FILE"
    exit 1
fi

if [[ ! -f "$MODIFIED_GTF" ]]; then
    echo "Error: Modified GTF file not found at $MODIFIED_GTF"
    exit 1
fi

# 1. PERIODICITY ANALYSIS
echo ""
echo "1. PERIODICITY ANALYSIS"
echo "======================="
echo "Analyzing periodicity for each sample..."

PERIODICITY_DIR="${OUTPUT_DIR}/periodicity"
mkdir -p "$PERIODICITY_DIR"

for sample in "${!BAM_FILES[@]}"; do
    bam_file="${BAM_DIR}/${BAM_FILES[$sample]}"
    output_prefix="${PERIODICITY_DIR}/${sample}_periodicity"
    
    echo "Processing $sample (${BAM_LEGENDS[$sample]})..."
    
    if [[ ! -f "$bam_file" ]]; then
        echo "Warning: BAM file not found: $bam_file"
        continue
    fi
    
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

# 2. READING FRAME DISTRIBUTION ANALYSIS
echo ""
echo "2. READING FRAME DISTRIBUTION ANALYSIS"
echo "======================================"
echo "Analyzing reads distribution among different reading frames..."

FRAMES_DIR="${OUTPUT_DIR}/reading_frames"
mkdir -p "$FRAMES_DIR"

output_prefix="${FRAMES_DIR}/reading_frames_distribution"

RiboDensityOfDiffFrames \
    -f "$ATTRIBUTES_FILE" \
    -c "$LONGEST_TRANSCRIPTS_INFO" \
    -o "$output_prefix" \
    --plot yes

if [[ $? -eq 0 ]]; then
    echo "✓ Reading frame distribution analysis completed"
else
    echo "✗ Reading frame distribution analysis failed"
    exit 1
fi

# 3. LENGTH DISTRIBUTION ANALYSIS
echo ""
echo "3. LENGTH DISTRIBUTION ANALYSIS"
echo "==============================="
echo "Analyzing length distribution for each sample..."

LENGTH_DIR="${OUTPUT_DIR}/length_distribution"
mkdir -p "$LENGTH_DIR"

for sample in "${!BAM_FILES[@]}"; do
    bam_file="${BAM_DIR}/${BAM_FILES[$sample]}"
    output_prefix="${LENGTH_DIR}/${sample}_length_distribution"
    
    echo "Processing $sample (${BAM_LEGENDS[$sample]})..."
    
    if [[ ! -f "$bam_file" ]]; then
        echo "Warning: BAM file not found: $bam_file"
        continue
    fi
    
    LengthDistribution \
        -i "$bam_file" \
        -o "$output_prefix" \
        -f bam
    
    if [[ $? -eq 0 ]]; then
        echo "✓ Length distribution analysis completed for $sample"
    else
        echo "✗ Length distribution analysis failed for $sample"
        exit 1
    fi
done

# 4. READS LENGTH IN SPECIFIC REGIONS ANALYSIS
echo ""
echo "4. READS LENGTH IN SPECIFIC REGIONS ANALYSIS"
echo "============================================"
echo "Analyzing length distribution for reads in different regions..."

REGIONS_DIR="${OUTPUT_DIR}/reads_length_regions"
mkdir -p "$REGIONS_DIR"

# Define regions to analyze
REGIONS=("CDS" "5UTR" "3UTR")

for sample in "${!BAM_FILES[@]}"; do
    bam_file="${BAM_DIR}/${BAM_FILES[$sample]}"
    
    echo "Processing $sample (${BAM_LEGENDS[$sample]})..."
    
    if [[ ! -f "$bam_file" ]]; then
        echo "Warning: BAM file not found: $bam_file"
        continue
    fi
    
    for region in "${REGIONS[@]}"; do
        output_prefix="${REGIONS_DIR}/${sample}_${region}_length"
        
        echo "  Analyzing $region region..."
        
        ReadsLengthOfSpecificRegions \
            -i "$bam_file" \
            -c "$LONGEST_TRANSCRIPTS_INFO" \
            -o "$output_prefix" \
            --type "$region"
        
        if [[ $? -eq 0 ]]; then
            echo "    ✓ $region analysis completed for $sample"
        else
            echo "    ✗ $region analysis failed for $sample"
            exit 1
        fi
    done
done

# 5. DNA CONTAMINATION ANALYSIS
echo ""
echo "5. DNA CONTAMINATION ANALYSIS"
echo "============================"
echo "Analyzing DNA contamination for each sample..."

DNA_CONTAM_DIR="${OUTPUT_DIR}/dna_contamination"
mkdir -p "$DNA_CONTAM_DIR"

for sample in "${!BAM_FILES[@]}"; do
    genome_bam_file="${GENOME_BAM_DIR}/${BAM_FILES[$sample]}"
    output_prefix="${DNA_CONTAM_DIR}/${sample}_dna_contamination"
    
    echo "Processing $sample (${BAM_LEGENDS[$sample]})..."
    
    if [[ ! -f "$genome_bam_file" ]]; then
        echo "Warning: Genome BAM file not found: $genome_bam_file"
        continue
    fi
    
    StatisticReadsOnDNAsContam \
        -i "$genome_bam_file" \
        -g "$MODIFIED_GTF" \
        -o "$output_prefix"
    
    if [[ $? -eq 0 ]]; then
        echo "✓ DNA contamination analysis completed for $sample"
    else
        echo "✗ DNA contamination analysis failed for $sample"
        exit 1
    fi
done

# Create comprehensive summary report
echo ""
echo "Creating comprehensive summary report..."
echo "======================================="

SUMMARY_FILE="${OUTPUT_DIR}/comprehensive_qc_summary.txt"
cat > "$SUMMARY_FILE" << EOF
RiboMiner Comprehensive Quality Control Analysis Summary
=======================================================
Generated on: $(date)
Analysis Directory: $WORK_DIR
BAM Files Directory: $BAM_DIR
Genome BAM Files Directory: $GENOME_BAM_DIR
Output Directory: $OUTPUT_DIR

Configuration Files:
- Longest Transcripts Info: $LONGEST_TRANSCRIPTS_INFO
- Attributes File: $ATTRIBUTES_FILE
- Modified GTF: $MODIFIED_GTF

Samples Analyzed:
EOF

for sample in "${!BAM_FILES[@]}"; do
    echo "- $sample: ${BAM_FILES[$sample]} (${BAM_LEGENDS[$sample]})" >> "$SUMMARY_FILE"
done

cat >> "$SUMMARY_FILE" << EOF

Analyses Performed:
==================

1. Periodicity Analysis (periodicity/):
   - Assesses 3-nt periodicity of ribosome footprints
   - Files: *_periodicity.pdf, *_periodicity.txt
   - Purpose: Determine optimal read lengths and P-site offsets

2. Reading Frame Distribution (reading_frames/):
   - Analyzes read densities across different reading frames
   - Files: reading_frames_distribution.pdf, reading_frames_distribution.txt
   - Purpose: Assess translation fidelity and frame preference

3. Length Distribution (length_distribution/):
   - Analyzes length distribution of ribosome footprints
   - Files: *_length_distribution.pdf, *_length_distribution.txt
   - Purpose: Evaluate quality of ribosome footprint preparation

4. Reads Length in Specific Regions (reads_length_regions/):
   - Analyzes length distribution in CDS, 5'UTR, and 3'UTR regions
   - Files: *_CDS_length.pdf, *_5UTR_length.pdf, *_3UTR_length.pdf
   - Purpose: Compare footprint characteristics across transcript regions

5. DNA Contamination Analysis (dna_contamination/):
   - Assesses DNA contamination in ribosome profiling data
   - Files: *_dna_contamination.pdf, *_dna_contamination.txt
   - Purpose: Evaluate RNA-seq contamination and data quality

Quality Control Expectations:
============================
- Good periodicity: Strong 3-nt periodicity in CDS regions
- Frame distribution: Majority of reads in frame 0 (correct reading frame)
- Length distribution: Typical range 25-35 nt with clear modal lengths
- Region-specific lengths: CDS reads typically shorter than UTR reads
- DNA contamination: Low contamination indicates good rRNA depletion

Next Steps:
===========
1. Review all generated plots and statistics
2. Use periodicity results to optimize P-site offset parameters
3. Assess data quality based on frame distribution and contamination levels
4. Proceed with metagene analysis if quality metrics are satisfactory

Output Directory Structure:
===========================
$OUTPUT_DIR/
├── periodicity/                    # Periodicity analysis results
├── reading_frames/                 # Reading frame distribution
├── length_distribution/            # Length distribution analysis
├── reads_length_regions/           # Region-specific length analysis
├── dna_contamination/              # DNA contamination assessment
└── comprehensive_qc_summary.txt    # This summary file
EOF

echo ""
echo "Comprehensive quality control analysis completed!"
echo "================================================="
echo "Summary report: $SUMMARY_FILE"
echo "Output files are located in: $OUTPUT_DIR"
echo ""
echo "Analysis breakdown:"
echo "  - Periodicity analysis: $PERIODICITY_DIR"
echo "  - Reading frame distribution: $FRAMES_DIR"
echo "  - Length distribution: $LENGTH_DIR"
echo "  - Region-specific lengths: $REGIONS_DIR"
echo "  - DNA contamination: $DNA_CONTAM_DIR"
echo ""
echo "Please review all generated plots and statistics before proceeding."
