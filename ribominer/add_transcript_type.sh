#!/bin/bash
# Script to add transcript_type "protein_coding" to GTF file for RiboMiner compatibility

# Check if input file is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <input_gtf_file> [output_gtf_file]"
    echo "Example: $0 Athaliana_447_Araport11.gene.gtf Athaliana_447_Araport11.gene.modified.gtf"
    exit 1
fi

INPUT_GTF="$1"
OUTPUT_GTF="${2:-${INPUT_GTF%.gtf}.modified.gtf}"

# Check if input file exists
if [ ! -f "$INPUT_GTF" ]; then
    echo "Error: Input file '$INPUT_GTF' not found!"
    exit 1
fi

echo "Processing GTF file: $INPUT_GTF"
echo "Output file: $OUTPUT_GTF"

# Add transcript_type "protein_coding" to each line
# This adds the attribute before the final semicolon (if present) or at the end
awk '{
    # Remove trailing whitespace and semicolon if present
    gsub(/[[:space:]]*;?[[:space:]]*$/, "")
    
    # Add transcript_type attribute
    print $0 "; transcript_type \"protein_coding\";"
}' "$INPUT_GTF" > "$OUTPUT_GTF"

echo "Successfully created modified GTF file: $OUTPUT_GTF"
echo ""
echo "Sample of modified file:"
head -5 "$OUTPUT_GTF"
echo ""
echo "File is ready for RiboMiner!"
