#!/bin/bash

# Reorder GTF features so genes come first, then transcripts, then exons before CDS
# Usage: ./reorder_gtf.sh input.gtf output.gtf

INPUT_FILE="$1"
OUTPUT_FILE="$2"

if [ $# -ne 2 ]; then
    echo "Usage: $0 input.gtf output.gtf"
    exit 1
fi

# Sort by chromosome, gene_id, feature priority, start position
awk -F'\t' 'BEGIN{OFS="\t"} {
    # Extract gene_id (this is consistent across gene/transcript/exon/CDS)
    match($9, /gene_id "([^"]+)"/, gene_match)
    gene_id = gene_match[1]
    
    # Extract transcript_id for sub-sorting
    match($9, /transcript_id "([^"]+)"/, trans_match)
    transcript_id = trans_match[1]
    
    # Assign sort priority: gene=0, transcript=1, exon=2, CDS=3
    if ($3 == "gene") priority = 0
    else if ($3 == "transcript") priority = 1  
    else if ($3 == "exon") priority = 2
    else if ($3 == "CDS") priority = 3
    else priority = 4
    
    # Create sort key: chromosome:gene_id:priority:transcript_id:start_position
    sort_key = $1 ":" gene_id ":" priority ":" transcript_id ":" sprintf("%010d", $4)
    
    # Store line with sort key
    lines[sort_key] = $0
}
END {
    # Sort by keys and print
    n = asorti(lines, sorted_keys)
    for (i = 1; i <= n; i++) {
        print lines[sorted_keys[i]]
    }
}' "$INPUT_FILE" > "$OUTPUT_FILE"

echo "Reordering complete!"
echo "Original lines: $(wc -l < "$INPUT_FILE")"
echo "Reordered lines: $(wc -l < "$OUTPUT_FILE")"
