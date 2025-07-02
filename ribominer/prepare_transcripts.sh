#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=prepare_transcripts.log

set -euo pipefail

echo "Starting prepare_transcripts job on $(hostname)"
date

# Paths
ORIGINAL_GTF="/global/scratch/users/enricocalvane/riboseq/Araport11_GTF_genes_transposons.current.filtered.gtf"
MODIFIED_GTF="/global/scratch/users/enricocalvane/riboseq/Araport11_GTF_genes_ribominer.gtf"
FASTA="/global/scratch/users/enricocalvane/riboseq/Athaliana_447_TAIR10.fa"
OUTDIR="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/prepared_transcripts"
SIF="ribocode_ribominer_latest.sif"

# Create output directory
mkdir -p "$OUTDIR"

# Step 1: First, let's analyze the original GTF file to understand the issue
echo "Step 1: Analyzing the original GTF file..."
echo "Feature types in the original GTF:"
awk '{print $3}' "$ORIGINAL_GTF" | sort | uniq -c

echo ""
echo "Lines missing transcript_id (excluding gene lines):"
awk '$3 != "gene" && !/transcript_id/ {print NR ": " $0}' "$ORIGINAL_GTF" | head -10

echo ""
echo "Step 2: Creating modified GTF file..."
if [ ! -f "$MODIFIED_GTF" ]; then
    echo "Creating modified GTF file with transcript_type information..."
    awk '{
        # Remove trailing whitespace and semicolon if present
        gsub(/[[:space:]]*;?[[:space:]]*$/, "")
        
        if ($3 == "gene") {
            # For gene lines, keep only gene_id (remove transcript_id if present)
            attributes = $9
            for (i = 10; i <= NF; i++) {
                attributes = attributes " " $i
            }
            
            # Extract gene_id from attributes
            if (match(attributes, /gene_id "[^"]*"/)) {
                gene_id_part = substr(attributes, RSTART, RLENGTH)
                print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" gene_id_part ";"
            } else {
                print $0 ";"
            }
        } else {
            # For non-gene lines, check if transcript_id exists
            if (!/transcript_id/) {
                # Skip lines without transcript_id as they will cause errors
                print "# SKIPPED LINE WITHOUT transcript_id: " $0 > "/dev/stderr"
                next
            }
            # Add transcript_type to lines that have transcript_id
            print $0 "; transcript_type \"protein_coding\";"
        }
    }' "$ORIGINAL_GTF" > "$MODIFIED_GTF"
    
    echo "Modified GTF file created: $MODIFIED_GTF"
    echo "Sample of modified file:"
    head -5 "$MODIFIED_GTF"
else
    echo "Modified GTF file already exists: $MODIFIED_GTF"
fi

# Step 2: Run prepare_transcripts with the modified GTF
echo ""
echo "Step 2: Running prepare_transcripts with Singularity..."
echo "Using modified GTF: $MODIFIED_GTF"
echo "Using FASTA: $FASTA"
echo "Output directory: $OUTDIR"

singularity exec "$SIF" \
  /root/miniconda3/bin/prepare_transcripts \
  -g "$MODIFIED_GTF" \
  -f "$FASTA" \
  -o "$OUTDIR"
date
echo "Job completed on $(hostname)"
