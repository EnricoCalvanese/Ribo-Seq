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
ORIGINAL_GTF="/global/scratch/users/enricocalvane/riboseq/Athaliana_447_Araport11.gene.gtf"
FIXED_GTF="/global/scratch/users/enricocalvane/riboseq/Athaliana_447_Araport11.gene.ribominer_fixed.gtf"
FASTA="/global/scratch/users/enricocalvane/riboseq/Athaliana_447_TAIR10.fa"
OUTDIR="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/prepared_transcripts"
SIF="ribocode_ribominer_latest.sif"

# Create output directory
mkdir -p "$OUTDIR"

# Step 1: Fix GTF file for RiboCode compatibility
echo "Step 1: Fixing GTF file for RiboCode compatibility..."
if [ ! -f "$FIXED_GTF" ]; then
    echo "Creating RiboCode-compatible GTF file..."
    awk '{
        # Remove trailing whitespace and semicolon if present
        gsub(/[[:space:]]*;?[[:space:]]*$/, "")
        
        # Extract and modify gene_id
        if (match($0, /gene_id "[^"]*"/)) {
            gene_id_full = substr($0, RSTART, RLENGTH)
            # Extract just the gene ID part (e.g., AT1G01010 from AT1G01010.Araport11.447)
            if (match(gene_id_full, /"([^.]+)/)) {
                base_gene_id = substr(gene_id_full, RSTART+1, RLENGTH-1)
                # Replace the full gene_id with the base gene_id
                gsub(/gene_id "[^"]*"/, "gene_id \"" base_gene_id "\"")
            }
        }
        
        # Add transcript_type if not present
        if ($0 !~ /transcript_type/) {
            print $0 "; transcript_type \"protein_coding\";"
        } else {
            print $0 ";"
        }
    }' "$ORIGINAL_GTF" > "$FIXED_GTF"
    
    echo "Fixed GTF file created: $FIXED_GTF"
    echo "Sample of fixed file:"
    head -3 "$FIXED_GTF"
    
    echo ""
    echo "Verification - gene_id format:"
    grep -m 1 "gene_id" "$FIXED_GTF" | sed 's/.*\(gene_id "[^"]*"\).*/\1/'
else
    echo "Fixed GTF file already exists: $FIXED_GTF"
fi

# Step 2: Run prepare_transcripts with the fixed GTF
echo ""
echo "Step 2: Running prepare_transcripts with Singularity..."
echo "Using fixed GTF: $FIXED_GTF"
echo "Using FASTA: $FASTA"
echo "Output directory: $OUTDIR"

singularity exec "$SIF" \
  /root/miniconda3/bin/prepare_transcripts \
  -g "$FIXED_GTF" \
  -f "$FASTA" \
  -o "$OUTDIR"

echo ""
echo "prepare_transcripts completed successfully."
echo "Output files in: $OUTDIR"
ls -la "$OUTDIR"

date
echo "Job completed on $(hostname)"
