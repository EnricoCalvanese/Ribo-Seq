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

# Step 1: Modify GTF file to add transcript_type information ONLY to transcript-related lines
echo "Step 1: Adding transcript_type information to GTF file..."
if [ ! -f "$MODIFIED_GTF" ]; then
    echo "Creating modified GTF file with transcript_type information..."
    awk '{
        # Remove trailing whitespace and semicolon if present
        gsub(/[[:space:]]*;?[[:space:]]*$/, "")
        
        if ($3 == "gene") {
            # For gene lines, remove transcript_id if present and keep only gene_id
            # Split the attributes part and rebuild without transcript_id
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
            # For transcript, exon, CDS, etc. lines, add transcript_type
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
  
echo ""
echo "prepare_transcripts completed successfully."
echo "Output files in: $OUTDIR"
ls -la "$OUTDIR"

date
echo "Job completed on $(hostname)"
