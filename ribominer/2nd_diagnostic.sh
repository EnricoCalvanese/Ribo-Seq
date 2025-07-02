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

# Step 1: Create a properly formatted GTF with exon features
echo "Step 1: Creating GTF with required exon features..."
if [ ! -f "$MODIFIED_GTF" ]; then
    echo "The issue is that RiboCode requires exon features, but your GTF only has CDS."
    echo "Creating modified GTF file with exon features derived from CDS regions..."
    
    awk '{
        # Remove trailing whitespace and semicolon if present
        gsub(/[[:space:]]*;?[[:space:]]*$/, "")
        
        if ($3 == "gene") {
            # For gene lines, keep only gene_id
            attributes = $9
            for (i = 10; i <= NF; i++) {
                attributes = attributes " " $i
            }
            
            # Extract gene_id from attributes
            if (match(attributes, /gene_id "[^"]*"/)) {
                gene_id_part = substr(attributes, RSTART, RLENGTH)
                print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" gene_id_part ";"
            }
        } else if ($3 == "transcript") {
            # For transcript lines, add transcript_type
            print $0 "; transcript_type \"protein_coding\";"
        } else if ($3 == "CDS") {
            # For each CDS, create a corresponding exon feature first, then the CDS
            # Extract transcript_id and gene_id for the exon line
            attributes = $9
            for (i = 10; i <= NF; i++) {
                attributes = attributes " " $i
            }
            
            # Create exon line with same coordinates as CDS
            exon_line = $1 "\t" $2 "\texon\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t.\t" attributes "; transcript_type \"protein_coding\";"
            print exon_line
            
            # Then print the CDS line
            print $0 "; transcript_type \"protein_coding\";"
        }
    }' "$ORIGINAL_GTF" > "$MODIFIED_GTF"
    
    echo "Modified GTF file created: $MODIFIED_GTF"
    echo "Sample of modified file:"
    head -10 "$MODIFIED_GTF"
    echo ""
    echo "Feature types in modified GTF:"
    awk '{print $3}' "$MODIFIED_GTF" | sort | uniq -c
else
    echo "Modified GTF file already exists: $MODIFIED_GTF"
fi

echo ""
echo "Testing with a corrected minimal GTF first..."
# Create a minimal test GTF with exon features
TEST_GTF="/tmp/test_with_exons.gtf"
cat > "$TEST_GTF" << 'EOF'
Chr1	Araport11	gene	3631	5899	.	+	.	gene_id "AT1G01010";
Chr1	Araport11	transcript	3631	5899	.	+	.	transcript_id "AT1G01010.1"; gene_id "AT1G01010"; transcript_type "protein_coding";
Chr1	Araport11	exon	3760	3913	.	+	.	transcript_id "AT1G01010.1"; gene_id "AT1G01010"; transcript_type "protein_coding";
Chr1	Araport11	CDS	3760	3913	.	+	0	transcript_id "AT1G01010.1"; gene_id "AT1G01010"; transcript_type "protein_coding";
EOF

echo "Testing minimal GTF with exon features:"
echo "Content:"
cat "$TEST_GTF"
echo ""
echo "Running test..."
singularity exec "$SIF" \
  /root/miniconda3/bin/prepare_transcripts \
  -g "$TEST_GTF" \
  -f "$FASTA" \
  -o "/tmp/test_output_with_exons"

if [ $? -eq 0 ]; then
    echo ""
    echo "Test successful! Now running with full modified GTF..."

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
