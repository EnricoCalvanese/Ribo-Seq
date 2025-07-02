#!/bin/bash

# Paths
ORIGINAL_GTF="/global/scratch/users/enricocalvane/riboseq/Araport11_GTF_genes_transposons.current.filtered.gtf"
MODIFIED_GTF="/global/scratch/users/enricocalvane/riboseq/Araport11_GTF_genes_ribominer.gtf"

# Let's examine the exact structure and look for potential parsing issues
echo "=== Checking for problematic lines ==="

echo "1. Lines that might have parsing issues (malformed attributes):"
awk 'NF < 9 {print "Line " NR ": " $0}' "$MODIFIED_GTF" | head -5

echo ""
echo "2. Checking attribute format consistency:"
echo "Lines where attributes don't end with semicolon:"
awk '$9 !~ /;$/ {print "Line " NR ": " $9}' "$MODIFIED_GTF" | head -5

echo ""
echo "3. Lines with unusual characters or formatting:"
awk '/[^[:print:]\t]/ {print "Line " NR ": contains non-printable chars"}' "$MODIFIED_GTF" | head -5

echo ""
echo "4. Check if any non-gene lines are missing transcript_id:"
awk '$3 != "gene" && $0 !~ /transcript_id/ {print "Line " NR ": " $0}' "$MODIFIED_GTF" | head -5

echo ""
echo "5. Let's see lines around where the error might occur (first 20 lines):"
head -20 "$MODIFIED_GTF" | nl

echo ""
echo "6. Let's try a minimal test - create a very simple GTF with just a few lines:"

# Create a minimal test GTF
TEST_GTF="/tmp/test_minimal.gtf"
cat > "$TEST_GTF" << 'EOF'
Chr1	Araport11	gene	3631	5899	.	+	.	gene_id "AT1G01010";
Chr1	Araport11	transcript	3631	5899	.	+	.	transcript_id "AT1G01010.1"; gene_id "AT1G01010"; transcript_type "protein_coding";
Chr1	Araport11	CDS	3760	3913	.	+	0	transcript_id "AT1G01010.1"; gene_id "AT1G01010"; transcript_type "protein_coding";
EOF

echo "Created minimal test GTF:"
cat "$TEST_GTF"

echo ""
echo "7. Test with the minimal GTF:"
singularity exec ribocode_ribominer_latest.sif \
  /root/miniconda3/bin/prepare_transcripts \
  -g "$TEST_GTF" \
  -f "/global/scratch/users/enricocalvane/riboseq/Athaliana_447_TAIR10.fa" \
  -o "/tmp/test_output"
