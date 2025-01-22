#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=0:10:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=create_gtf
#SBATCH --output=create_gtf_%j.out
#SBATCH --error=create_gtf_%j.err

# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/systemPipeR/attempt2

# Load required module
module load r

# Create R script
cat << 'EOF' > create_morf_gtf.R
# Load required libraries
library(rtracklayer)
library(GenomicRanges)

print("Starting GTF file creation with proper mORF definitions...")

# Load uORF predictions
print("Loading uORF predictions...")
uorf_data <- read.table("results/reference_analysis/final_matched_uorfs.txt",
                       header=TRUE, stringsAsFactors=FALSE)

# Get list of genes with uORFs
genes_with_uorfs <- unique(sub("_.*", "", uorf_data$uORF_ID))
print(paste("Found", length(genes_with_uorfs), "genes with uORFs"))

# Process uORFs (keeping previous correct implementation)
# ... [previous uORF processing code stays the same] ...

# New function to process mORFs correctly
create_complete_morf <- function(cds_entries) {
    # Sort CDS entries by position
    cds_entries <- cds_entries[order(cds_entries$start),]
    
    # Create single mORF entry
    data.frame(
        seqname = cds_entries$seqname[1],
        source = "systemPipeR",
        feature = "mORF",
        start = min(cds_entries$start),    # Most upstream start
        end = max(cds_entries$end),        # Most downstream end
        score = ".",
        strand = cds_entries$strand[1],
        frame = ".",
        attribute = cds_entries$attribute[1]
    )
}

# Process GTF for mORFs
print("Processing mORF entries...")
gtf$gene_id <- sub('.*gene_id "([^"]+)".*', "\\1", gtf$V9)
gtf$gene_id <- sub("gene:", "", gtf$gene_id)
gtf$transcript_id <- sub('.*transcript_id "([^"]+)".*', "\\1", gtf$V9)
gtf$transcript_id <- sub("transcript:", "", gtf$transcript_id)

# Filter for CDS entries of genes with uORFs
cds_entries <- gtf[gtf$V3 == "CDS" & gtf$gene_id %in% genes_with_uorfs,]

# Process each gene to create complete mORFs
morf_gtf <- data.frame()
for(gene in genes_with_uorfs) {
    # Get all CDS entries for this gene
    gene_cds <- cds_entries[cds_entries$gene_id == gene,]
    if(nrow(gene_cds) == 0) next
    
    # Find primary transcript (prefer .1)
    transcripts <- unique(gene_cds$transcript_id)
    primary_transcript <- transcripts[grep("\\.1$", transcripts)]
    if(length(primary_transcript) == 0) {
        primary_transcript <- transcripts[1]
    }
    
    # Get CDS entries for primary transcript
    transcript_cds <- gene_cds[gene_cds$transcript_id == primary_transcript,]
    
    # Create complete mORF entry
    morf_entry <- create_complete_morf(transcript_cds)
    morf_gtf <- rbind(morf_gtf, morf_entry)
}

print(paste("Created", nrow(morf_gtf), "complete mORF entries"))
print(paste("Average CDS entries per gene before combining:", 
            round(nrow(cds_entries)/length(genes_with_uorfs), 2)))
print(paste("mORFs after combining (should be ≤ genes):", nrow(morf_gtf)))

# Write GTF files
print("Writing GTF files...")
write.table(uorf_gtf, "results/counts/predicted_uorfs.gtf",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(morf_gtf, "results/counts/morfs.gtf",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

print("GTF file creation complete.")
print(paste("Created uORF GTF with", nrow(uorf_gtf), "entries"))
print(paste("Created mORF GTF with", nrow(morf_gtf), "entries"))

# Print example entries for verification
print("\nFirst few entries of mORF GTF:")
print(head(morf_gtf))
EOF

# Run R script
Rscript create_morf_gtf.R
